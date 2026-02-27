#' Circular Voronoi Treemap for microbiome composition (interactive widget)
#'
#' Draws an interactive circular Voronoi treemap using the voronoiTreemap
#' package. Tile area represents relative abundance and color represents
#' the parent phylum. Output is an HTML widget viewable in RStudio Viewer
#' or a web browser.
#'
#' @param data A list with elements `otu_table`, `tax_table`, and `sample_data`,
#'   as in [demo_data].
#' @param sample Name of the sample to display. If `NULL`, uses the mean
#'   abundance across all samples.
#' @param level Taxonomic level for inner tiles. Default is `"Family"`.
#' @param color_by Taxonomic level for coloring. Default is `"Phylum"`.
#' @param top_n Number of top taxa to show. Default is `30`.
#' @param min_abundance Minimum abundance to include. Default is `0.005`.
#' @param palette A named character vector of colors. If `NULL`, uses built-in.
#' @param legend Logical. Whether to show the built-in legend. Default `TRUE`.
#'
#' @return An HTML widget object.
#' @export
#'
#' @examples
#' data(demo_data)
#' plot_voronoi_widget(demo_data)
#' plot_voronoi_widget(demo_data, level = "Genus", top_n = 20)
plot_voronoi_widget <- function(data,
                                sample        = NULL,
                                level         = "Family",
                                color_by      = "Phylum",
                                top_n         = 30,
                                min_abundance = 0.005,
                                palette       = NULL,
                                legend        = TRUE) {

  if (!requireNamespace("voronoiTreemap", quietly = TRUE)) {
    stop("Package 'voronoiTreemap' is required. ",
         "Install with: install.packages('voronoiTreemap')")
  }

  # 1. 提取数据
  otu <- as.data.frame(data$otu_table)
  tax <- data$tax_table

  # 2. 选取样本或取均值
  if (!is.null(sample)) {
    if (!sample %in% rownames(otu)) stop("Sample not found in otu_table.")
    abu <- otu[sample, , drop = FALSE]
  } else {
    abu <- as.data.frame(t(colMeans(otu)))
  }

  # 3. 合并丰度与分类
  abu_vec <- unlist(abu)
  df <- data.frame(
    OTU       = names(abu_vec),
    Abundance = as.numeric(abu_vec),
    stringsAsFactors = FALSE
  )
  df <- merge(df, tax[, c("OTU", level, color_by)], by = "OTU")

  # 4. 聚合到目标层级
  agg <- stats::aggregate(Abundance ~ get(level) + get(color_by),
                          data = df, FUN = sum)
  colnames(agg)[1:2] <- c(level, color_by)
  agg <- agg[agg$Abundance >= min_abundance, ]
  if (nrow(agg) == 0) stop("No taxa above min_abundance threshold.")
  agg <- agg[order(agg$Abundance, decreasing = TRUE), ]
  agg <- head(agg, top_n)

  # 5. 配色
  all_groups <- unique(agg[[color_by]])
  default_pal <- c(
    "#4E79A7","#F28E2B","#E15759","#76B7B2","#59A14F",
    "#EDC948","#B07AA1","#FF9DA7","#9C755F","#BAB0AC"
  )
  if (is.null(palette)) {
    palette <- stats::setNames(
      rep(default_pal, length.out = length(all_groups)),
      all_groups
    )
  }
  agg$color <- unname(palette[match(agg[[color_by]], names(palette))])

  # 6. 构建 voronoiTreemap 输入格式
  agg$weight <- agg$Abundance / sum(agg$Abundance) * 100
  vt_df <- data.frame(
    h1     = "Root",
    h2     = agg[[color_by]],
    h3     = agg[[level]],
    color  = agg$color,
    weight = agg$weight,
    codes  = agg[[level]],
    stringsAsFactors = FALSE
  )

  # 7. 构建树结构并渲染
  vt_tree <- voronoiTreemap::vt_input_from_df(vt_df)

  # 添加标签（截断过长的名称）
  add_labels <- function(node) {
    if (!is.null(node$children)) {
      node$children <- lapply(node$children, add_labels)
    } else {
      nm <- node$name
      node$label <- ifelse(nchar(nm) > 15,
                           paste0(substr(nm, 1, 13), "..."),
                           nm)
    }
    return(node)
  }
  vt_tree <- add_labels(vt_tree)

  # 8. 输出 widget
  vt_json   <- voronoiTreemap::vt_export_json(vt_tree)
  widget    <- voronoiTreemap::vt_d3(vt_json, legend = legend)

  return(widget)
}
