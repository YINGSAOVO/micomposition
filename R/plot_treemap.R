#' @importFrom rlang .data
utils::globalVariables(c("area", "fill_col", "subgroup", "label_text", "pct"))

#' Rectangular treemap for microbiome composition
#'
#' Draws a rectangular treemap where tile area represents relative abundance.
#' Tiles are grouped and colored by a parent taxonomic level (e.g. Phylum),
#' with subdivisions showing a child level (e.g. Family).
#'
#' @param data A list with elements `otu_table`, `tax_table`, and `sample_data`,
#'   as in [demo_data].
#' @param sample Name of the sample to display. If `NULL`, uses the mean
#'   abundance across all samples.
#' @param level Taxonomic level for the inner tiles. Default is `"Family"`.
#' @param group_by Taxonomic level for grouping and coloring. Default is
#'   `"Phylum"`.
#' @param top_n Number of top taxa (at `level`) to display. Default is `30`.
#' @param min_abundance Minimum abundance to include a taxon. Default is `0.005`.
#' @param palette A named character vector of colors. If `NULL`, uses built-in.
#' @param label_min_pct Minimum percentage to show a tile label. Default `1`.
#' @param show_group_label Logical. Whether to show group (Phylum) labels.
#'   Default `TRUE`.
#'
#' @return A ggplot object.
#' @export
#'
#' @examples
#' data(demo_data)
#' plot_treemap(demo_data)
#' plot_treemap(demo_data, level = "Genus", top_n = 20)
#' plot_treemap(demo_data, group_by = "Phylum", level = "Family")
plot_treemap <- function(data,
                         sample           = NULL,
                         level            = "Family",
                         group_by         = "Phylum",
                         top_n            = 30,
                         min_abundance    = 0.005,
                         palette          = NULL,
                         label_min_pct    = 1,
                         show_group_label = TRUE) {

  if (!requireNamespace("treemapify", quietly = TRUE)) {
    stop("Package 'treemapify' is required. ",
         "Install with: install.packages('treemapify')")
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
  df <- merge(df, tax[, c("OTU", level, group_by)], by = "OTU")

  # 4. 聚合
  agg <- stats::aggregate(Abundance ~ get(level) + get(group_by),
                          data = df, FUN = sum)
  colnames(agg)[1:2] <- c(level, group_by)
  agg <- agg[agg$Abundance >= min_abundance, ]
  if (nrow(agg) == 0) stop("No taxa above min_abundance threshold.")
  agg <- agg[order(agg$Abundance, decreasing = TRUE), ]
  agg <- head(agg, top_n)
  agg$pct <- agg$Abundance / sum(agg$Abundance) * 100

  # 5. 配色
  all_groups <- unique(agg[[group_by]])
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
  agg$fill_col <- unname(palette[match(agg[[group_by]], names(palette))])

  # 6. 标签：只显示足够大的格子
  agg$label_text <- ifelse(
    agg$pct >= label_min_pct,
    paste0(agg[[level]], "\n", round(agg$pct, 1), "%"),
    ""
  )

  # 7. 画图
  p <- ggplot2::ggplot(
    agg,
    ggplot2::aes(
      area     = pct,
      fill     = fill_col,
      subgroup = .data[[group_by]],
      label    = label_text
    )
  ) +
    treemapify::geom_treemap(
      color     = "white",
      linewidth = 0.8
    ) +
    treemapify::geom_treemap_subgroup_border(
      color     = "white",
      linewidth = 2.5
    ) +
    ggplot2::scale_fill_identity() +
    treemapify::geom_treemap_text(
      color      = "white",
      fontface   = "bold",
      size       = 10,
      place      = "centre",
      grow       = FALSE,
      reflow     = TRUE,
      min.size   = 6
    )

  # 8. 加 Phylum 标签（左上角）
  if (show_group_label) {
    p <- p +
      treemapify::geom_treemap_subgroup_text(
        color    = "white",
        fontface = "italic",
        size     = 13,
        place    = "topleft",
        alpha    = 0.6,
        grow     = FALSE
      )
  }

  # 9. 图例
  p <- p +
    ggplot2::geom_point(
      data = data.frame(
        x   = rep(-9999, length(all_groups)),
        y   = rep(-9999, length(all_groups)),
        grp = all_groups,
        stringsAsFactors = FALSE
      ),
      ggplot2::aes(x = x, y = y, color = grp),
      inherit.aes = FALSE,
      size = 0
    ) +
    ggplot2::scale_color_manual(
      name   = group_by,
      values = palette,
      guide  = ggplot2::guide_legend(
        override.aes = list(size = 5, shape = 15),
        ncol = 1
      )
    ) +
    ggplot2::theme_void(base_size = 11) +
    ggplot2::theme(
      legend.position = "right",
      legend.title    = ggplot2::element_text(face = "bold", size = 10),
      legend.text     = ggplot2::element_text(size = 9)
    )

  return(p)
}
