#' Static circular packed bubble chart for microbiome composition
#'
#' Draws a static circle packing chart where circle area represents relative
#' abundance and color represents the parent phylum. Suitable for publication.
#'
#' @param data A list with elements `otu_table`, `tax_table`, and `sample_data`,
#'   as in [demo_data].
#' @param sample Name of the sample to display. If `NULL`, uses the mean
#'   abundance across all samples.
#' @param level Taxonomic level to display. Default is `"Family"`.
#' @param color_by Taxonomic level for coloring. Default is `"Phylum"`.
#' @param top_n Number of top taxa to show. Default is `30`.
#' @param min_abundance Minimum abundance to include. Default is `0.005`.
#' @param palette A named character vector of colors. If `NULL`, uses built-in.
#' @param label_min_pct Minimum percentage (0-100) to show a label. Default `1.5`.
#' @param alpha Fill transparency, 0-1. Default `0.85`.
#'
#' @return A ggplot object.
#' @export
#'
#' @examples
#' data(demo_data)
#' plot_voronoi_static(demo_data)
#' plot_voronoi_static(demo_data, level = "Genus", top_n = 20)
plot_voronoi_static <- function(data,
                                sample        = NULL,
                                level         = "Family",
                                color_by      = "Phylum",
                                top_n         = 30,
                                min_abundance = 0.005,
                                palette       = NULL,
                                label_min_pct = 1.5,
                                alpha         = 0.85) {

  if (!requireNamespace("packcircles", quietly = TRUE)) {
    stop("Package 'packcircles' is required. ",
         "Install with: install.packages('packcircles')")
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

  # 4. 聚合
  agg <- stats::aggregate(Abundance ~ get(level) + get(color_by),
                          data = df, FUN = sum)
  colnames(agg)[1:2] <- c(level, color_by)
  agg <- agg[agg$Abundance >= min_abundance, ]
  if (nrow(agg) == 0) stop("No taxa above min_abundance threshold.")
  agg <- agg[order(agg$Abundance, decreasing = TRUE), ]
  agg <- head(agg, top_n)
  agg$pct <- agg$Abundance / sum(agg$Abundance) * 100

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
  agg$fill_color <- unname(palette[match(agg[[color_by]], names(palette))])

  # 6. 用 packcircles 计算圆的位置和半径
  #    面积正比于丰度 → 半径正比于 sqrt(丰度)
  packing <- packcircles::circleProgressiveLayout(
    agg$pct,
    sizetype = "area"
  )
  agg$x      <- packing$x
  agg$y      <- packing$y
  agg$radius <- packing$radius

  # 7. 生成每个圆的多边形坐标（用于 geom_polygon）
  circle_df <- packcircles::circleLayoutVertices(packing, npoints = 100)
  circle_df$fill_color <- rep(agg$fill_color, each = 101)
  circle_df$taxon      <- rep(agg[[level]],   each = 101)
  circle_df$group_var  <- rep(agg[[color_by]], each = 101)

  # 8. 标签数据
  label_df <- agg[agg$pct >= label_min_pct, ]
  label_df$short_label <- ifelse(
    nchar(label_df[[level]]) > 12,
    paste0(substr(label_df[[level]], 1, 11), "..."),
    label_df[[level]]
  )
  label_df$pct_label <- paste0(label_df$short_label, "\n",
                               round(label_df$pct, 1), "%")

  # 9. 画图
  # 给 circle_df 加上 color_by 列用于图例
  circle_df$phylum_group <- rep(agg[[color_by]], each = 101)

  p <- ggplot2::ggplot() +
    ggplot2::geom_polygon(
      data = circle_df,
      ggplot2::aes(x = x, y = y, group = id),
      fill  = rep(agg$fill_color, each = 101),
      color = "white",
      linewidth = 0.4,
      alpha = alpha
    ) +
    ggrepel::geom_text_repel(
      data = label_df,
      ggplot2::aes(x = x, y = y, label = pct_label),
      size         = 2.8,
      color        = "white",
      fontface     = "bold",
      bg.color     = "transparent",
      max.overlaps = Inf,
      seed         = 42
    ) +
    ggplot2::geom_point(
      data = data.frame(
        x   = rep(-9999, length(all_groups)),
        y   = rep(-9999, length(all_groups)),
        grp = all_groups,
        stringsAsFactors = FALSE
      ),
      ggplot2::aes(x = x, y = y, color = grp),
      size = 0
    ) +
    ggplot2::scale_color_manual(
      name   = color_by,
      values = palette,
      guide  = ggplot2::guide_legend(
        override.aes = list(size = 5, shape = 15),
        ncol = 1
      )
    ) +
    ggplot2::coord_fixed(
      xlim = c(min(circle_df$x) - 1, max(circle_df$x) + 1),
      ylim = c(min(circle_df$y) - 1, max(circle_df$y) + 1)
    ) +
    ggplot2::theme_void(base_size = 11) +
    ggplot2::theme(
      legend.position = "right",
      legend.title    = ggplot2::element_text(face = "bold", size = 10),
      legend.text     = ggplot2::element_text(size = 9)
    )

  return(p)
}
