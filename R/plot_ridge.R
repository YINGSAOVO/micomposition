#' @importFrom rlang .data
utils::globalVariables(c("Abundance", "SampleID", "taxon"))

#' Ridge plot for microbiome composition across samples
#'
#' Draws a ridge (joy) plot showing the abundance distribution of top taxa
#' across multiple samples. Each ridge represents one taxon, colored by
#' parent phylum.
#'
#' @param data A list with elements `otu_table`, `tax_table`, and `sample_data`,
#'   as in [demo_data].
#' @param level Taxonomic level to display. Default is `"Family"`.
#' @param color_by Taxonomic level for coloring ridges. Default is `"Phylum"`.
#' @param top_n Number of top taxa to display. Default is `12`.
#' @param group_by Column in `sample_data` to color sample points. If `NULL`,
#'   all points are the same color.
#' @param palette A named character vector of colors for `color_by`. If `NULL`,
#'   uses built-in palette.
#' @param scale Height scaling of ridges (overlap). Default is `2`.
#' @param alpha Ridge fill transparency. Default is `0.8`.
#'
#' @return A ggplot object.
#' @export
#'
#' @examples
#' data(demo_data)
#' plot_ridge(demo_data)
#' plot_ridge(demo_data, level = "Genus", top_n = 10)
#' plot_ridge(demo_data, group_by = "Group")
plot_ridge <- function(data,
                       level    = "Family",
                       color_by = "Phylum",
                       top_n    = 12,
                       group_by = NULL,
                       palette  = NULL,
                       scale    = 2,
                       alpha    = 0.8) {

  if (!requireNamespace("ggridges", quietly = TRUE)) {
    stop("Package 'ggridges' is required. ",
         "Install with: install.packages('ggridges')")
  }

  # 1. 提取数据
  otu  <- as.data.frame(data$otu_table)
  tax  <- data$tax_table
  meta <- data$sample_data

  # 2. 合并丰度与分类
  otu$SampleID <- rownames(otu)
  otu_long <- tidyr::pivot_longer(
    otu, cols = -SampleID,
    names_to = "OTU", values_to = "Abundance"
  )
  otu_long <- merge(otu_long, tax[, c("OTU", level, color_by)], by = "OTU")

  # 3. 按目标层级聚合
  agg <- stats::aggregate(Abundance ~ SampleID + get(level) + get(color_by),
                          data = otu_long, FUN = sum)
  colnames(agg)[2:3] <- c(level, color_by)

  # 4. 找 top_n（按全局均值丰度）
  mean_abu <- stats::aggregate(Abundance ~ get(level), data = agg, FUN = mean)
  colnames(mean_abu)[1] <- level
  mean_abu <- mean_abu[order(mean_abu$Abundance, decreasing = TRUE), ]
  top_taxa <- mean_abu[[level]][1:min(top_n, nrow(mean_abu))]
  agg <- agg[agg[[level]] %in% top_taxa, ]

  # 5. 合并样本元数据
  if (!is.null(group_by) && group_by %in% colnames(meta)) {
    agg <- merge(agg, meta[, c("SampleID", group_by)], by = "SampleID")
  }

  # 6. 排序：按均值丰度从高到低（ridge 图从下到上）
  taxa_order <- rev(top_taxa[top_taxa %in% unique(agg[[level]])])
  agg[[level]] <- factor(agg[[level]], levels = taxa_order)

  # 7. 配色（按 color_by）
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

  # 8. 画图
  p <- ggplot2::ggplot(
    agg,
    ggplot2::aes(
      x      = Abundance,
      y      = .data[[level]],
      fill   = .data[[color_by]]
    )
  ) +
    ggridges::geom_density_ridges2(
      scale          = scale,
      alpha          = alpha,
      color          = "white",
      linewidth      = 0.4,
      rel_min_height = 0.01
    ) +
    ggplot2::scale_fill_manual(
      name   = color_by,
      values = palette
    ) +
    ggplot2::scale_x_continuous(
      labels = scales::percent_format(accuracy = 1),
      expand = c(0, 0)
    ) +
    ggplot2::labs(
      x = "Relative Abundance",
      y = NULL
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      axis.text.y     = ggplot2::element_text(size = 10),
      axis.text.x     = ggplot2::element_text(size = 9),
      legend.position = "right",
      legend.title    = ggplot2::element_text(face = "bold", size = 10),
      legend.text     = ggplot2::element_text(size = 9),
      panel.grid.major.x = ggplot2::element_line(color = "grey90")
    )

  # 9. 如果有 group_by，加样本点
  if (!is.null(group_by) && group_by %in% colnames(agg)) {
    p <- p +
      ggplot2::geom_rug(
        ggplot2::aes(color = .data[[group_by]]),
        sides = "b", alpha = 0.6,
        length = ggplot2::unit(0.03, "npc")
      )
  }

  return(p)
}
