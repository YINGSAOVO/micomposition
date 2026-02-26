#' Stacked bar chart for microbiome composition
#'
#' @param data A list with elements `otu_table`, `tax_table`, and `sample_data`,
#'   as in [demo_data].
#' @param level Taxonomic level to display. One of `"Phylum"`, `"Class"`,
#'   `"Order"`, `"Family"`, or `"Genus"`. Default is `"Phylum"`.
#' @param top_n Number of top taxa to show. Others are merged into "Others".
#'   Default is `8`.
#' @param group_by Column name in `sample_data` to group samples on the x-axis.
#'   Default is `NULL` (show individual samples).
#' @param palette A character vector of colors. If `NULL`, uses a built-in palette.
#'
#' @return A ggplot object.
#' @export
#'
#' @examples
#' data(demo_data)
#' plot_bar_stack(demo_data)
#' plot_bar_stack(demo_data, level = "Genus", top_n = 10)
#' plot_bar_stack(demo_data, group_by = "Group")
plot_bar_stack <- function(data,
                           level    = "Phylum",
                           top_n    = 8,
                           group_by = NULL,
                           palette  = NULL) {

  # 1. 提取三张表
  otu  <- as.data.frame(data$otu_table)
  tax  <- data$tax_table
  meta <- data$sample_data

  # 2. 把 OTU 丰度表 合并上分类信息
  otu$SampleID <- rownames(otu)
  otu_long <- tidyr::pivot_longer(
    otu,
    cols      = -SampleID,
    names_to  = "OTU",
    values_to = "Abundance"
  )
  otu_long <- merge(otu_long, tax[, c("OTU", level)], by = "OTU")

  # 3. 按分类水平加总丰度
  otu_agg <- stats::aggregate(
    Abundance ~ SampleID + get(level),
    data = otu_long,
    FUN  = sum
  )
  colnames(otu_agg)[2] <- level

  # 4. 找出 top_n 个丰度最高的分类，其余归为 Others
  mean_abu <- stats::aggregate(
    Abundance ~ get(level),
    data = otu_agg,
    FUN  = mean
  )
  colnames(mean_abu)[1] <- level
  mean_abu <- mean_abu[order(mean_abu$Abundance, decreasing = TRUE), ]
  top_taxa <- mean_abu[[level]][1:min(top_n, nrow(mean_abu))]

  otu_agg[[level]] <- ifelse(
    otu_agg[[level]] %in% top_taxa,
    otu_agg[[level]],
    "Others"
  )

  otu_agg <- stats::aggregate(
    Abundance ~ SampleID + get(level),
    data = otu_agg,
    FUN  = sum
  )
  colnames(otu_agg)[2] <- level

  # 5. 如果指定了 group_by，合并元数据并按组平均
  if (!is.null(group_by)) {
    otu_agg <- merge(otu_agg, meta[, c("SampleID", group_by)], by = "SampleID")
    otu_agg <- stats::aggregate(
      Abundance ~ get(group_by) + get(level),
      data = otu_agg,
      FUN  = mean
    )
    colnames(otu_agg)[1:2] <- c(group_by, level)
    x_var <- group_by
  } else {
    x_var <- "SampleID"
  }

  # 6. 配色
  n_colors <- length(unique(otu_agg[[level]]))
  if (is.null(palette)) {
    palette <- c(
      "#4E79A7","#F28E2B","#E15759","#76B7B2","#59A14F",
      "#EDC948","#B07AA1","#FF9DA7","#9C755F","#BAB0AC",
      "#D3D3D3"
    )
  }
  # Others 始终用灰色
  taxa_levels <- c(top_taxa[top_taxa %in% unique(otu_agg[[level]])], "Others")
  otu_agg[[level]] <- factor(otu_agg[[level]], levels = taxa_levels)
  names(palette) <- taxa_levels

  # 7. 画图
  p <- ggplot2::ggplot(
    otu_agg,
    ggplot2::aes(
      x    = .data[[x_var]],
      y    = Abundance,
      fill = .data[[level]]
    )
  ) +
    ggplot2::geom_bar(stat = "identity", width = 0.7) +
    ggplot2::scale_fill_manual(values = palette) +
    ggplot2::scale_y_continuous(
      labels = scales::percent_format(accuracy = 1),
      expand = c(0, 0)
    ) +
    ggplot2::labs(
      x    = NULL,
      y    = "Relative Abundance",
      fill = level
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      axis.text.x  = ggplot2::element_text(angle = 45, hjust = 1),
      legend.title = ggplot2::element_text(face = "bold")
    )

  return(p)
}
