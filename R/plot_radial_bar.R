#' @importFrom rlang .data
utils::globalVariables(c("SampleID", "Abundance", "taxon", "fill_col",
                         "x_pos", "y_start", "y_end"))

#' Radial (circular) stacked bar chart for microbiome composition
#'
#' Draws a polar coordinate stacked bar chart where each bar represents a
#' sample and segments represent taxa. A visually striking alternative to
#' the standard stacked bar chart.
#'
#' @param data A list with elements `otu_table`, `tax_table`, and `sample_data`,
#'   as in [demo_data].
#' @param level Taxonomic level to display. Default is `"Phylum"`.
#' @param top_n Number of top taxa to show, others merged into "Others".
#'   Default is `8`.
#' @param group_by Column name in `sample_data` to order samples. If `NULL`,
#'   uses original sample order.
#' @param palette A named character vector of colors. If `NULL`, uses built-in.
#' @param inner_radius Size of the inner hole (0-1). Default is `0.3`.
#'
#' @return A ggplot object.
#' @export
#'
#' @examples
#' data(demo_data)
#' plot_radial_bar(demo_data)
#' plot_radial_bar(demo_data, level = "Family", top_n = 10)
#' plot_radial_bar(demo_data, group_by = "Group")
plot_radial_bar <- function(data,
                            level        = "Phylum",
                            top_n        = 8,
                            group_by     = NULL,
                            palette      = NULL,
                            inner_radius = 0.3) {

  # 1. 提取数据
  otu  <- as.data.frame(data$otu_table)
  tax  <- data$tax_table
  meta <- data$sample_data

  # 2. 合并丰度与分类信息
  otu$SampleID <- rownames(otu)
  otu_long <- tidyr::pivot_longer(
    otu, cols = -SampleID,
    names_to = "OTU", values_to = "Abundance"
  )
  otu_long <- merge(otu_long, tax[, c("OTU", level)], by = "OTU")

  # 3. 按分类层级聚合
  agg <- stats::aggregate(Abundance ~ SampleID + get(level),
                          data = otu_long, FUN = sum)
  colnames(agg)[2] <- level

  # 4. 找 top_n 分类
  mean_abu <- stats::aggregate(Abundance ~ get(level), data = agg, FUN = mean)
  colnames(mean_abu)[1] <- level
  mean_abu <- mean_abu[order(mean_abu$Abundance, decreasing = TRUE), ]
  top_taxa <- mean_abu[[level]][1:min(top_n, nrow(mean_abu))]

  agg[[level]] <- ifelse(agg[[level]] %in% top_taxa, agg[[level]], "Others")
  agg <- stats::aggregate(Abundance ~ SampleID + get(level),
                          data = agg, FUN = sum)
  colnames(agg)[2] <- level

  # 5. 归一化为相对丰度
  total <- stats::aggregate(Abundance ~ SampleID, data = agg, FUN = sum)
  agg   <- merge(agg, total, by = "SampleID", suffixes = c("", "_total"))
  agg$Abundance <- agg$Abundance / agg$Abundance_total

  # 6. 排序：Others 最后，其余按均值丰度
  taxa_order <- c(top_taxa[top_taxa %in% unique(agg[[level]])], "Others")
  agg[[level]] <- factor(agg[[level]], levels = taxa_order)

  # 7. 样本排序
  if (!is.null(group_by) && group_by %in% colnames(meta)) {
    meta_ord  <- meta[order(meta[[group_by]]), ]
    samp_order <- meta_ord$SampleID
  } else {
    samp_order <- unique(agg$SampleID)
  }
  samp_order <- samp_order[samp_order %in% agg$SampleID]
  agg$SampleID <- factor(agg$SampleID, levels = samp_order)

  # 8. 计算每个样本内的累积丰度（用于极坐标堆叠）
  agg <- agg[order(agg$SampleID, agg[[level]]), ]
  agg$y_end   <- unlist(lapply(split(agg$Abundance, agg$SampleID), cumsum))
  agg$y_start <- agg$y_end - agg$Abundance

  # 9. 配色
  all_taxa <- levels(agg[[level]])
  default_pal <- c(
    "#4E79A7","#F28E2B","#E15759","#76B7B2","#59A14F",
    "#EDC948","#B07AA1","#FF9DA7","#9C755F","#BAB0AC",
    "#D3D3D3"
  )
  if (is.null(palette)) {
    palette <- stats::setNames(
      rep(default_pal, length.out = length(all_taxa)),
      all_taxa
    )
  }
  if (!"Others" %in% names(palette)) palette["Others"] <- "#D3D3D3"

  # 10. x 轴位置（每个样本一个位置）
  n_samp   <- length(samp_order)
  agg$x_num <- as.numeric(agg$SampleID)

  # 11. 画图
  p <- ggplot2::ggplot(agg) +
    ggplot2::geom_rect(
      ggplot2::aes(
        xmin = x_num - 0.4,
        xmax = x_num + 0.4,
        ymin = y_start + inner_radius,
        ymax = y_end   + inner_radius,
        fill = .data[[level]]
      )
    ) +
    ggplot2::scale_fill_manual(
      name   = level,
      values = palette,
      guide  = ggplot2::guide_legend(ncol = 1)
    ) +
    # 样本标签
    # 样本标签放到每个扇形外侧
    ggplot2::geom_text(
      data = data.frame(
        x     = seq_len(n_samp),
        y     = 1 + inner_radius + 0.08,
        label = samp_order,
        stringsAsFactors = FALSE
      ),
      ggplot2::aes(x = x, y = y, label = label),
      size  = 3,
      color = "grey30",
      inherit.aes = FALSE
    ) +
    ggplot2::coord_polar(start = 0) +
    ggplot2::ylim(0, 1 + inner_radius + 0.25) +
    ggplot2::theme_void(base_size = 11) +
    ggplot2::theme(
      legend.position = "right",
      legend.title    = ggplot2::element_text(face = "bold", size = 10),
      legend.text     = ggplot2::element_text(size = 9)
    )

  return(p)
}
