#' @importFrom rlang .data
utils::globalVariables(c("col", "row", "fill_col", "taxon"))

#' Waffle chart for microbiome composition
#'
#' Draws a waffle (square pie) chart where each small square represents a
#' fixed percentage of the community. Color represents taxonomic group.
#'
#' @param data A list with elements `otu_table`, `tax_table`, and `sample_data`,
#'   as in [demo_data].
#' @param sample Name of the sample to display. If `NULL`, uses the mean
#'   abundance across all samples.
#' @param level Taxonomic level to display. Default is `"Phylum"`.
#' @param top_n Number of top taxa to show, others merged into "Others".
#'   Default is `8`.
#' @param n_rows Number of rows in the waffle grid. Default is `10`.
#' @param palette A named character vector of colors. If `NULL`, uses built-in.
#' @param legend_title Title for the legend. Default uses `level`.
#'
#' @return A ggplot object.
#' @export
#'
#' @examples
#' data(demo_data)
#' plot_waffle(demo_data)
#' plot_waffle(demo_data, level = "Family", top_n = 10, n_rows = 10)
plot_waffle <- function(data,
                        sample       = NULL,
                        level        = "Phylum",
                        top_n        = 8,
                        n_rows       = 10,
                        palette      = NULL,
                        legend_title = NULL) {

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

  # 3. 合并丰度
  abu_vec <- unlist(abu)
  df <- data.frame(
    OTU       = names(abu_vec),
    Abundance = as.numeric(abu_vec),
    stringsAsFactors = FALSE
  )
  df <- merge(df, tax[, c("OTU", level)], by = "OTU")

  # 4. 聚合 + 归一化
  agg <- stats::aggregate(Abundance ~ get(level), data = df, FUN = sum)
  colnames(agg)[1] <- level
  agg <- agg[order(agg$Abundance, decreasing = TRUE), ]

  # 5. 前 top_n 保留，其余归 Others
  if (nrow(agg) > top_n) {
    top    <- agg[1:top_n, ]
    others <- data.frame(setNames(list("Others"), level),
                         Abundance = sum(agg$Abundance[(top_n+1):nrow(agg)]))
    agg <- rbind(top, others)
  }
  agg$pct <- agg$Abundance / sum(agg$Abundance) * 100

  # 6. 四舍五入到整数格（总共 100 格）
  agg$n_squares <- round(agg$pct)
  # 修正总和为 100
  diff <- 100 - sum(agg$n_squares)
  if (diff != 0) {
    # 把差值加到最大的那个
    idx <- which.max(agg$n_squares)
    agg$n_squares[idx] <- agg$n_squares[idx] + diff
  }
  agg <- agg[agg$n_squares > 0, ]

  # 7. 配色
  all_taxa <- agg[[level]]
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

  # 8. 生成格子坐标
  taxon_vec <- rep(agg[[level]], times = agg$n_squares)
  n_cols    <- ceiling(100 / n_rows)
  grid_df   <- data.frame(
    taxon = taxon_vec,
    col   = ((seq_along(taxon_vec) - 1) %% n_cols) + 1,
    row   = ((seq_along(taxon_vec) - 1) %/% n_cols) + 1,
    stringsAsFactors = FALSE
  )
  grid_df$fill_col <- palette[match(grid_df$taxon, names(palette))]
  grid_df$taxon    <- factor(grid_df$taxon, levels = agg[[level]])

  # 9. 画图
  if (is.null(legend_title)) legend_title <- level

  p <- ggplot2::ggplot(grid_df,
                       ggplot2::aes(x = col, y = row, fill = taxon)) +
    ggplot2::geom_tile(
      color     = "white",
      linewidth = 0.8
    ) +
    ggplot2::scale_fill_manual(
      name   = legend_title,
      values = palette,
      guide  = ggplot2::guide_legend(ncol = 1)
    ) +
    ggplot2::scale_y_reverse() +
    ggplot2::coord_fixed() +
    ggplot2::theme_void(base_size = 11) +
    ggplot2::theme(
      legend.position = "right",
      legend.title    = ggplot2::element_text(face = "bold", size = 10),
      legend.text     = ggplot2::element_text(size = 9)
    )

  return(p)
}
