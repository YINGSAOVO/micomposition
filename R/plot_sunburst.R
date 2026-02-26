#' @importFrom rlang .data
utils::globalVariables(c("xmin", "xmax", "ymin", "ymax", "label",
                         "angle", "hjust", "mid", "fill_color",
                         "text_angle", "text_hjust"))

#' Sunburst chart for microbiome composition
#'
#' Draws a multi-ring sunburst chart showing hierarchical taxonomic composition
#' for a single sample or the mean across all samples.
#'
#' @param data A list with elements `otu_table`, `tax_table`, and `sample_data`,
#'   as in [demo_data].
#' @param sample Name of the sample to display. If `NULL`, uses the mean
#'   abundance across all samples.
#' @param levels Character vector of taxonomic levels to display from inner to
#'   outer ring. Default is `c("Phylum", "Family")`.
#' @param min_abundance Minimum mean abundance to display a slice. Taxa below
#'   this threshold are merged into "Others". Default is `0.01`.
#' @param palette A named or unnamed character vector of colors. If `NULL`,
#'   uses a built-in palette.
#' @param show_labels Logical. Whether to show text labels on slices.
#'   Default is `TRUE`.
#'
#' @return A ggplot object.
#' @export
#'
#' @examples
#' data(demo_data)
#' plot_sunburst(demo_data)
#' plot_sunburst(demo_data, sample = "Sample1")
#' plot_sunburst(demo_data, levels = c("Phylum", "Class", "Family"))
plot_sunburst <- function(data,
                          sample        = NULL,
                          levels        = c("Phylum", "Family"),
                          min_abundance = 0.01,
                          palette       = NULL,
                          show_labels   = TRUE) {

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

  # 3. 合并丰度与分类信息
  abu_vec <- unlist(abu)
  df <- data.frame(
    OTU       = names(abu_vec),
    Abundance = as.numeric(abu_vec),
    stringsAsFactors = FALSE
  )
  df <- merge(df, tax, by = "OTU")

  # 4. 配色：按最内层分类着色
  all_top <- unique(df[[levels[1]]])
  default_pal <- c(
    "#4E79A7","#F28E2B","#E15759","#76B7B2","#59A14F",
    "#EDC948","#B07AA1","#FF9DA7","#9C755F","#BAB0AC"
  )
  if (is.null(palette)) {
    palette <- stats::setNames(
      rep(default_pal, length.out = length(all_top)),
      all_top
    )
  }

  # 辅助函数：根据分类名查颜色
  get_color <- function(taxon, alpha_val = 1) {
    if (is.na(taxon) || grepl("^Others", taxon)) return("#D3D3D3")
    col <- palette[taxon]
    if (is.na(col)) return("#D3D3D3")
    scales::alpha(unname(col), alpha_val)
  }

  # 5. 逐层构建扇形数据
  n_rings   <- length(levels)
  ring_list <- list()

  for (i in seq_along(levels)) {
    lv <- levels[i]

    if (i == 1) {
      # 最内层
      agg <- stats::aggregate(Abundance ~ get(lv), data = df, FUN = sum)
      colnames(agg)[1] <- lv
      agg[[lv]] <- ifelse(agg$Abundance >= min_abundance, agg[[lv]], "Others")
      agg <- stats::aggregate(Abundance ~ get(lv), data = agg, FUN = sum)
      colnames(agg)[1] <- lv
      agg$Abundance <- agg$Abundance / sum(agg$Abundance)
      non_others <- agg[agg[[lv]] != "Others", ]
      non_others <- non_others[order(non_others$Abundance, decreasing = TRUE), ]
      others     <- agg[agg[[lv]] == "Others", ]
      agg        <- rbind(non_others, others)

      agg$ymax <- cumsum(agg$Abundance)
      agg$ymin <- c(0, head(agg$ymax, -1))
      agg$mid  <- (agg$ymin + agg$ymax) / 2

      agg$fill_color <- sapply(agg[[lv]], get_color, USE.NAMES = FALSE)
      agg$ancestor   <- agg[[lv]]   # 最内层自己就是祖先
      agg$label      <- agg[[lv]]

    } else {
      # 外层：在父层弧度范围内分配子扇形
      parent_lv <- levels[i - 1]
      parent_df <- ring_list[[i - 1]]

      # 建立 lv -> 最内层祖先 的映射（不用 merge，用 match 避免行序乱）
      ancestor_map <- unique(df[, c(lv, levels[1]), drop = FALSE])

      all_rows <- list()

      for (p_idx in seq_len(nrow(parent_df))) {
        p_label <- parent_df$label[p_idx]
        p_ymin  <- parent_df$ymin[p_idx]
        p_ymax  <- parent_df$ymax[p_idx]
        p_span  <- p_ymax - p_ymin
        p_ancestor <- parent_df$ancestor[p_idx]

        # 找属于这个父节点的子节点
        if (grepl("^Others", p_label)) {
          # 父节点是 Others：子节点也统一为 Others
          child_row <- data.frame(
            label      = "Others",
            ymin       = p_ymin,
            ymax       = p_ymax,
            mid        = (p_ymin + p_ymax) / 2,
            fill_color = "#D3D3D3",
            ancestor   = "Others",
            stringsAsFactors = FALSE
          )
          all_rows[[length(all_rows) + 1]] <- child_row
          next
        }

        # 找 df 里父层级 == p_label 的行，聚合当前层
        sub_df <- df[df[[parent_lv]] == p_label, ]
        if (nrow(sub_df) == 0) next

        children <- stats::aggregate(Abundance ~ get(lv), data = sub_df, FUN = sum)
        colnames(children)[1] <- lv

        # 低丰度归 Others
        thresh <- min_abundance * sum(df$Abundance)
        children[[lv]] <- ifelse(children$Abundance >= thresh, children[[lv]], "Others")
        children <- stats::aggregate(Abundance ~ get(lv), data = children, FUN = sum)
        colnames(children)[1] <- lv

        # 排序：非 Others 按丰度降序，Others 排最后
        non_oth <- children[children[[lv]] != "Others", ]
        oth     <- children[children[[lv]] == "Others", ]
        non_oth <- non_oth[order(non_oth$Abundance, decreasing = TRUE), ]
        children <- rbind(non_oth, oth)

        # 在父扇形范围内按比例分配坐标
        children$prop <- children$Abundance / sum(children$Abundance)
        children$ymax <- p_ymin + cumsum(children$prop) * p_span
        children$ymin <- c(p_ymin, head(children$ymax, -1))
        children$mid  <- (children$ymin + children$ymax) / 2

        # 取色：通过最内层祖先
        children$fill_color <- sapply(children[[lv]], function(x) {
          if (is.na(x) || x == "Others") return("#D3D3D3")
          # 查这个 taxon 属于哪个最内层分类
          ancestor_val <- ancestor_map[[levels[1]]][match(x, ancestor_map[[lv]])]
          get_color(ancestor_val, alpha_val = 0.95 - 0.2 * (i - 1))
        }, USE.NAMES = FALSE)

        children$ancestor <- sapply(children[[lv]], function(x) {
          if (is.na(x) || x == "Others") return("Others")
          ancestor_map[[levels[1]]][match(x, ancestor_map[[lv]])]
        }, USE.NAMES = FALSE)

        child_rows <- data.frame(
          label      = children[[lv]],
          ymin       = children$ymin,
          ymax       = children$ymax,
          mid        = children$mid,
          fill_color = children$fill_color,
          ancestor   = children$ancestor,
          stringsAsFactors = FALSE
        )
        all_rows[[length(all_rows) + 1]] <- child_rows
      }

      agg <- do.call(rbind, all_rows)
    }

    ring_width <- 0.85
    gap        <- 0.08

    ring_df <- data.frame(
      xmin       = (i - 1) * (ring_width + gap) + 0.5,
      xmax       = (i - 1) * (ring_width + gap) + 0.5 + ring_width,
      ymin       = agg$ymin,
      ymax       = agg$ymax,
      mid        = agg$mid,
      label      = agg$label,
      fill_color = agg$fill_color,
      ancestor   = agg$ancestor,
      ring       = i,
      stringsAsFactors = FALSE
    )
    ring_df[[lv]]  <- agg$label
    ring_list[[i]] <- ring_df
  }

  # 6. 合并所有层数据
  plot_df <- do.call(rbind, lapply(ring_list, function(d) {
    data.frame(
      xmin       = d$xmin,
      xmax       = d$xmax,
      ymin       = d$ymin,
      ymax       = d$ymax,
      mid        = d$mid,
      label      = d$label,
      fill_color = d$fill_color,
      ring       = d$ring,
      stringsAsFactors = FALSE
    )
  }))

  max_x <- max(plot_df$xmax)

  # 7. 计算标签角度
  plot_df$text_angle <- 90 - 360 * plot_df$mid
  plot_df$text_hjust <- ifelse(
    plot_df$text_angle < -90 | plot_df$text_angle > 90, 1, 0
  )
  plot_df$text_angle <- ifelse(
    plot_df$text_angle < -90,
    plot_df$text_angle + 180,
    plot_df$text_angle
  )

  # 8. 画扇形
  p <- ggplot2::ggplot(plot_df) +
    ggplot2::geom_rect(
      ggplot2::aes(
        xmin = xmin, xmax = xmax,
        ymin = ymin, ymax = ymax,
        fill = label
      ),
      color = "white", linewidth = 0.5
    ) +
    ggplot2::scale_fill_manual(
      values = stats::setNames(plot_df$fill_color, plot_df$label),
      breaks = plot_df$label[plot_df$ring == 1]
    ) +
    ggplot2::coord_polar(theta = "y", start = 0, clip = "off") +
    ggplot2::xlim(c(0, max_x + 1.8)) +
    ggplot2::theme_void(base_size = 11) +
    ggplot2::theme(
      legend.position  = "right",
      legend.title     = ggplot2::element_blank(),
      legend.text      = ggplot2::element_text(size = 9),
      legend.key.size  = ggplot2::unit(0.45, "cm"),
      legend.spacing.y = ggplot2::unit(0.15, "cm")
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(ncol = 1, byrow = TRUE)
    )

  # 9. 标签
  if (show_labels) {
    outer_df <- plot_df[
      plot_df$ring == n_rings & (plot_df$ymax - plot_df$ymin) > 0.035, ]
    if (nrow(outer_df) > 0) {
      p <- p + ggplot2::geom_text(
        data = outer_df,
        ggplot2::aes(
          x     = xmax + 0.25,
          y     = mid,
          label = label,
          angle = text_angle,
          hjust = text_hjust
        ),
        size = 2.8, color = "grey20"
      )
    }

    if (n_rings >= 3) {
      mid_df <- plot_df[
        plot_df$ring > 1 & plot_df$ring < n_rings &
          (plot_df$ymax - plot_df$ymin) > 0.07, ]
      if (nrow(mid_df) > 0) {
        p <- p + ggplot2::geom_text(
          data = mid_df,
          ggplot2::aes(
            x     = (xmin + xmax) / 2,
            y     = mid,
            label = label,
            angle = text_angle
          ),
          size = 2.5, color = "white", fontface = "bold"
        )
      }
    }
  }

  return(p)
}
