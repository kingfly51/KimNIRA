# -----------------------------------------------------------------------------
# kimPlotNIRA() -- Mean projection plot with significance annotation
# -----------------------------------------------------------------------------

#' Plot KIM NodeIdentifyR Algorithm (NIRA) Node Means with CI and Significance Stars
#'
#' @description
#' Produces a publication-ready plot matching the style of
#' \code{NIRApost::plotMeanNIRA} / Lunansky et al. (2022) Figure 3:
#' mean projected activation rate per node, 95\% CI error bars,
#' significance stars from corrected permutation p-values.
#'
#' Sorting:
#' \itemize{
#'   \item \code{"alleviating"}: ascending (lowest = most effective).
#'   \item \code{"aggravating"}: descending (highest = strongest).
#' }
#' \code{"baseline"} always leftmost.
#'
#' Stars: \code{***} p<0.001, \code{**} p<0.01, \code{*} p<0.05.
#'
#' @param statresult Output of \code{\link{kimPermutationNIRAtest}}.
#' @param perturbation_type \code{"alleviating"} or \code{"aggravating"}.
#' @param node_labels Named character vector for custom labels.
#'   Default \code{NULL}.
#'
#' @return A \code{ggplot2} object.
#'
#' @examples
#' \dontrun{
#' res <- runKIMNIRA(kim_fit, intervention_type = "alleviating")
#' stat_result <- kimPermutationNIRAtest(res)
#' kimPlotNIRA(stat_result, "alleviating")
#'
#' # With custom labels
#' kimPlotNIRA(stat_result, "alleviating",
#'             node_labels = c(Dep = "Depression", Anx = "Anxiety"))
#' }
#'
#' @seealso \code{\link{kimPermutationNIRAtest}}, \code{\link{kimPermutationTest}}
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_errorbar
#'   geom_text scale_x_continuous labs theme_minimal theme
#'   element_text element_blank
#' @export
kimPlotNIRA <- function(statresult,
                        perturbation_type,
                        node_labels = NULL) {

  if (!perturbation_type %in% c("aggravating", "alleviating"))
    stop("perturbation_type must be 'aggravating' or 'alleviating'.")
  if (!all(c("plot_data", "stat") %in% names(statresult)))
    stop("statresult must be produced by kimPermutationNIRAtest().")
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required.")

  pd   <- statresult$plot_data
  stat <- statresult$stat

  # -- Replace node labels -----------------------------------------------------
  if (!is.null(node_labels)) {
    for (nm in names(node_labels)) {
      pd$node[pd$node == nm] <- node_labels[nm]
      if (nm %in% rownames(stat))
        rownames(stat)[rownames(stat) == nm] <- node_labels[nm]
    }
  }

  # -- Sort --------------------------------------------------------------------
  is_bl     <- pd$node == "baseline"
  bl_row    <- pd[is_bl, , drop = FALSE]
  other_row <- pd[!is_bl, , drop = FALSE]

  if (perturbation_type == "alleviating") {
    other_row <- other_row[order(other_row$mean), ]
  } else {
    other_row <- other_row[order(other_row$mean, decreasing = TRUE), ]
  }

  od       <- rbind(bl_row, other_row)
  od$newid <- seq(0L, nrow(od) - 1L)
  od$node  <- factor(od$node, levels = od$node)

  # -- Merge significance ------------------------------------------------------
  sig_df <- data.frame(
    node     = rownames(stat),
    p.adjust = stat$p.adjust,
    stringsAsFactors = FALSE
  )
  if (!is.null(node_labels)) {
    for (nm in names(node_labels))
      sig_df$node[sig_df$node == nm] <- node_labels[nm]
  }

  od <- merge(od, sig_df, by = "node", all.x = TRUE)
  od <- od[order(od$newid), ]

  od$stars <- ifelse(is.na(od$p.adjust), "",
                     ifelse(od$p.adjust < 0.001, "***",
                            ifelse(od$p.adjust < 0.01,  "**",
                                   ifelse(od$p.adjust < 0.05,  "*", ""))))

  y_range  <- max(od$ciUpper, na.rm = TRUE) - min(od$ciLower, na.rm = TRUE)
  od$star_y <- od$ciUpper + y_range * 0.04

  # -- Plot --------------------------------------------------------------------
  gg <- ggplot2::ggplot(od, ggplot2::aes(x = .data$newid,
                                         y = .data$mean,
                                         group = 1)) +
    ggplot2::geom_line(color = "steelblue", linewidth = 1) +
    ggplot2::geom_point(size = 3, color = "tomato") +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = .data$ciLower, ymax = .data$ciUpper),
      width = 0.2, color = "gray40"
    ) +
    ggplot2::geom_text(
      ggplot2::aes(label = .data$stars, y = .data$star_y),
      color = "black", size = 5, vjust = 0.5
    ) +
    ggplot2::scale_x_continuous(
      breaks = od$newid,
      labels = as.character(od$node)
    ) +
    ggplot2::labs(
      x     = "Intervention Target",
      y     = "Mean Projected Activation Rate",
      title = ""
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.title.x       = ggplot2::element_text(face = "bold", size = 12),
      axis.title.y       = ggplot2::element_text(face = "bold", size = 12),
      axis.text.x        = ggplot2::element_text(angle = 45, hjust = 1,
                                                 face = "bold", size = 10),
      axis.text.y        = ggplot2::element_text(face = "bold", size = 10),
      panel.grid.major.x = ggplot2::element_blank()
    )

  gg
}
