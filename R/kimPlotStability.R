#' Stacked Bar Plot of Node Ranking Stability from KIM NodeIdentifyR Algorithm (NIRA) Stability Test
#'
#' @description
#' Visualises the output of \code{\link{kimStabilityNIRAtest}} as a horizontal
#' stacked bar plot.  Each bar represents one intervention node; stacked
#' segments show the \strong{proportion of NodeIdentifyR Algorithm (NIRA) runs}
#' in which that node occupied each rank position (from rank \code{low} to
#' rank \code{high}).
#'
#' This directly mirrors \code{plotMaxN()} from the \pkg{nodeIdentifyR}
#' ecosystem, adapted for the KIM framework.
#'
#' The plot communicates \emph{stability} at a glance:
#' \itemize{
#'   \item A node with a long dark segment (rank 1) that dominates all others
#'     is a \strong{highly stable} intervention target -- it consistently
#'     produced the largest simulated intervention effect across independent
#'     NodeIdentifyR Algorithm (NIRA) runs.
#'   \item A node whose colour segments are spread roughly equally across
#'     ranks is \strong{unstable} -- its position varies widely across runs.
#' }
#'
#' Bars are ordered so that the most stable node (highest proportion at rank 1,
#' ties broken by rank 2, rank 3, ?) appears at the \strong{bottom}.
#' Colours follow the Nature Publishing Group (NPG) palette from \pkg{ggsci}.
#'
#' @param stability_result A \code{"kim_stability_nira_result"} object from
#'   \code{\link{kimStabilityNIRAtest}}.
#' @param low Integer. First (lowest) rank position to display. Default:
#'   \code{1}.
#' @param high Integer. Last (highest) rank position to display. Must satisfy
#'   \code{high >= low} and \code{high <= p}. Default: \code{min(3, p)}.
#' @param label_ranks Logical. If \code{TRUE} (default), legend labels are
#'   shown as \code{"Rank 1"}, \code{"Rank 2"}, etc.
#'
#' @return A \code{\link[ggplot2]{ggplot}} object.
#'
#' @examples
#' \dontrun{
#' fit  <- fitKIM(data, vars = paste0("V", 1:6), id = "id", time = "time")
#' stab <- kimStabilityNIRAtest(fit, n_runs = 20)
#'
#' # Default: show top-3 ranks
#' kimPlotStability(stab)
#'
#' # Show all ranks (p = 6)
#' kimPlotStability(stab, high = 6)
#'
#' # Customise further
#' library(ggplot2)
#' p <- kimPlotStability(stab, high = 6) +
#'   labs(title = "KIM-NIRA Intervention Target Stability") +
#'   theme(plot.title = element_text(face = "bold", size = 14))
#' ggsave("Figure_stability.pdf", p, width = 7, height = 5, dpi = 300)
#' }
#'
#' @seealso
#' \itemize{
#'   \item \code{\link{kimStabilityNIRAtest}} -- produces the required input.
#'   \item \code{\link[ggsci]{scale_fill_npg}} -- the NPG colour palette.
#'   \item \code{\link[ggplot2]{ggsave}} -- save the returned plot.
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_bar scale_x_continuous labs
#'   theme_minimal theme element_text element_blank
#' @importFrom tidyr pivot_longer
#' @importFrom ggsci scale_fill_npg
#' @export
kimPlotStability <- function(stability_result,
                              low         = 1L,
                              high        = NULL,
                              label_ranks = TRUE) {

  # -- Parameter validation ----------------------------------------------------
  if (!inherits(stability_result, "kim_stability_nira_result"))
    stop("`stability_result` must be produced by kimStabilityNIRAtest() ",
         "(class = 'kim_stability_nira_result').")

  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required.")
  if (!requireNamespace("tidyr", quietly = TRUE))
    stop("Package 'tidyr' is required.")
  if (!requireNamespace("ggsci", quietly = TRUE))
    stop("Package 'ggsci' is required.")

  vars   <- stability_result$kim_fit$vars
  p      <- length(vars)
  n_runs <- nrow(stability_result$run_rankings)   # number of valid runs

  low  <- as.integer(low)
  if (is.null(high)) high <- min(3L, p)
  high <- as.integer(high)

  if (low  < 1L || low  > p) stop("`low` must be between 1 and p (", p, ").")
  if (high < low)             stop("`high` must be >= `low`.")
  if (high > p)               stop("`high` must be <= p (", p, ").")

  # -- Extract proportion matrix, select columns low:high ---------------------
  # rank_proportions[j, r] = proportion of runs where node j was ranked r
  rank_mat  <- stability_result$rank_proportions   # p x p, named rows & cols

  pct_cols  <- paste0("rank_", low:high)
  sub       <- as.data.frame(rank_mat[, pct_cols, drop = FALSE])
  sub$node  <- rownames(sub)

  # -- Node ordering (lexicographic strategy mirroring plotMaxN) --------------
  # Primary key: rank_low descending; secondary key: rank_(low+1), etc.
  # The lexicographically largest node is placed last in factor levels
  # -> appears at the bottom of the chart (most prominent position)
  sort_data <- sub[, pct_cols, drop = FALSE]
  sort_args <- lapply(rev(pct_cols), function(col) sort_data[[col]])
  ord_idx   <- do.call(order, c(sort_args, list(decreasing = FALSE)))
  ord_labels <- sub$node[ord_idx]

  sub$node <- factor(sub$node, levels = ord_labels)

  # -- Convert to long format --------------------------------------------------
  dl      <- tidyr::pivot_longer(sub,
                                  cols      = pct_cols,
                                  names_to  = "Rank",
                                  values_to = "Pct")
  dl$Pct  <- dl$Pct * 100   # convert to percentage

  # Human-readable rank labels: "rank_1" -> "Rank 1"
  if (label_ranks) {
    dl$Rank <- gsub("rank_(\\d+)", "Rank \\1", dl$Rank)
    dl$Rank <- factor(dl$Rank, levels = paste("Rank", low:high))
  } else {
    dl$Rank <- factor(dl$Rank, levels = pct_cols)
  }

  # -- Subtitle: report run count and top-1 node -------------------------------
  top1_node   <- stability_result$stable_ranking[1]
  top1_pct    <- round(rank_mat[top1_node, "rank_1"] * 100, 1)
  subtitle_txt <- ""

  # -- Plot --------------------------------------------------------------------
  ggplot2::ggplot(dl, ggplot2::aes(x = Pct, y = node, fill = Rank)) +
    ggplot2::geom_bar(stat = "identity", position = "stack") +
    ggsci::scale_fill_npg() +
    ggplot2::scale_x_continuous(
      labels = function(x) paste0(x, "%"),
      expand = c(0, 0),
      limits = c(0, NA)
    ) +
    ggplot2::labs(
      x        = "Cumulative occurrence percentage across NIRA runs",
      y        = "Intervention node",
      fill     = "Rank",
      subtitle = subtitle_txt
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      text               = ggplot2::element_text(size = 12),
      legend.position    = "right",
      panel.grid.major.y = ggplot2::element_blank(),
      plot.subtitle      = ggplot2::element_text(size = 9, colour = "grey40")
    )
}

# Suppress R CMD check NOTE for ggplot2 non-standard evaluation variables
utils::globalVariables(c("Pct", "node", "Rank"))
