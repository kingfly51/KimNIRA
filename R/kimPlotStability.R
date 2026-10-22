#' Stacked Bar Plot of Node Ranking Stability (sampling bootstrap)
#'
#' @description
#' Visualises the output of \code{\link{kimStabilityBootstrap}} as a horizontal
#' stacked bar plot.  Each bar represents one intervention node; stacked
#' segments show the \strong{proportion of person-bootstrap resamples} in which
#' that node occupied each rank position (from rank \code{low} to rank
#' \code{high}). This is the \emph{sampling} stability of the target ranking.
#'
#' This directly mirrors \code{plotMaxN()} from the \pkg{nodeIdentifyR}
#' ecosystem, adapted for the KIM framework.
#'
#' The plot communicates \emph{stability} at a glance:
#' \itemize{
#'   \item A node with a long first-rank segment that dominates all others
#'     is a \strong{highly stable} intervention target -- it is ranked 1st in
#'     most person-bootstrap resamples.
#'   \item A node whose colour segments are spread roughly equally across
#'     ranks is \strong{unstable} -- its position varies across resamples.
#' }
#'
#' Bars are ordered so that the most stable node (highest proportion at rank 1,
#' ties broken by rank 2, rank 3, and so on) appears at the \strong{bottom}.
#' Rank segments use a sequential palette (dark = rank 1, lighter = lower ranks),
#' so colour order matches rank order.
#'
#' @param stability_result The output of \code{\link{kimStabilityBootstrap}}
#'   (a \code{"kim_stability_bootstrap"} object, engine-agnostic).
#' @param low Integer. First (lowest) rank position to display. Default:
#'   \code{1}.
#' @param high Integer. Last (highest) rank position to display. Must satisfy
#'   \code{high >= low} and \code{high <= p}. Default: \code{min(3, p)}.
#' @param label_ranks Logical. If \code{TRUE} (default), legend labels are
#'   shown as \code{"Rank 1"}, \code{"Rank 2"}, etc.
#' @param colors Optional colours, unnamed or named rank_1, rank_2, etc.
#'   Short palettes are extended automatically.
#'
#' @return A \code{\link[ggplot2]{ggplot}} object.
#'
#' @examples
#' \dontrun{
#' fit  <- fitKIM(data, vars = paste0("V", 1:6), id = "id", time = "time")
#' res  <- runKIMNIRA(fit, method = "exact")     # or method = "mc"
#' boot <- kimStabilityBootstrap(res, data)
#'
#' # Default: show top-3 ranks
#' kimPlotStability(boot)
#'
#' # Show all ranks (p = 6)
#' kimPlotStability(boot, high = 6)
#'
#' # Customise further
#' library(ggplot2)
#' p <- kimPlotStability(boot, high = 6) +
#'   labs(title = "KIM-NIRA Intervention Target Stability") +
#'   theme(plot.title = element_text(face = "bold", size = 14))
#' ggsave("Figure_stability.pdf", p, width = 7, height = 5, dpi = 300)
#' }
#'
#' @seealso
#' \itemize{
#'   \item \code{\link{kimStabilityBootstrap}} -- produces the required input.
#'   \item \code{\link[ggplot2]{ggsave}} -- save the returned plot.
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_bar scale_x_continuous labs
#'   theme_minimal theme element_text element_blank
#' @importFrom tidyr pivot_longer
#' @export
kimPlotStability <- function(stability_result,
                              low         = 1L,
                              high        = NULL,
                              label_ranks = TRUE,
                              colors = NULL) {

  # -- Parameter validation ----------------------------------------------------
  # Accept the SAMPLING bootstrap output (kimStabilityBootstrap, engine exact
  # or mc) and convert it into the same kim_stability_nira_result structure
  # (rank proportions per resample) that kimStabilityBootstrap() produces, so
  # the drawing code below is reused unchanged (identical figure style).
  if (inherits(stability_result, "kim_stability_bootstrap")) {
    fb <- stability_result$kim_fit
    effM <- stability_result$resample_effects          # B x p
    vars <- fb$vars; p <- length(vars)
    sign <- if(stability_result$parameters$intervention_type=="alleviating")1 else -1
    rr <- t(apply(effM, 1, function(r) rank(sign*r, ties.method = "average")))
    colnames(rr) <- vars
    prop <- stability_result$rank_proportions
    weight <- (p - rr + 1L)
    tot <- colSums(weight)
    stability_result <- list(
      kim_fit          = fb,
      run_rankings     = rr,
      rank_proportions = prop,
      stable_ranking   = vars[order(tot, decreasing = TRUE)],
      weighted_scores  = tot,
      score_df         = data.frame(node = vars, total_score = tot))
    class(stability_result) <- "kim_stability_nira_result"
  }
  if (!inherits(stability_result, "kim_stability_nira_result"))
    stop("`stability_result` must be the output of kimStabilityBootstrap().")

  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required.")
  if (!requireNamespace("tidyr", quietly = TRUE))
    stop("Package 'tidyr' is required.")

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
  # The lexicographically largest node is first in factor levels,
  # which places it at the bottom of a discrete y axis.
  sort_data <- sub[, pct_cols, drop = FALSE]
  sort_args <- lapply(pct_cols, function(col) sort_data[[col]])
  ord_idx   <- do.call(order, c(sort_args, list(decreasing = TRUE)))
  ord_labels <- sub$node[ord_idx]

  sub$node <- factor(sub$node, levels = ord_labels)

  # -- Convert to long format --------------------------------------------------
  dl      <- tidyr::pivot_longer(sub,
                                  cols      = -"node",
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

  # -- Plot --------------------------------------------------------------------
  ggplot2::ggplot(dl, ggplot2::aes(x = Pct, y = node, fill = Rank)) +
    ggplot2::geom_bar(stat = "identity", position = "stack") +
    ggplot2::scale_fill_manual(values = if (length(colors)) unname(.kim_plot_colors(
      paste0("rank_", seq_len(p)), colors, character(0))[low:high]) else
      grDevices::hcl.colors(high - low + 2L, "Blues 3")[seq_len(high - low + 1L)]) +
    ggplot2::scale_x_continuous(
      labels = function(x) paste0(x, "%"),
      expand = c(0, 0),
      limits = c(0, NA)
    ) +
    ggplot2::labs(
      x        = "Share of person-bootstrap resamples (model re-estimated in each)",
      y        = "Intervention node",
      fill     = "Rank"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      text               = ggplot2::element_text(size = 12),
      legend.position    = "right",
      panel.grid.major.y = ggplot2::element_blank()
    )
}

# Suppress R CMD check NOTE for ggplot2 non-standard evaluation variables
utils::globalVariables(c("Pct", "node", "Rank"))
