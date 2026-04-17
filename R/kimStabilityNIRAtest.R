# =============================================================================
# kimStabilityNIRAtest.R  ----  KIM NodeIdentifyR Algorithm (NIRA) Stability Test
#
# Functions:
#   kimStabilityNIRAtest()           Repeat NodeIdentifyR Algorithm (NIRA) runs to assess node ranking stability
#   print.kim_stability_nira_result  Print method
#
# Design notes:
#   Closely mirrors the logic of stabilityNIRAtest() / findMaxN():
#     - Accepts kim_fit (class "dynaNIRA_kim") as input
#     - Runs runKIMNIRA() multiple times with different random seeds
#     - Tabulates the proportion of runs in which each node occupied each
#       rank position (analogous to findMaxN() output)
#         -> does NOT output Kendall's tau; proportions are the stability metric
#     - Aggregates rankings across all nodes using a position-weighted scheme
#       (to derive stable_ranking):
#         if there are p nodes, rank 1 scores p points, rank 2 scores p-1, ...,
#         rank p scores 1 point
#     - Output: per-node rank-proportion matrix, stable ranking, weighted scores
# =============================================================================


#' Stability Test for KIM NodeIdentifyR Algorithm (NIRA) Node Rankings
#'
#' @description
#' Runs \code{\link{runKIMNIRA}} \code{n_runs} times with independent random
#' seeds and assesses the stability of node rankings by computing, for each
#' node, the \strong{proportion of runs} in which it occupied each rank
#' position (rank 1 = best intervention target, rank p = worst).
#'
#' This mirrors the \code{stabilityNIRAtest} / \code{findMaxN} workflow from
#' the \pkg{nodeIdentifyR} ecosystem, adapted for the KIM framework.
#'
#' The consensus ranking (returned as \code{stable_ranking}) is derived via a
#' position-weighted score: a node ranked 1st scores \eqn{p} points, rank 2nd
#' scores \eqn{p-1}, ..., rank \eqn{p} scores 1 point; scores are summed across
#' runs and nodes are sorted descending.
#'
#' @param kim_fit A \code{"dynaNIRA_kim"} object from \code{\link{fitKIM}}.
#' @param intervention_type Character, \code{"alleviating"} (default) or
#'   \code{"aggravating"}.
#' @param n_runs Positive integer. Number of independent NodeIdentifyR
#'   Algorithm (NIRA) runs. Default \code{10L}. More runs yield more reliable
#'   frequency estimates; 10-20 is usually sufficient.
#' @param n_sim Positive integer. Simulations per NodeIdentifyR Algorithm
#'   (NIRA) run. Default \code{5000L}.
#' @param n_steps Positive integer. Trajectory steps per simulation. Default
#'   \code{30L}.
#' @param intervention_sd Numeric. Intervention strength multiplier. Default
#'   \code{2}.
#' @param base_seed Integer. Base random seed; run \eqn{k} uses
#'   \code{base_seed + (k - 1) * 1000L}. Default \code{2025L}.
#' @param verbose Logical. Print progress. Default \code{TRUE}.
#'
#' @return A named \code{list} of class \code{"kim_stability_nira_result"}:
#' \describe{
#'   \item{\code{stable_ranking}}{Character vector of node names ordered from
#'     most to least effective intervention target (by aggregated weighted score).}
#'   \item{\code{weighted_scores}}{Named numeric vector of aggregated position
#'     weights (higher = more consistently ranked effective).}
#'   \item{\code{rank_proportions}}{Matrix (\code{p x p}) where entry
#'     \code{[j, r]} is the proportion of runs in which node \code{j} was
#'     ranked \code{r}. Row names = node names; column names = \code{"rank_1"},
#'     \code{"rank_2"}, ..., \code{"rank_p"}.
#'     This is the KIM analogue of \code{findMaxN()} output.}
#'   \item{\code{score_df}}{Data frame with columns \code{node},
#'     \code{total_score}, \code{mean_score}, \code{mean_rank}, \code{sd_rank},
#'     and \code{final_rank}.}
#'   \item{\code{run_rankings}}{Matrix (\code{n_valid x p}) of node ranks
#'     (1 = best target) for each NodeIdentifyR Algorithm (NIRA) run,
#'     with node names as column names.}
#'   \item{\code{run_effects}}{Matrix (\code{n_valid x p}) of raw mean
#'     intervention effects per run.}
#'   \item{\code{nira_results}}{List of raw \code{kim_nira_result} objects
#'     from each NodeIdentifyR Algorithm (NIRA) run.}
#'   \item{\code{kim_fit}}{Reference to input \code{kim_fit}.}
#'   \item{\code{params}}{List of call parameters.}
#' }
#'
#' @examples
#' \dontrun{
#' fit <- fitKIM(data, vars = c("V1","V2","V3","V4","V5","V6"),
#'               id = "id", time = "time")
#' stab <- kimStabilityNIRAtest(fit, n_runs = 10)
#' print(stab)
#' kimPlotStability(stab, high = 6)
#' }
#'
#' @seealso \code{\link{kimPlotStability}}, \code{\link{runKIMNIRA}}
#' @importFrom stats sd setNames
#' @export
kimStabilityNIRAtest <- function(kim_fit,
                                 intervention_type = c("alleviating", "aggravating"),
                                 n_runs            = 10L,
                                 n_sim             = 5000L,
                                 n_steps           = 30L,
                                 intervention_sd   = 2,
                                 base_seed         = 2025L,
                                 verbose           = TRUE) {

  # -- Parameter validation ----------------------------------------------------
  # Accept "dynaNIRA_kim" (class name returned by fitKIM) or "kim_fit" (alias)
  if (!inherits(kim_fit, "dynaNIRA_kim") && !inherits(kim_fit, "kim_fit"))
    stop("`kim_fit` must be produced by fitKIM() ",
         "(class 'dynaNIRA_kim' or 'kim_fit').")

  intervention_type <- match.arg(intervention_type)

  n_runs  <- as.integer(n_runs)
  n_sim   <- as.integer(n_sim)
  n_steps <- as.integer(n_steps)

  if (n_runs < 2L)
    stop("`n_runs` must be >= 2 (need at least 2 runs to assess stability).")

  vars <- kim_fit$vars
  p    <- length(vars)

  if (verbose) {
    message("-- KIM NodeIdentifyR Algorithm (NIRA) Stability Test --")
    message("  nodes = ", p, "  |  n_runs = ", n_runs,
            "  |  n_sim = ", n_sim, "  |  n_steps = ", n_steps)
    message("  intervention_type = ", intervention_type,
            "  |  base_seed = ", base_seed)
    message("  Scoring: position weights p=", p, " down to 1")
  }

  # -- Storage for run results -------------------------------------------------
  run_rankings <- matrix(NA_real_, nrow = n_runs, ncol = p,
                         dimnames = list(paste0("run", seq_len(n_runs)), vars))
  run_effects  <- matrix(NA_real_, nrow = n_runs, ncol = p,
                         dimnames = list(paste0("run", seq_len(n_runs)), vars))
  nira_results <- vector("list", n_runs)

  # -- Run NodeIdentifyR Algorithm (NIRA) iteratively -------------------------
  for (k in seq_len(n_runs)) {

    run_seed <- base_seed + (k - 1L) * 1000L

    if (verbose)
      message("  [", k, "/", n_runs, "] seed = ", run_seed)

    res <- tryCatch(
      runKIMNIRA(kim_fit,
                 intervention_type = intervention_type,
                 n_sim             = n_sim,
                 n_steps           = n_steps,
                 intervention_sd   = intervention_sd,
                 seed              = run_seed,
                 verbose           = FALSE),
      error = function(e) {
        warning("Run ", k, " failed: ", conditionMessage(e))
        NULL
      }
    )

    if (is.null(res)) next

    nira_results[[k]] <- res

    # Get per-node mean intervention effects, aligned to vars order
    eff           <- res$summary$mean_diff
    names(eff)    <- res$summary$node
    eff_aligned   <- eff[vars]

    run_effects[k, ]  <- eff_aligned
    # Rankings: absolute effect from largest to smallest; rank 1 = strongest.
    # ties.method = "first" ensures integer ranks with no ties (fractional
    # ranks from the default "average" method would cause sum(node_ranks == r)
    # to return 0 for the tied positions, making row sums of rank_proportions
    # fall below 1 and the stacked bar chart fail to reach 100%).
    run_rankings[k, ] <- rank(-abs(eff_aligned), ties.method = "first")
  }

  # -- Remove failed runs ------------------------------------------------------
  valid_runs <- which(apply(run_rankings, 1, function(r) !any(is.na(r))))
  if (length(valid_runs) < 2L)
    stop("Fewer than 2 runs succeeded. Cannot assess stability.")

  if (length(valid_runs) < n_runs && verbose)
    message("  Warning: ", n_runs - length(valid_runs), " run(s) failed; ",
            "using ", length(valid_runs), " valid runs.")

  run_rankings_valid <- run_rankings[valid_runs, , drop = FALSE]
  run_effects_valid  <- run_effects[valid_runs, , drop = FALSE]
  nira_results_valid <- nira_results[valid_runs]
  n_valid            <- length(valid_runs)

  # -- Per-node rank-proportion matrix (core statistic, analogous to findMaxN)
  #
  # rank_proportions[j, r] = proportion of runs where node j was ranked r

  rank_proportions <- matrix(0.0, nrow = p, ncol = p,
                             dimnames = list(vars,
                                            paste0("rank_", seq_len(p))))
  for (j in seq_len(p)) {
    node_ranks <- run_rankings_valid[, j]
    for (r in seq_len(p)) {
      rank_proportions[j, r] <- sum(node_ranks == r, na.rm = TRUE) / n_valid
    }
  }

  # -- Position-weighted aggregation (derive stable_ranking) ------------------
  #
  # weight for rank r = p - r + 1
  # (rank 1 -> p points, rank p -> 1 point)
  weight_matrix <- (p - run_rankings_valid + 1)  # n_valid x p
  total_scores  <- colSums(weight_matrix)
  mean_scores   <- colMeans(weight_matrix)

  mean_ranks <- colMeans(run_rankings_valid)
  sd_ranks   <- apply(run_rankings_valid, 2, sd)

  # Stable ranking: descending total score
  stable_order   <- order(total_scores, decreasing = TRUE)
  stable_ranking <- vars[stable_order]

  # -- score_df ----------------------------------------------------------------
  score_df <- data.frame(
    node        = vars,
    total_score = total_scores,
    mean_score  = round(mean_scores, 3),
    mean_rank   = round(mean_ranks, 3),
    sd_rank     = round(sd_ranks, 3),
    stringsAsFactors = FALSE
  )
  score_df$final_rank <- rank(score_df$mean_rank)
  score_df <- score_df[order(score_df$final_rank), ]
  rownames(score_df) <- NULL

  if (verbose) {
    message("  Stable ranking: ", paste(stable_ranking, collapse = " > "))
    message("  Top-1 proportions: ",
            paste(paste0(vars[stable_order], "=",
                         round(rank_proportions[stable_order, "rank_1"] * 100, 1), "%"),
                  collapse = ", "))
  }

  # -- Return object -----------------------------------------------------------
  result <- list(
    stable_ranking   = stable_ranking,
    weighted_scores  = setNames(total_scores, vars),
    rank_proportions = rank_proportions,
    score_df         = score_df,
    run_rankings     = run_rankings_valid,
    run_effects      = run_effects_valid,
    nira_results     = nira_results_valid,
    kim_fit          = kim_fit,
    params           = list(
      intervention_type = intervention_type,
      n_runs            = n_runs,
      n_valid           = n_valid,
      n_sim             = n_sim,
      n_steps           = n_steps,
      intervention_sd   = intervention_sd,
      base_seed         = base_seed
    )
  )
  class(result) <- "kim_stability_nira_result"
  result
}


# -----------------------------------------------------------------------------
# print method
# -----------------------------------------------------------------------------

#' @export
print.kim_stability_nira_result <- function(x, ...) {
  cat("\n-- KIM NodeIdentifyR Algorithm (NIRA) Stability Test Result --\n\n")

  p <- length(x$kim_fit$vars)
  cat("  Nodes (p)     :", p, "\n")
  cat("  Runs (valid)  :", x$params$n_valid, "/", x$params$n_runs, "\n")
  cat("  n_sim         :", x$params$n_sim, "\n")
  cat("  n_steps       :", x$params$n_steps, "\n")
  cat("  Type          :", x$params$intervention_type, "\n")
  cat("  Scoring       : position weights p=", p, " down to 1\n\n")

  cat("  Stable Ranking:\n")
  cat("   ", paste(seq_along(x$stable_ranking), x$stable_ranking,
                   sep = ". ", collapse = "\n    "), "\n\n")

  cat("  Rank Proportions (proportion of runs at each rank):\n")
  rp <- round(x$rank_proportions * 100, 1)
  colnames(rp) <- sub("rank_", "Rank ", colnames(rp))
  print(as.data.frame(rp))

  cat("\n  Score Table (sorted by mean rank):\n")
  df <- x$score_df
  df$total_score <- round(df$total_score, 1)
  print(df, row.names = FALSE)
  cat("\n")
  invisible(x)
}
