# =============================================================================
# runKIMNIRA.R  ----  KIM NodeIdentifyR Algorithm (NIRA) Simulated Intervention Analysis
#
# This file can be sourced standalone without loading other files.
# Dependency: glmnet (already installed by fitKIM)
#
# Functions:
#   kim_sigmoid()              Sigmoid activation function
#   kim_compute_delta()        Compute intercept shift delta
#   kim_delta_prob_table()     Table of per-node probability changes
#   kim_init_states()          Sample initial states from independent distributions
#   kim_step()                 KIM single-step state update (independent sampling)
#   kim_trajectory()           Full trajectory simulation
#   runKIMNIRA()               Main function
#   print.kim_nira_result      Print method
#
# Design notes:
#   In KIM, the temporal network J captures cross-time influences. During
#   simulation, given s(t), nodes at t+1 are conditionally independent:
#
#     P(s_i(t+1) = 1 | s(t)) = sigma( sum_j J_ij * s_j(t) + h_i )
#
#   Sampling strategy: independent Bernoulli draws (nodes are conditionally
#   independent given s(t)).
#
#   Paired design: all random numbers are pre-generated and shared; the only
#   difference between baseline and intervention conditions is the intercept shift.
# =============================================================================


# -----------------------------------------------------------------------------
# Basic helper functions
# -----------------------------------------------------------------------------

#' Sigmoid Function
#'
#' Directly wraps R's built-in \code{plogis} for numerical stability and speed.
#' @keywords internal
kim_sigmoid <- function(x) {
  stats::plogis(x)
}


#' Compute Intervention Intercept Shift (delta)
#'
#' Follows Lunansky et al. (2022): delta = +/- intervention_sd * sd(h).
#' When sd(h) is degenerate (near zero), the median of |h| is used as a
#' fallback scale and a warning is issued.
#' @export
kim_compute_delta <- function(h,
                              intervention_sd  = 2,
                              type = c("alleviating", "aggravating")) {
  type <- match.arg(type)
  h_sd <- sd(h)

  if (h_sd < 1e-6) {
    # Degenerate fallback: use median of absolute intercepts
    fallback <- median(abs(h))
    # Final lower bound
    if (fallback < 1e-6) fallback <- 0.5

    h_sd <- fallback
    warning("Intercept vector SD is near zero. Using fallback scale ",
            round(h_sd, 4), " to compute delta.",
            "\n  Suggestion: control intervention strength directly via ",
            "`intervention_sd`, or specify delta manually.")
  }

  delta <- intervention_sd * h_sd
  if (type == "alleviating") -delta else delta
}


#' Compute Per-Node Activation Probability Change Table
#' @export
kim_delta_prob_table <- function(h, delta, vars) {
  data.frame(
    node         = vars,
    h_original   = round(h, 4),
    h_intervened = round(h + delta, 4),
    p_original   = round(kim_sigmoid(h), 4),
    p_intervened = round(kim_sigmoid(h + delta), 4),
    delta_p      = round(kim_sigmoid(h + delta) - kim_sigmoid(h), 4),
    stringsAsFactors = FALSE
  )
}


# -----------------------------------------------------------------------------
# Initial state sampling
# -----------------------------------------------------------------------------

#' Sample Initial States from Independent Distributions
#'
#' Each node is sampled independently: P(s_i = 1) = sigma(h_i)
#' @export
kim_init_states <- function(h, n_sim, vars = NULL, seed = NULL) {
  n_vars <- length(h)
  if (!is.null(seed)) set.seed(seed)

  init_probs <- kim_sigmoid(h)
  u     <- matrix(runif(n_sim * n_vars), nrow = n_sim, ncol = n_vars)
  p_mat <- matrix(init_probs, nrow = n_sim, ncol = n_vars, byrow = TRUE)
  states <- matrix(as.integer(u < p_mat), nrow = n_sim, ncol = n_vars)

  colnames(states) <- if (!is.null(vars)) vars else paste0("V", seq_len(n_vars))
  storage.mode(states) <- "integer"
  states
}


# -----------------------------------------------------------------------------
# KIM single-step state update
# -----------------------------------------------------------------------------

#' KIM Single-Step State Update (Independent Sampling)
#'
#' Given s(t), each node at t+1 is conditionally independent:
#'   P(s_i(t+1) = 1 | s(t)) = sigma( sum_j J_ij * s_j(t) + h_i )
#'
#' @export
kim_step <- function(s, J, h_eff, U_step = NULL) {
  n_sim  <- nrow(s)
  n_vars <- ncol(s)

  h_mat <- matrix(h_eff, nrow = n_sim, ncol = n_vars, byrow = TRUE)
  eta   <- s %*% t(J) + h_mat

  probs <- kim_sigmoid(eta)
  U <- if (!is.null(U_step)) U_step else
    matrix(runif(n_sim * n_vars), nrow = n_sim, ncol = n_vars)
  result <- matrix(as.integer(U < probs), nrow = n_sim, ncol = n_vars)
  colnames(result) <- colnames(s)
  result
}


# -----------------------------------------------------------------------------
# Pre-generate paired random numbers
# -----------------------------------------------------------------------------

#' @keywords internal
.kim_pregenerate_U <- function(n_sim, n_vars, n_steps, seed) {
  set.seed(seed)

  lapply(seq_len(n_steps), function(step) {
    matrix(runif(n_sim * n_vars), nrow = n_sim, ncol = n_vars)
  })
}


# -----------------------------------------------------------------------------
# Trajectory simulation
# -----------------------------------------------------------------------------

#' Simulate a Full KIM Trajectory
#'
#' @description
#' Runs a complete forward simulation of the Kinetic Ising Model for
#' \code{n_steps} time steps, starting from \code{init_states}, and returns
#' the mean projected activation rate (sum score averaged over steps and nodes)
#' for each simulated individual. This is the core trajectory primitive used
#' by \code{\link{runKIMNIRA}}.
#'
#' @param init_states Integer matrix (\code{n_sim x p}) of initial binary
#'   states, as produced by \code{\link{kim_init_states}}.
#' @param J Numeric matrix (\code{p x p}), the temporal edge weight matrix
#'   from \code{\link{fitKIM}} (\code{J[i,j]} = effect of j on i).
#' @param h Numeric vector (length p) of node thresholds (intercepts).
#' @param n_steps Positive integer. Number of time steps to simulate.
#' @param target_node Integer or \code{NULL}. Index of the intervention node.
#'   \code{NULL} (default) runs the baseline (no intervention).
#' @param delta Numeric. Threshold shift applied to \code{target_node}.
#'   Default \code{0} (no shift).
#' @param U_list List of pre-generated uniform random matrices (one per step),
#'   as produced by \code{.kim_pregenerate_U}. \code{NULL} generates fresh
#'   random numbers at each step.
#'
#' @return Numeric vector of length \code{n_sim}: each element is the mean
#'   projected activation rate (proportion of node-steps activated) for one
#'   simulated trajectory.
#'
#' @seealso \code{\link{runKIMNIRA}}, \code{\link{kim_step}},
#'   \code{\link{kim_init_states}}
#' @importFrom stats runif
#' @export
kim_trajectory <- function(init_states, J, h, n_steps,
                           target_node    = NULL,
                           delta          = 0,
                           U_list         = NULL) {
  n_sim  <- nrow(init_states)
  n_vars <- ncol(init_states)

  h_eff <- h
  if (!is.null(target_node) && delta != 0)
    h_eff[target_node] <- h[target_node] + delta

  s         <- init_states
  cum_score <- numeric(n_sim)

  for (step in seq_len(n_steps)) {
    U_step    <- if (!is.null(U_list)) U_list[[step]] else NULL
    s         <- kim_step(s, J, h_eff, U_step = U_step)
    cum_score <- cum_score + rowSums(s)
  }

  cum_score / (n_steps * n_vars)
}


# -----------------------------------------------------------------------------
# Main function
# -----------------------------------------------------------------------------

#' Run the NodeIdentifyR Algorithm (NIRA) on a Fitted KIM
#'
#' @description
#' Implements the NodeIdentifyR Algorithm (NIRA) via paired Monte Carlo
#' simulation. For each of the \eqn{p} nodes, an alleviating or aggravating
#' threshold shift \eqn{\delta} is applied to that node's intercept while all
#' temporal edge weights \eqn{J} remain fixed. The mean difference in projected
#' network activation rate between the intervention and baseline conditions is
#' used to rank nodes as intervention targets.
#'
#' All random numbers are pre-generated and shared between the baseline and
#' every intervention condition (paired design), so that differences reflect
#' only the threshold shift and not Monte Carlo noise.
#'
#' @param kim_fit A \code{"dynaNIRA_kim"} object from \code{\link{fitKIM}}.
#' @param intervention_type Character. \code{"alleviating"} (default) applies
#'   a negative threshold shift (reducing activation); \code{"aggravating"}
#'   applies a positive shift (increasing activation).
#' @param n_sim Positive integer. Number of Monte Carlo trajectories.
#'   Default \code{5000L}.
#' @param n_steps Positive integer. Simulation steps per trajectory. Should
#'   exceed \code{kim_fit$mixing_time}. Default \code{14L}.
#' @param intervention_sd Numeric. Intervention strength as a multiple of
#'   \code{SD(h)}. Default \code{2}. Passed to \code{\link{kim_compute_delta}}.
#' @param seed Integer random seed for reproducibility. Default \code{2025L}.
#' @param verbose Logical. Print progress messages. Default \code{TRUE}.
#'
#' @return A named list of class \code{"kim_nira_result"}:
#' \describe{
#'   \item{\code{summary}}{Data frame ranked by \code{|mean_diff|}: columns
#'     \code{node}, \code{mean_diff}, \code{cohen_d}, \code{cohen_d_se},
#'     \code{delta_p}, \code{rank}.}
#'   \item{\code{ranked_targets}}{Character vector of node names ordered from
#'     most to least effective intervention target.}
#'   \item{\code{baseline_scores}}{Numeric vector (length \code{n_sim}) of
#'     baseline projected activation rates.}
#'   \item{\code{intervention_scores}}{Named list of length \code{p}; each
#'     element is a numeric vector of intervention projected activation rates.}
#'   \item{\code{delta}}{The computed threshold shift \eqn{\delta}.}
#'   \item{\code{delta_prob_table}}{Data frame of direct marginal probability
#'     changes from \code{\link{kim_delta_prob_table}}.}
#'   \item{\code{params}}{List of call parameters.}
#' }
#'
#' @examples
#' \dontrun{
#' fit  <- fitKIM(data, vars = c("dep","anx","sleep"), id="id", time="time")
#' nira <- runKIMNIRA(fit, intervention_type = "alleviating")
#' print(nira)
#' }
#'
#' @seealso \code{\link{fitKIM}}, \code{\link{kimPermutationNIRAtest}},
#'   \code{\link{kimPlotTrajectory}}, \code{\link{kim_compute_delta}}
#' @importFrom stats sd
#' @export
runKIMNIRA <- function(kim_fit,
                       intervention_type = c("alleviating", "aggravating"),
                       n_sim             = 5000L,
                       n_steps           = 14L,
                       intervention_sd   = 2,
                       seed              = 2025L,
                       verbose           = TRUE) {

  if (!inherits(kim_fit, "dynaNIRA_kim"))
    stop("`kim_fit` must be produced by fitKIM() (class = 'dynaNIRA_kim').")

  intervention_type <- match.arg(intervention_type)

  J      <- kim_fit$J
  h      <- kim_fit$h
  vars   <- kim_fit$vars
  n_vars <- length(vars)

  # Mixing time warning: compare n_steps against the rho_star-based estimate.
  # mixing_time is derived from 1 / (1 - rho_star), the weighted-Jacobian
  # spectral radius that correctly accounts for sigmoid compression.
  mt <- kim_fit$mixing_time
  if (is.null(mt)) mt <- kim_fit$mixing_time_approx
  if (!is.null(mt) && is.finite(mt) && n_steps < mt)
    warning(
      "n_steps (", n_steps, ") is less than the approximate mixing time (",
      mt, ") derived from the weighted Jacobian spectral radius rho_star = ",
      round(kim_fit$spectral_radius_star %||% kim_fit$spectral_radius, 4), ".\n",
      "  Intervention effects may not have fully propagated. ",
      "Consider increasing n_steps to at least ", mt, "."
    )

  # Compute intervention strength
  delta <- kim_compute_delta(h, intervention_sd, intervention_type)
  dpt   <- kim_delta_prob_table(h, delta, vars)

  if (verbose) {
    message("-- KIM NodeIdentifyR Algorithm (NIRA) --")
    message("  sampling     = independent (nodes conditionally independent)")
    message("  paired       = TRUE")
    message("  intervention = ", intervention_type)
    message("  n_sim        = ", n_sim,
            "  |  n_steps = ", n_steps)
    message("  delta        = ", round(delta, 4),
            "  (", intervention_sd, " x SD(h) = ",
            round(sd(h), 4), ")")
    message("  mean |delta_p| = ", round(mean(abs(dpt$delta_p)), 4))
  }

  # Step 1: Sample initial states
  if (verbose) message("  [1/3] Sampling initial states (n = ", n_sim, ")...")
  init_states <- kim_init_states(
    h = h, n_sim = n_sim, vars = vars, seed = seed
  )

  # Step 2: Pre-generate shared random numbers
  if (verbose) message("  Pre-generating ", n_steps, " steps of shared random numbers (paired design)...")
  U_list <- .kim_pregenerate_U(
    n_sim = n_sim, n_vars = n_vars, n_steps = n_steps,
    seed = seed + 1000L
  )

  # Step 3: Baseline simulation
  if (verbose) message("  [2/3] Baseline simulation...")
  baseline_scores <- kim_trajectory(
    init_states = init_states, J = J, h = h,
    n_steps = n_steps,
    target_node = NULL, delta = 0, U_list = U_list
  )

  # Step 4: Intervention simulation
  if (verbose) message("  [3/3] Intervention simulation (", n_vars, " nodes)...")
  intervention_scores <- vector("list", n_vars)
  names(intervention_scores) <- vars

  for (j in seq_len(n_vars)) {
    if (verbose) message("    Node ", j, "/", n_vars, ": ", vars[j])
    intervention_scores[[j]] <- kim_trajectory(
      init_states = init_states, J = J, h = h,
      n_steps = n_steps,
      target_node = j, delta = delta, U_list = U_list
    )
  }

  # Compute effect sizes
  summary_df <- do.call(rbind, lapply(seq_len(n_vars), function(j) {
    sc          <- intervention_scores[[j]]
    diff_scores <- sc - baseline_scores
    md          <- mean(diff_scores)
    sd_diff     <- sd(diff_scores)
    d           <- if (sd_diff > 1e-10) md / sd_diff else 0
    se_d        <- sqrt((1 / n_sim) + (d^2 / (2 * n_sim)))
    data.frame(
      node       = vars[j],
      mean_diff  = round(md, 4),
      cohen_d    = round(d,  4),
      cohen_d_se = round(se_d, 4),
      delta_p    = round(dpt$delta_p[j], 4),
      stringsAsFactors = FALSE
    )
  }))

  summary_df      <- summary_df[order(abs(summary_df$mean_diff),
                                      decreasing = TRUE), ]
  summary_df$rank <- seq_len(n_vars)
  rownames(summary_df) <- NULL

  result <- list(
    summary             = summary_df,
    ranked_targets      = summary_df$node,
    baseline_scores     = baseline_scores,
    intervention_scores = intervention_scores,
    init_states         = init_states,
    delta_prob_table    = dpt,
    delta               = delta,
    paired              = TRUE,
    kim_fit             = kim_fit,
    params              = list(
      intervention_type  = intervention_type,
      n_sim              = n_sim,
      n_steps            = n_steps,
      intervention_sd    = intervention_sd,
      seed               = seed
    )
  )

  class(result) <- "kim_nira_result"

  if (verbose)
    message("  Done. Best intervention target: [", result$ranked_targets[1], "]")

  result
}


# -----------------------------------------------------------------------------
# S3 methods
# -----------------------------------------------------------------------------

#' @export
print.kim_nira_result <- function(x, ...) {
  p <- x$params
  cat("\n-- KIM NodeIdentifyR Algorithm (NIRA) Result --\n")
  cat("Intervention : ", p$intervention_type, "\n")
  cat("Nodes        : ", paste(x$kim_fit$vars, collapse = ", "), "\n")
  cat("n_sim =", p$n_sim, "  |  n_steps =", p$n_steps, "\n")
  cat("Sampling     : independent (nodes conditionally independent)\n")
  cat("Paired design: yes\n")
  cat("delta        =", round(x$delta, 4), "\n\n")
  cat("Probability change table (delta_p):\n")
  print(x$delta_prob_table[, c("node", "p_original",
                               "p_intervened", "delta_p")],
        row.names = FALSE)
  cat("\nTarget ranking (sorted by |mean_diff| descending):\n")
  cols <- c("rank", "node", "mean_diff", "cohen_d", "cohen_d_se", "delta_p")
  cols <- intersect(cols, names(x$summary))
  print(x$summary[, cols], row.names = FALSE)
  cat("\nBest intervention target: [", x$ranked_targets[1], "]\n")
  cat("\nNote: ranking criterion = |mean_diff| (change in mean activation rate after intervention).\n")
  cat("      delta_p   = direct marginal effect of intercept shift only (no network propagation).\n")
  cat("      mean_diff = total effect after network propagation, captured by Monte Carlo simulation.\n")
  cat("      cohen_d   = mean_diff / sd(diff); reflects effect stability.\n")
  cat("      cohen_d_se = Monte Carlo approximate SE of Cohen's d.\n")
  cat("\n")
  invisible(x)
}
