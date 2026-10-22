# =============================================================================
# kimDose.R -- intervention scale, dose-response and structural decomposition
#
# Two reviewer questions are answered here with exact calculations:
#
# (1) Does the ranking depend on how the intervention is defined? The default
#     NIRA intervention shifts every target's log-odds by the same delta. Because
#     the logistic link is nonlinear, an equal logit shift is not an equal change
#     in activation probability: a node whose activation probability is near 0.5
#     moves more than one near 0 or 1. kimDoseResponse() therefore (a) repeats
#     the ranking over a grid of logit doses and (b) re-ranks after calibrating
#     each node's dose so that its OWN stationary activation probability drops by
#     the same amount ("equal-probability" intervention).
#
# (2) Which structural feature makes a node a strong target? kimPathwayDecomposition()
#     recomputes each target's exact effect in knock-out versions of the network
#     (no outgoing edges, no edges feeding back into the target, no
#     self-persistence) and reports how much of the effect each feature carries.
#     These are model-based sensitivity contrasts, not causal mediation effects,
#     and they are not additive.
# =============================================================================


# Population-averaged stationary marginal of every node, optionally with node j
# shifted by d and with a transformed J.
.kim_pop_marginals <- function(fit, j = NULL, d = 0, Jfun = NULL) {
  pat <- .kim_patterns(fit); bs <- pat$baselines; p <- length(fit$vars)
  acc <- numeric(p)
  for (kk in seq_along(pat$index)) {
    k <- pat$index[kk]; Jk <- .kim_J_for(fit, k); b <- bs[k, ]
    if (!is.null(Jfun)) Jk <- Jfun(Jk)
    if (!is.null(j)) b[j] <- b[j] + d
    tr <- .kim_transition(Jk, b)
    acc <- acc + pat$weights[kk] * as.numeric(.kim_stationary(tr$P)$pi %*% tr$states)
  }
  stats::setNames(acc, fit$vars)
}

# Per-node doses giving an equal change in activation probability.
#   on = "stationary": the population-averaged STATIONARY activation probability
#        of the target itself changes by target_dp (solved by root finding on
#        the exact engine; feedback onto the target is included).
#   on = "intercept":  the conditional probability at the typical intercept,
#        sigma(h_j + delta_j) - sigma(h_j), equals +/- target_dp (closed form).
.kim_calibrate_delta <- function(fit, target_dp, intervention_type = "alleviating",
                                 on = c("stationary", "intercept"), max_abs = 10) {
  on <- match.arg(on)
  vars <- fit$vars; p <- length(vars)
  if (length(target_dp) != 1L || !is.finite(target_dp) || target_dp <= 0 || target_dp >= 1)
    stop("target_dp must be one number in (0, 1): the probability change per target.")
  sgn <- if (intervention_type == "alleviating") -1 else 1
  if (on == "intercept") {
    p0 <- stats::plogis(fit$h); p1 <- p0 + sgn * target_dp
    d <- ifelse(p1 > 0 & p1 < 1, stats::qlogis(p1) - fit$h, NA_real_)
  } else {
    .kim_check_exact_p(p)
    base <- .kim_pop_marginals(fit)
    d <- vapply(seq_len(p), function(j) {
      g <- function(x) .kim_pop_marginals(fit, j, x)[j] - base[j] - sgn * target_dp
      lim <- sgn * max_abs
      if (sign(g(lim)) == sign(g(0)) ) return(NA_real_)
      stats::uniroot(g, sort(c(0, lim)), tol = 1e-8)$root
    }, numeric(1))
  }
  d <- stats::setNames(as.numeric(d), vars)
  if (anyNA(d))
    stop("A probability change of ", target_dp, " is not attainable for node(s) ",
         paste(vars[is.na(d)], collapse = ", "), " (baseline probability too close to ",
         if (sgn < 0) "0" else "1", "). Choose a smaller target_dp.", call. = FALSE)
  attr(d, "calibration") <- list(target_dp = target_dp, on = on,
                                 baseline_probability = if (on == "intercept")
                                   stats::plogis(fit$h) else .kim_pop_marginals(fit))
  d
}


#' Dose-response and intervention-scale sensitivity of the NIRA ranking
#'
#' Re-runs the (exact by default) NIRA analysis over a grid of logit doses and,
#' optionally, under equal-probability calibration, and reports how the target
#' ranking changes.
#'
#' @param fit A `"dynaNIRA_kim"` fit.
#' @param deltas Signed logit doses (same sign as the intervention direction).
#'   Default: 0.25, 0.5, 1, 1.5, 2 and 3 logits.
#' @param target_dp Optional probability changes for equal-probability doses, e.g.
#'   `c(0.05, 0.10)`. Each node then receives its own logit dose such that its
#'   own activation probability changes by exactly that amount.
#' @param calibrate_on `"stationary"` (default; population-averaged stationary
#'   probability of the target) or `"intercept"` (conditional probability at the
#'   typical intercept; closed form).
#' @param intervention_type `"alleviating"` or `"aggravating"`.
#' @param method Engine passed to [runKIMNIRA()]; `"exact"` by default.
#' @param ... Further arguments for [runKIMNIRA()] (e.g. `n_sim` for `"mc"`).
#' @return A `"kim_dose_response"` list: `effects` (long table: dose, scale, node,
#'   effect, benefit, rank, logit dose used), `top1` (top-ranked node per dose),
#'   `rank_agreement` (Kendall tau between the rankings of every pair of doses),
#'   `top1_constant` and the per-dose `runs`.
#' @export
kimDoseResponse <- function(fit, deltas = NULL, target_dp = NULL,
                            calibrate_on = c("stationary", "intercept"),
                            intervention_type = c("alleviating", "aggravating"),
                            method = c("exact", "mc"), ...) {
  intervention_type <- match.arg(intervention_type)
  calibrate_on <- match.arg(calibrate_on)
  method <- match.arg(method)
  if (!inherits(fit, "dynaNIRA_kim")) stop("fit must come from fitKIM().")
  sgn <- if (intervention_type == "alleviating") -1 else 1
  if (is.null(deltas)) deltas <- sgn * c(0.25, 0.5, 1, 1.5, 2, 3)
  if (any(sign(deltas) != sgn)) stop("All deltas must have the sign of the intervention direction.")
  runs <- list(); rows <- list()
  add <- function(label, scale, r) {
    s <- r$summary
    rows[[label]] <<- data.frame(dose = label, scale = scale, node = s$node,
      effect = s$mean_diff, benefit = s$benefit, rank = s$rank,
      logit_delta = unname(r$delta[s$node]), stringsAsFactors = FALSE)
    runs[[label]] <<- r
  }
  for (d in deltas) add(sprintf("logit %+.2f", d), "logit",
    runKIMNIRA(fit, intervention_type = intervention_type, method = method,
               delta = d, verbose = FALSE, trajectory = FALSE, diagnose = FALSE, ...))
  for (dp in target_dp) add(sprintf("prob %.3f (%s)", dp, calibrate_on), "probability",
    runKIMNIRA(fit, intervention_type = intervention_type, method = method,
               dose_scale = "probability", target_dp = dp, calibrate_on = calibrate_on,
               verbose = FALSE, trajectory = FALSE, diagnose = FALSE, ...))
  eff <- do.call(rbind, rows); rownames(eff) <- NULL
  labels <- names(runs); vars <- fit$vars
  R <- vapply(labels, function(l) { e <- eff[eff$dose == l, ]; e$benefit[match(vars, e$node)] },
              numeric(length(vars)))
  tau <- matrix(NA_real_, length(labels), length(labels), dimnames = list(labels, labels))
  for (a in seq_along(labels)) for (b in seq_along(labels))
    tau[a, b] <- suppressWarnings(stats::cor(R[, a], R[, b], method = "kendall"))
  top1 <- data.frame(dose = labels, top1 = vapply(labels, function(l) {
    e <- eff[eff$dose == l, ]; paste(e$node[e$rank == 1], collapse = "/") }, ""),
    stringsAsFactors = FALSE)
  out <- list(effects = eff, top1 = top1, rank_agreement = tau,
              top1_constant = length(unique(top1$top1)) == 1L,
              runs = runs, intervention_type = intervention_type, method = method)
  class(out) <- "kim_dose_response"
  out
}

#' @export
print.kim_dose_response <- function(x, ...) {
  cat("\n-- NIRA dose-response / intervention-scale sensitivity --\n")
  print(x$top1, row.names = FALSE)
  cat("Top-ranked target constant across all doses:", x$top1_constant, "\n")
  cat("Minimum Kendall tau between any two rankings:",
      round(min(x$rank_agreement, na.rm = TRUE), 3), "\n")
  invisible(x)
}

#' Plot the dose-response of NIRA benefits
#' @param x A `"kim_dose_response"` object.
#' @return A ggplot object (logit doses on the x axis; equal-probability doses,
#'   if any, are shown as a separate facet of ranks).
#' @export
kimPlotDoseResponse <- function(x) {
  if (!inherits(x, "kim_dose_response")) stop("x must come from kimDoseResponse().")
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required.")
  e <- x$effects
  e$x <- ifelse(e$scale == "logit", abs(e$logit_delta), NA_real_)
  lg <- e[e$scale == "logit", , drop = FALSE]
  ggplot2::ggplot(lg, ggplot2::aes(x = .data$x, y = .data$rank, colour = .data$node)) +
    ggplot2::geom_line(linewidth = 0.9) + ggplot2::geom_point(size = 2) +
    ggplot2::scale_y_reverse(breaks = seq_len(length(unique(lg$node)))) +
    ggplot2::labs(x = "|delta| (logit units)", y = "Rank (1 = largest benefit)", colour = NULL) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "bottom", panel.grid.minor = ggplot2::element_blank())
}


#' Structural (knock-out) decomposition of NIRA effects
#'
#' For each target node j, the exact stationary effect of the same sustained
#' intervention is recomputed in modified networks, with BOTH the baseline and
#' the intervention condition evaluated in the modified network:
#'
#' * `isolated`: outgoing edges of j removed (J[-j, j] = 0). Nothing propagates
#'   from j, so the effect is the target's own change (including its
#'   self-persistence). `propagation = full - isolated` is the part of the
#'   effect carried by j's outgoing pathways.
#' * `no_feedback`: edges into j from other nodes removed (J[j, -j] = 0). Paths
#'   that leave j and return to it are cut. `feedback = full - no_feedback` is
#'   the part carried by cycles through j.
#' * `no_self`: the autoregressive edge of j removed (J[j, j] = 0).
#'   `self_persistence = full - no_self` is the amplification by the target's
#'   own persistence.
#'
#' The table also reports the outcome-node decomposition of the full effect
#' (`direct` = the target's own change / p, `spillover` = change in the other
#' nodes / p), the structural indices a centrality analysis would use
#' (out-strength and out-expected-influence excluding the diagonal, the
#' autoregressive weight), the target's baseline stationary activation
#' probability, and, optionally, the number of transitions until half of the
#' stationary effect is reached from the baseline stationary state.
#'
#' Knock-outs change the baseline as well as the response, and the contrasts are
#' not additive; they are sensitivity analyses that show which features of the
#' fitted network an effect depends on, not causal mediation estimates.
#'
#' @param x A `"kim_nira_result"` (its fit and dose are used) or a fit.
#' @param delta Dose (scalar or per node) when `x` is a fit.
#' @param timing Compute the half-effect time.
#' @param horizon Maximum number of transitions for the half-effect time.
#' @return A data frame with one row per target, ordered by the full effect, and
#'   attribute `rank_correlation` (Kendall tau between the NIRA benefit and
#'   out-strength / out-EI).
#' @export
kimPathwayDecomposition <- function(x, delta = NULL, timing = TRUE, horizon = 50L) {
  if (inherits(x, "kim_nira_result")) {
    fit <- x$kim_fit; if (is.null(delta)) delta <- x$delta
    type <- x$params$intervention_type
  } else if (inherits(x, "dynaNIRA_kim")) {
    fit <- x; type <- if (!is.null(delta) && all(delta < 0)) "alleviating" else "aggravating"
  } else stop("x must be a runKIMNIRA() result or a fitKIM() fit.")
  if (is.null(delta)) stop("Supply delta.")
  vars <- fit$vars; p <- length(vars)
  .kim_check_exact_p(p)
  delta <- .kim_expand_delta(delta, vars)
  sgn <- if (type == "alleviating") -1 else 1
  mu <- function(j, d, Jfun) mean(.kim_pop_marginals(fit, j, d, Jfun))
  effect_under <- function(j, Jfun) mu(j, delta[j], Jfun) - mu(NULL, 0, Jfun)

  full <- .kim_exact_effects(fit, delta)
  rows <- lapply(seq_len(p), function(j) {
    iso  <- effect_under(j, function(M) { M[-j, j] <- 0; M })
    nofb <- effect_under(j, function(M) { M[j, -j] <- 0; M })
    nos  <- effect_under(j, function(M) { M[j, j] <- 0; M })
    tot  <- unname(full$effect[j])
    data.frame(node = vars[j], delta = delta[j], effect_full = tot,
      direct = full$node_changes[j, j] / p,
      spillover = (sum(full$node_changes[j, ]) - full$node_changes[j, j]) / p,
      effect_isolated = iso, propagation = tot - iso,
      effect_no_feedback = nofb, feedback = tot - nofb,
      effect_no_self = nos, self_persistence = tot - nos,
      stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  Jo <- fit$J; diag(Jo) <- 0
  out$out_strength <- colSums(abs(Jo)); out$out_ei <- colSums(Jo)
  out$autoregressive <- diag(fit$J)
  out$baseline_probability <- unname(full$node_baseline)
  out$benefit <- sgn * out$effect_full
  if (timing) {
    resp <- .kim_response(fit, delta, horizon = horizon)$mean_changes  # (h+1) x p (targets)
    out$half_time <- vapply(seq_len(p), function(j) {
      path <- resp[, j]; target <- 0.5 * out$effect_full[j]
      hit <- which(abs(path) >= abs(target))
      if (length(hit)) as.integer(hit[1] - 1L) else NA_integer_
    }, integer(1))
  }
  out$rank_nira <- rank(-out$benefit, ties.method = "min")
  out$rank_out_strength <- rank(-out$out_strength, ties.method = "min")
  out <- out[order(out$rank_nira), , drop = FALSE]; rownames(out) <- NULL
  attr(out, "rank_correlation") <- c(
    out_strength = suppressWarnings(stats::cor(out$benefit, out$out_strength, method = "kendall")),
    out_ei = suppressWarnings(stats::cor(out$benefit, out$out_ei, method = "kendall")))
  attr(out, "note") <- paste("Knock-out contrasts are model-based sensitivity analyses;",
    "they are not additive and are not causal mediation effects.")
  out
}


#' Chance levels for target-recovery metrics in simulation studies
#'
#' Expected values of common recovery metrics when the estimated ranking is a
#' uniformly random permutation of p nodes.
#'
#' @param p Number of nodes.
#' @param k Size of the top set for overlap metrics.
#' @return A named numeric vector: probability that the true top-1 node is
#'   ranked first (1/p); probability that the true top-1 node is in the
#'   estimated top-k (k/p); expected proportion overlap of the true and
#'   estimated top-k sets (k/p); expected Spearman / Kendall correlation (0).
#' @export
kimChanceLevel <- function(p, k = 3L) {
  if (k > p) stop("k cannot exceed p.")
  c(top1_exact = 1 / p, top1_in_topk = k / p, topk_overlap = k / p,
    spearman = 0, kendall = 0)
}
