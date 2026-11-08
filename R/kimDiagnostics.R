# =============================================================================
# kimDiagnostics.R  ----  mixing, stationarity and Monte Carlo precision
#
# These helpers support the Monte Carlo engine's claim that it estimates the SAME
# stationary contrast as the exact engine:
#
#   Delta_j_stat = E_{pi_j}[ g(S) ] - E_{pi_0}[ g(S) ],   g(S) = mean_i S_i
#
# Three distinct questions are answered separately, and must not be conflated:
#
#   .kim_rhat()      is the warm-up long enough?  Chains started from dispersed
#                    states must agree; disagreement means warm-up bias, which is
#                    NOT visible in the Monte Carlo standard error.
#   .kim_ess_geyer() does the retained window hold enough independent information
#                    for its time average to have settled?
#   .kim_mc_se()     how precise is the reported effect?  The n_sim trajectories
#                    are independent of one another, so the across-trajectory
#                    standard error is valid and needs no ESS correction.
#
# Nothing here is exported.
# =============================================================================


#' Local relaxation time of the fitted chain
#'
#' Reuses the weighted-Jacobian spectral radius fitKIM() already reports (at the
#' mean-field fixed point, maximised over persons), so the heuristic is defined
#' in exactly one place. It is a linearisation; kimMixing() gives the exact
#' relaxation time -1/log|lambda_2| for p within the exact-engine limit. Falls back to computing it from
#' (J, h) when a hand-built fit does not carry it.
#'
#' @param fit A `"dynaNIRA_kim"` object.
#' @return Numeric relaxation time in steps, or `Inf` when the chain is not
#'   contractive (rho_star >= 1) and no finite time scale exists.
#' @keywords internal
.kim_relaxation <- function(fit) {
  mt <- fit$mixing_time
  if (!is.null(mt) && length(mt) == 1L && is.finite(mt)) return(as.numeric(mt))
  rho <- fit$spectral_radius_star_max %||% fit$spectral_radius_star
  if (is.null(rho) || length(rho) != 1L || !is.finite(rho)) {
    J <- fit$J; h <- fit$h
    s <- .kim_mf_fixed_point(J, h)
    rho <- max(abs(eigen(diag(s * (1 - s), length(s)) %*% J, only.values = TRUE)$values))
  }
  if (!is.finite(rho) || rho >= 1) return(Inf)
  max(1, 1 / (1 - rho))
}


#' Geyer initial-positive-sequence effective sample size
#'
#' Autocovariances are obtained by FFT, then summed in adjacent pairs
#' (lag 2k plus lag 2k+1) and truncated at the first non-positive pair. This is
#' the estimator used by Stan and coda::effectiveSize, and it is
#' deliberately conservative under strong autocorrelation, where naive
#' autocorrelation truncation overstates the effective sample size.
#'
#' A constant series has no defined effective sample size and returns `NA`: a
#' near-absorbing chain can sit in one state for a whole window, and reporting
#' "n independent draws" for a series that never moved would be a false comfort.
#'
#' @param x Numeric series from a single chain (not a set of independent draws).
#' @return Effective sample size in `[1, length(x)]`, or `NA_real_` when the
#'   series has (numerically) zero variance.
#' @keywords internal
.kim_ess_geyer <- function(x) {
  x <- as.numeric(x)
  n <- length(x)
  if (n < 2L) return(NA_real_)
  xc <- x - mean(x)
  g0 <- sum(xc * xc) / n
  if (!is.finite(g0) || g0 <= 0) return(NA_real_)
  if (n < 4L) return(as.numeric(n))

  nfft <- 2^ceiling(log2(2 * n))
  fx   <- stats::fft(c(xc, rep(0, nfft - n)))
  acov <- Re(stats::fft(fx * Conj(fx), inverse = TRUE))[seq_len(n)] / (nfft * n)

  # Pair adjacent lags starting at lag 1, so sum(Gamma) covers every lag >= 1.
  K <- floor((n - 1) / 2)
  if (K < 1L) return(as.numeric(n))
  Gamma <- acov[2 * seq_len(K)] + acov[2 * seq_len(K) + 1L]
  bad   <- which(Gamma <= 0)
  Kuse  <- if (length(bad)) min(bad) - 1L else K

  # tau = sigma^2 / gamma_0 = 1 + 2 * sum_{k>=1} rho_k, and sum(Gamma) is exactly
  # sum_{k>=1} gamma_k / gamma_0. Geyer truncation is mildly conservative: it can
  # overshoot tau, which understates ESS rather than overstating it.
  tau <- 1 + 2 * sum(Gamma[seq_len(Kuse)]) / g0
  if (!is.finite(tau) || tau < 1) return(as.numeric(n))
  min(as.numeric(n), n / tau)
}


#' Gelman-Rubin R-hat across dispersed-start chains
#'
#' @param chains List of numeric vectors, one per chain, holding the retained
#'   window means (or any scalar summary) from that chain.
#' @return R-hat, or `NA_real_` when there are too few chains or draws.
#' @keywords internal
.kim_rhat <- function(chains) {
  if (length(chains) < 2L) return(NA_real_)
  chains <- lapply(chains, as.numeric)
  n <- min(vapply(chains, length, numeric(1)))
  if (!is.finite(n) || n < 2L) return(NA_real_)
  x <- vapply(chains, function(z) z[seq_len(n)], numeric(n))   # n x m
  m <- ncol(x)
  if (m < 2L) return(NA_real_)
  w <- mean(apply(x, 2, stats::var))
  b <- n * stats::var(colMeans(x))
  if (!is.finite(w) || !is.finite(b)) return(NA_real_)
  if (w <= 0) return(if (b > 0) Inf else 1)
  sqrt(((n - 1) / n * w + b / n) / w)
}


#' Monte Carlo standard error of a trajectory-average
#'
#' The trajectories are independent of one another, so the across-trajectory
#' standard error is the right quantity here; ESS is NOT a correction for it.
#'
#' @param x Numeric vector of per-trajectory window means.
#' @return Standard error, or `NA_real_` for fewer than two values.
#' @keywords internal
.kim_mc_se <- function(x) {
  x <- as.numeric(x)
  if (length(x) < 2L) return(NA_real_)
  stats::sd(x) / sqrt(length(x))
}


#' Approximate two-sided person-bootstrap test
#'
#' The default compares absolute centred bootstrap errors with the observed
#' distance from the null. This is a bootstrap approximation, not an exact test;
#' regularization, selection, finite panels and MC bias can affect calibration.
#' Percentile mode is retained as an explicit alternative, with inclusive tails
#' so ties at the null cannot become significant. Neither method is guaranteed
#' to agree with BC intervals or multiplicity-adjusted decisions.
#' @param z Bootstrap replicates of paired intervention-minus-baseline effects.
#' @param point Original-sample effect; required for centred mode.
#' @param null Null value.
#' @param method centred or percentile.
#' @keywords internal
.kim_boot_p <- function(z, point = NULL, null = 0,
                        method = c("centered", "percentile")) {
  method <- match.arg(method)
  z <- as.numeric(z)
  if (length(null) != 1L || !is.finite(null)) stop("Invalid null value.")
  if (any(!is.finite(z))) stop("Bootstrap replicates must all be finite.")
  b <- length(z)
  if (b < 2L) return(NA_real_)
  if (method == "percentile")
    return(min(1, 2 * min((1 + sum(z <= null))/(b+1),
                          (1 + sum(z >= null))/(b+1))))
  if (length(point) != 1L || !is.finite(point))
    stop("Centred bootstrap testing requires the original point estimate.")
  observed <- abs(point - null)
  tol <- 32 * .Machine$double.eps * max(1, abs(point), abs(null), abs(z))
  if (observed <= tol) return(1)
  # A degenerate resampling distribution cannot estimate nonzero uncertainty.
  if (max(z)-min(z) <= tol) return(NA_real_)
  (1 + sum(abs(z-point) >= observed-tol))/(b+1)
}


#' Significance label from an adjusted p-value
#' @keywords internal
.kim_star_label <- function(p) {
  ifelse(is.na(p), "", ifelse(p < .001, "***", ifelse(p < .01, "**",
         ifelse(p < .05, "*", ""))))
}


# .kim_stationary_init() used to draw the Monte Carlo chain's initial states from
# the baseline stationary distribution wherever pi_0 could be solved for, which
# tied the Monte Carlo engine's start to the dense engine's feasibility. Chains
# now start from an independent Bernoulli draw at every p and rely on the warm-up
# diagnostics instead, so the helper had no callers left and was removed.
