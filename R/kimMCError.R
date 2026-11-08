# =============================================================================
# kimMCError.R  --  MONTE-CARLO numerical uncertainty (method = "mc" ONLY)
#
# Sampling (person-bootstrap) error is handled by kimStabilityBootstrap().
# This function quantifies the OTHER source of uncertainty -- how much the
# finite-trajectory Monte-Carlo estimate wobbles across independent seeds --
# by repeating runKIMNIRA(method = "mc") n_runs times and plotting the
# mean effect across runs with across-run error bars. Not applicable to the
# exact engine (which has zero MC error).
# =============================================================================

#' Monte-Carlo numerical uncertainty of the MC engine (repeat seeds + plot)
#'
#' @param obj A `"kim_nira_result"` from [runKIMNIRA()] with `method = "mc"`.
#' @param n_runs Number of independent MC runs (different seeds). Default 20.
#' @param n_sim,n_steps Overrides for the MC runs (NULL = reuse `obj` settings).
#' @param seed Base seed; run r uses seed + (r-1)*7919.
#' @return list($effects (n_runs x p), $stat, $plot). `$stat$mc_sd` is the
#'   run-to-run standard deviation of one run; `run_range_lower/upper` (mean +/- 2
#'   SD) describe the spread of single-run results, NOT a confidence interval for
#'   the intervention effect. Sampling uncertainty: kimStabilityBootstrap().
#' @export
kimMCError <- function(obj, n_runs = 20L, n_sim = NULL, n_steps = NULL,
                       seed = 2025L) {
  if (!inherits(obj, "kim_nira_result"))
    stop("`obj` must be the output of runKIMNIRA().")
  if (is.null(obj$method) || obj$method != "mc")
    stop("kimMCError() quantifies MONTE-CARLO error and requires method = 'mc' ",
         "(the exact engine has zero MC error).")
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required.")
  fit    <- obj$kim_fit
  vars   <- fit$vars; p <- length(vars)
  type   <- obj$params$intervention_type
  aSD    <- obj$params$intervention_sd
  n_sim  <- if (is.null(n_sim)) obj$params$n_sim else n_sim
  n_steps<- if (is.null(n_steps)) obj$params$n_steps else n_steps
  n_runs <- as.integer(n_runs)
  if(n_runs<2)stop("n_runs must be at least two.")

  effM <- matrix(NA_real_, n_runs, p, dimnames = list(NULL, vars))
  for (r in seq_len(n_runs)) {
    rr <- tryCatch(runKIMNIRA(fit, intervention_type = type,
                              method = "mc", n_sim = n_sim, n_steps = n_steps,
                              intervention_sd = aSD, delta=obj$delta,
                              burn_in=obj$params$burn_in, paired=isTRUE(obj$paired),
                              trajectory = FALSE,   # only rr$summary is used here
                              diagnose   = FALSE,   # repeated runs are the diagnostic here
                              seed = seed + (r - 1L) * 7919L,
                              verbose = FALSE),
                   error = function(e) NULL)
    if (is.null(rr)) next
    md <- stats::setNames(rr$summary$mean_diff, rr$summary$node)
    effM[r, ] <- as.numeric(md[vars])
  }
  effM <- effM[rowSums(is.na(effM)) == 0, , drop = FALSE]
  if (nrow(effM) < 2) stop("Fewer than two valid MC runs.")
  if(nrow(effM)<n_runs)warning(n_runs-nrow(effM)," MC runs failed.")

  m  <- colMeans(effM)
  sd <- apply(effM, 2, stats::sd)
  # mc_sd is the run-to-run SD of ONE Monte Carlo run: the numerical precision a
  # single runKIMNIRA() call has. +/- 2 SD is the range of single-run results, not
  # a confidence interval for the effect (which needs the person bootstrap).
  st <- data.frame(node = vars, effect = m, mc_sd = sd, mc_se_mean = sd/sqrt(nrow(effM)),
                   run_range_lower = m - 2 * sd, run_range_upper = m + 2 * sd)
  st <- st[order(st$effect, decreasing = type=="aggravating"), ]
  st$node <- factor(st$node, levels = st$node)
  gg <- ggplot2::ggplot(st, ggplot2::aes(x = .data$node, y = .data$effect)) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2, colour = "grey40") +
    ggplot2::geom_col(fill = ifelse(st$effect >= 0, "#D95F02", "#1B9E77"),
                      width = 0.65, alpha = 0.85) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$run_range_lower,
                                        ymax = .data$run_range_upper),
                           width = 0.2) +
    ggplot2::labs(x = NULL, y = "effect across MC runs",
                  title = "Monte Carlo variability of a single run (mean across seeds +/- 2 run-to-run SD; not a CI)") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  list(effects = effM, stat = st, plot = gg, n_valid = nrow(effM))
}
