# =============================================================================
# runKIMNIRA.R  ----  KIM NodeIdentifyR Algorithm (NIRA) Simulated Intervention Analysis
#
# This file can be sourced standalone without loading other files.
#
# Functions:
#   kim_sigmoid()              Sigmoid activation function
#   kim_compute_delta()        Compute intercept shift delta
#   kim_delta_prob_table()     Table of per-node probability changes
#   kim_init_states()          Sample initial states from independent distributions
#   kim_step()                 KIM single-step state update (independent sampling)
#   kim_trajectory()           Full trajectory simulation (optionally with burn-in)
#   runKIMNIRA()               Main function (paired MC; steady-state oriented)
#   print.kim_nira_result      Print method
#
# Revision notes (v2 / reviewer-driven):
#   * `burn_in` discards the initial transient so that the summary statistic is
#     oriented towards the chain's long-run (steady-state) activation, rather
#     than being contaminated by the chosen starting distribution.
#   * `scope` selects the target of the intervention simulation:
#       - "typical"  : one representative individual at the typical-person
#                      intercept h (the classic NIRA choice).
#       - "population": the population-average (marginal) intervention effect,
#                      i.e. the mean over the estimation-sample persons of each
#                      person's own effect (baseline h + person_shifts[k,], and
#                      person_J[k,,] under "res"). The exact engine computes it
#                      person by person; the Monte Carlo engine estimates the
#                      same finite-person mixture by assigning each trajectory
#                      one randomly drawn person. For "pooled" there are no
#                      person shifts, so population == typical.
#   * `delta` may be one value or one value per target. `dose_scale =
#     "probability"` calibrates per-target doses so that every target's own
#     activation probability changes by the same `target_dp` (see
#     kimDoseResponse()).
#   Paired (common-random-number) design is kept: the same initial states and
#   uniform stream are shared between baseline and every intervention condition
#   within a person.
#   `trajectory = TRUE` records the per-step path of those same chains, so the
#   stored path and the reported effects come from one run rather than two. The
#   recorded per-step means averaged over the retained window (burn_in, n_steps]
#   reproduce `effects`; see tests/regression_checks.R.
# =============================================================================


# -----------------------------------------------------------------------------
# Basic helpers
# -----------------------------------------------------------------------------

#' Sigmoid Function
#' @keywords internal
kim_sigmoid <- function(x) stats::plogis(x)


#' Compute Intervention Intercept Shift (delta)
#'
#' Follows Lunansky et al. (2022): delta = +/- intervention_sd * sd(h).
#' @param h Numeric vector of baseline log odds.
#' @param intervention_sd Nonnegative multiplier of the intercept standard deviation.
#' @param type Direction of the perturbation.
#' @export
kim_compute_delta <- function(h, intervention_sd = 2,
                              type = c("alleviating", "aggravating")) {
  type <- match.arg(type)
  h_sd <- sd(h)
  if (!is.finite(h_sd) || h_sd < 1e-6) {
    fallback <- median(abs(h))
    if (fallback < 1e-6) fallback <- 0.5
    h_sd <- fallback
    warning("Intercept vector SD is near zero. Using fallback scale ",
            round(h_sd, 4), ".")
  }
  delta <- intervention_sd * h_sd
  if (type == "alleviating") -delta else delta
}


#' Compute Per-Node Activation Probability Change Table
#' @description Conditional probabilities at the all-zero preceding state;
#'   these are not marginal stationary probabilities.
#' @param h Numeric vector of baseline log odds.
#' @param delta Signed additive log-odds change.
#' @param vars Node names.
#' @export
kim_delta_prob_table <- function(h, delta, vars) {
  data.frame(
    node         = vars,
    h_original   = round(h, 6),
    h_intervened = round(h + delta, 6),
    p_original   = round(kim_sigmoid(h), 6),
    p_intervened = round(kim_sigmoid(h + delta), 6),
    delta_p      = round(kim_sigmoid(h + delta) - kim_sigmoid(h), 6),
    stringsAsFactors = FALSE
  )
}

# Decompose the change in the mean activation score into the contribution from
# the intervention target itself and the change in all other outcome nodes.
.kim_node_effect_decomposition <- function(node_changes, intervention_type) {
  p <- ncol(node_changes)
  self_change <- diag(node_changes)
  other_change <- rowSums(node_changes) - self_change
  direction <- if (intervention_type == "alleviating") -1 else 1
  data.frame(
    direct_contribution = self_change / p,
    spillover_contribution = other_change / p,
    spillover_mean_diff = if (p > 1L) other_change / (p - 1L) else NA_real_,
    direct_benefit = direction * self_change / p,
    # Same denominator as direct_benefit, so the two add up to `benefit`
    # (the change in the mean activation score). The average change among
    # the p - 1 non-target nodes is reported separately as spillover_mean_diff.
    spillover_benefit = if (p > 1L) direction * other_change / p else NA_real_,
    stringsAsFactors = FALSE
  )
}


# -----------------------------------------------------------------------------
# Initial state sampling
# -----------------------------------------------------------------------------

#' Sample Initial States from Independent Distributions
#'
#' Each node is sampled independently: P(s_i = 1) = sigma(h_i).
#' @param h Numeric vector of baseline log odds.
#' @param n_sim Number of trajectories.
#' @param vars Optional node names.
#' @param seed Optional random seed.
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
# KIM single-step update
# -----------------------------------------------------------------------------

#' KIM Single-Step State Update (independent sampling)
#'
#' P(s_i(t+1) = 1 | s(t)) = sigma( sum_j J_ij * s_j(t) + h_i )
#' @param s Matrix of current states, one trajectory per row.
#' @param J Directed coefficients, with entry i,j denoting j to i. Either one
#'   `p x p` matrix shared by every trajectory, or an `n_sim x p x p` array
#'   giving each trajectory its own transition matrix (random-slope models).
#' @param h_eff Baseline vector or a trajectory-by-node baseline matrix.
#' @param U_step Optional matrix of uniform random numbers matching s.
#' @export
kim_step <- function(s, J, h_eff, U_step = NULL) {
  n_sim  <- nrow(s)
  n_vars <- ncol(s)
  h_mat  <- if(is.matrix(h_eff)) h_eff else matrix(h_eff, nrow=n_sim,ncol=n_vars,byrow=TRUE)
  if(!identical(dim(h_mat),dim(s)))stop("Person baselines do not match trajectory dimensions.")
  eta    <- if (length(dim(J)) == 3L) {
    if (!identical(as.integer(dim(J)), as.integer(c(n_sim, n_vars, n_vars))))
      stop("A trajectory-specific J must have dimensions n_sim x p x p.")
    # eta[k,i] = sum_j J[k,i,j] * s[k,j]. R does not broadcast array dimensions
    # in `*`, so each destination node is summed explicitly.
    out <- h_mat
    for (i in seq_len(n_vars)) {
      Ji <- matrix(J[, i, , drop = TRUE], nrow = n_sim, ncol = n_vars)
      out[, i] <- out[, i] + rowSums(Ji * s)
    }
    out
  } else s %*% t(J) + h_mat
  probs  <- kim_sigmoid(eta)
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
# Trajectory simulation (with burn-in)
# -----------------------------------------------------------------------------

#' Simulate a Full KIM Trajectory
#'
#' Runs a forward simulation of the Kinetic Ising Model for `n_steps` time
#' steps, discards the first `burn_in` steps, and returns the mean projected
#' activation rate of the RETAINED steps for each simulated trajectory.
#'
#' @param init_states Integer matrix (n_sim x p) of initial binary states.
#' @param J Numeric (p x p) temporal edge matrix (J[i,j] = effect of j on i).
#' @param h Numeric vector (length p) of node thresholds.
#' @param n_steps Number of time steps to simulate (total, incl. burn-in).
#' @param burn_in Number of leading steps discarded (>=0, <= n_steps-1).
#' @param target_node Integer or NULL (NULL = baseline).
#' @param delta Threshold shift applied to `target_node`.
#' @param U_list Pre-generated uniform matrices (one per step).
#'
#' @return Numeric vector of length n_sim.
#' @export
kim_trajectory <- function(init_states, J, h, n_steps,
                           burn_in        = 0L,
                           target_node    = NULL,
                           delta          = 0,
                           U_list         = NULL) {
  n_sim  <- nrow(init_states)
  n_vars <- ncol(init_states)
  burn_in <- max(0L, as.integer(burn_in))
  if (burn_in >= n_steps) stop("burn_in must be < n_steps.")

  h_eff <- h
  if (!is.null(target_node) && delta != 0) {
    if(is.matrix(h_eff)) h_eff[,target_node]<-h_eff[,target_node]+delta
    else h_eff[target_node] <- h[target_node] + delta
  }

  s         <- init_states
  cum_score <- numeric(n_sim)

  for (step in seq_len(n_steps)) {
    U_step <- if (!is.null(U_list)) U_list[[step]] else NULL
    s      <- kim_step(s, J, h_eff, U_step = U_step)
    if (step > burn_in)
      cum_score <- cum_score + rowSums(s)
  }
  cum_score / ((n_steps - burn_in) * n_vars)
}


# -----------------------------------------------------------------------------
# Main function
# -----------------------------------------------------------------------------

#' Run the NodeIdentifyR Algorithm (NIRA) on a Fitted KIM
#'
#' @param kim_fit A `"dynaNIRA_kim"` object from [fitKIM()].
#' @param intervention_type `"alleviating"` (default) or `"aggravating"`.
#' @param method Monte Carlo or dense exact finite-state calculation. The exact
#'   engine is dense over `2^p` states, so its cost grows roughly as `2^(3p)` and
#'   memory as `4^p`; it fails on the allocation rather than at a fixed node count.
#'   Use `"mc"` where that is not affordable.
#' @param delta Fixed signed log-odds shift: one value for every target, or one
#'   value per target (named or in node order). Prefer specifying it explicitly;
#'   the default `intervention_sd x SD(h)` depends on the spread of the fitted
#'   intercepts and therefore changes with `between_person`.
#' @param dose_scale `"logit"` (default): every target receives the same logit
#'   shift. `"probability"`: each target receives its own logit shift, chosen so
#'   that its own activation probability changes by `target_dp`; this removes the
#'   dependence of the ranking on where each node sits on the logistic curve.
#' @param target_dp Probability change per target for `dose_scale = "probability"`.
#' @param calibrate_on `"stationary"` (population-averaged stationary probability of
#'   the target, solved exactly; requires p within the exact-engine limit) or
#'   `"intercept"` (conditional probability at the typical intercept; closed form).
#' @param paired Use common innovation random numbers in Monte Carlo.
#' @param n_sim Number of Monte Carlo trajectories. Default 5000.
#' @param n_steps Total steps simulated per trajectory (should exceed the
#'   mixing time; default 30). Kept for compatibility: it equals
#'   `warmup + sampling_steps` and is derived from those when either is supplied.
#' @param burn_in Leading steps discarded before the summary statistic is
#'   computed. NULL discards half the requested steps; validate burn-in separately.
#'   Kept for compatibility with `warmup`; see there.
#' @param warmup,sampling_steps Explicit split of `n_steps` into the discarded
#'   warm-up prefix and the retained sampling window. Supply either to derive the
#'   other (`n_steps` stays `warmup + sampling_steps` and `burn_in` becomes
#'   `warmup`). Left `NULL`, the existing `n_steps`/`burn_in` pair is used
#'   unchanged, so defaults are unaffected. `sampling_steps` is what the effects
#'   are averaged over, and it is the quantity compared against the chain's
#'   relaxation time.
#' @param diagnose Run Monte Carlo diagnostics (`method = "mc"` only) and store
#'   them in `$mc_diagnostics`: per-node Monte Carlo standard errors, the Geyer
#'   initial-positive-sequence effective sample size of the retained window, and
#'   a multi-start R-hat. Flags warm-up inadequacy with a warning rather than
#'   silently returning a biased estimate. Set `FALSE` in inner loops (the person
#'   bootstrap and [kimMCError()] do) where the extra pilot chains are not worth
#'   the cost. Note that `$diagnostics` is a different field: the exact engine
#'   uses it for stationary-solver residuals.
#' @param intervention_sd Intervention strength (multiple of SD(h)). Default 2.
#' @param seed Integer random seed.
#' @param verbose Logical.
#' @param trajectory Record the per-step path of the same Monte Carlo chains that
#'   produce `effects` and store it as `$trajectory`, so that
#'   [kimPlotTrajectory()] plots the run it describes instead of a second one.
#'   Set `FALSE` when only the effects are needed (the person bootstrap does).
#'   Ignored for `method = "exact"`: exact paths are deterministic propagation
#'   and [kimPlotTrajectory()] recomputes them on demand without any drift.
#'
#' @return A `"kim_nira_result"` list (summary, ranked_targets,
#'   baseline_scores, intervention_scores, delta, ...). With
#'   `trajectory = TRUE` and `method = "mc"`, `$trajectory` also holds the
#'   recorded path (`$baseline`, `$interventions`, `$n_steps`, `$burn_in`).
#'   Its per-step means averaged over the retained window `(burn_in, n_steps]`
#'   reproduce `$effects`; see `tests/regression_checks.R`. With
#'   `diagnose = TRUE` and `method = "mc"`, `$mc_diagnostics` holds `mc_se`,
#'   `ess`, `rhat`, `relaxation`, `initialization` and any `flags`.
#'
#'   The result also contains `node_changes` (target x outcome-node changes),
#'   `direct_effects` (the target node's contribution to the original mean score),
#'   and `spillover_effects` (the mean change among the other nodes). The summary
#'   separates `direct_contribution` and `spillover_contribution`, which add to
#'   `mean_diff`, from `spillover_mean_diff`, which averages over the p - 1
#'   non-target nodes. `direct_benefit` and `spillover_benefit` both divide by p,
#'   so they add up to `benefit`; `spillover_mean_diff` remains the per-other-node
#'   average. These are outcome-node decompositions, not path-specific causal
#'   mediation effects. The total-effect ranking remains in `rank`;
#'   `spillover_rank` ranks the effects after excluding the target node.
#'
#'   Two quantities are deliberately kept distinct. `$effects` estimates
#'   `Delta_j_stat`, the change in the stationary mean of g(S) under a sustained
#'   perturbation of node j -- the same contrast the exact engine computes.
#'   `$trajectory` describes `Delta_j(t)`, the response
#'   path starting from the baseline stationary state; its average over the
#'   retained window reproduces `$effects`, but individual steps do not.
#'
#' @details
#' Under `method = "mc"` the chains start from an independent Bernoulli draw at
#' the baseline probabilities and must reach the stationary law inside `warmup`.
#' Warm-up adequacy is therefore checked rather than assumed.
#' `$mc_diagnostics$rhat` compares chains started from dispersed states (all
#' inactive, all active, uniform) with the SAME persons assigned in every start,
#' so it isolates warm-up bias, which the Monte Carlo standard error cannot
#' cover. `$mc_diagnostics$relaxation` is the linearised relaxation time at the
#' slowest person; the exact value is available from [kimMixing()]. A warning is
#' issued when the retained window is shorter than five relaxation times, when
#' R-hat exceeds 1.05, or when the Geyer ESS indicates that the window holds
#' fewer than five autocorrelation times. Chains that never move inside the
#' window have no defined ESS and are excluded from that check.
#'
#' One step is one transition of the fitted model, i.e. one observation interval
#' of the data used in [fitKIM()] (see `$kim_fit$transition_types`).
#' @export
runKIMNIRA <- function(kim_fit,
                       intervention_type = c("alleviating", "aggravating"),
                       method            = c("mc", "exact"),
                       n_sim             = 5000L,
                       n_steps           = 30L,
                       burn_in           = NULL,
                       intervention_sd   = 2,
                       seed              = 2025L,
                       verbose           = TRUE,
                       delta             = NULL,
                       paired            = TRUE,
                       trajectory        = TRUE,
                       warmup            = NULL,
                       sampling_steps    = NULL,
                       diagnose          = TRUE,
                       dose_scale        = c("logit", "probability"),
                       target_dp         = NULL,
                       calibrate_on      = c("stationary", "intercept")) {

  if (!inherits(kim_fit, "dynaNIRA_kim"))
    stop("`kim_fit` must be produced by fitKIM() (class = 'dynaNIRA_kim').")

  intervention_type <- match.arg(intervention_type)
  method            <- match.arg(method)
  dose_scale        <- match.arg(dose_scale)
  calibrate_on      <- match.arg(calibrate_on)

  J      <- kim_fit$J
  h      <- kim_fit$h
  vars   <- kim_fit$vars
  n_vars <- length(vars)
  shifts <- kim_fit$person_shifts                 # n_persons x p (or NULL)

  # ---- warm-up / sampling split (steady-state orientation) -------------------
  # `n_steps`/`burn_in` remain the documented interface; `warmup`/`sampling_steps`
  # expose the split explicitly. Either pair derives the other, so a call that
  # supplies neither keeps exactly its previous meaning.
  if (!is.null(warmup) || !is.null(sampling_steps)) {
    if (is.null(warmup))         warmup         <- floor(n_steps/2)
    if (is.null(sampling_steps)) sampling_steps <- n_steps - warmup
    if (!is.null(burn_in) && isTRUE(burn_in != warmup))
      warning("`burn_in` = ", burn_in, " ignored: `warmup` = ", warmup,
              " takes precedence.", call. = FALSE)
    n_steps <- warmup + sampling_steps
    burn_in <- warmup
  } else {
    if (is.null(burn_in)) burn_in <- floor(n_steps/2)
    warmup         <- burn_in
    sampling_steps <- n_steps - burn_in
  }
  bad_split <- function(x) length(x)!=1 || !is.finite(x) || x != as.integer(x)
  if (method == "mc" && (bad_split(n_steps) || n_steps < 2 || bad_split(n_sim) || n_sim < 2))
    stop("Monte Carlo requires integer n_steps >= 2 and n_sim >= 2.")
  if (method == "mc" && (bad_split(warmup) || warmup < 0 || warmup >= n_steps ||
                         bad_split(sampling_steps) || sampling_steps < 1))
    stop("Require integer warmup in [0, n_steps-1] and integer sampling_steps >= 1.")

  # ---- which persons to simulate over ----------------------------------------
  # Population-average (marginal) effect over the estimation-sample persons when
  # person baselines exist (between_person = re/res); otherwise a typical-person
  # effect (between_person = pooled has zero shifts -> identical to typical).
  # Shared with the exact engine so the two cannot disagree about the estimand,
  # and TRUE even when every intercept shift is zero if persons have their own J.
  use_pop <- .kim_uses_population(kim_fit)
  if (use_pop && verbose)
    message("  population-average effect over ", nrow(shifts), " persons")

  # Intervention strength. Logit scale: one shift for every target (default
  # intervention_sd x SD(h), as in Lunansky et al., 2022) or a user-supplied
  # per-target vector. Probability scale: per-target shifts calibrated so that
  # every target's own activation probability changes by target_dp.
  if (dose_scale == "probability") {
    if (is.null(target_dp)) stop("dose_scale = 'probability' requires target_dp.")
    if (!is.null(delta)) warning("`delta` ignored: dose_scale = 'probability' calibrates it.",
                                 call. = FALSE)
    delta <- .kim_calibrate_delta(kim_fit, target_dp, intervention_type, on = calibrate_on)
  }
  if (is.null(delta)) delta <- kim_compute_delta(h, intervention_sd, intervention_type)
  calib <- attr(delta, "calibration")
  delta <- .kim_expand_delta(delta, vars)
  if ((intervention_type=="alleviating" && any(delta>0)) ||
      (intervention_type=="aggravating" && any(delta<0))) stop("delta sign conflicts with intervention_type.")
  common_dose <- length(unique(delta)) == 1L
  dpt   <- kim_delta_prob_table(h, delta, vars)

  if (verbose) {
    message("-- KIM NodeIdentifyR Algorithm (NIRA) --")
    message("  method       = ", method)
    message("  effect       = ", if (use_pop) "population-average (marginal)" else "typical-person")
    message("  intervention = ", intervention_type)
    if (method == "mc") {
      message("  paired       = ",paired)
      message("  n_sim = ", n_sim, " | n_steps = ", n_steps, " | burn_in = ", burn_in)
    }
    if (dose_scale == "probability")
      message("  dose         = equal-probability (", calibrate_on, ", dp = ", target_dp,
              "); logit doses ", paste(round(delta, 3), collapse = ", "))
    else if (common_dose)
      message("  delta = ", round(delta[1], 4), " logits (", intervention_sd, " x SD(h) = ",
              round(sd(h), 4), " if not supplied)")
    else message("  delta (per target) = ", paste(round(delta, 3), collapse = ", "))
  }

  # ---- EXACT engine: zero-MC-error steady-state effects -----------------------
  if (method == "exact") {
    pt  <- .kim_exact_effects(kim_fit, delta)
    eff <- pt$effect
    summary_df <- data.frame(
      node       = vars,
      mean_diff  = unname(eff),
      delta_p    = dpt$delta_p,
      stringsAsFactors = FALSE)
    node_changes <- pt$node_changes
    dimnames(node_changes) <- list(vars, vars)
    decomposition <- .kim_node_effect_decomposition(node_changes, intervention_type)
    summary_df <- cbind(summary_df, decomposition)
    summary_df$benefit <- if(intervention_type=="alleviating") -summary_df$mean_diff else summary_df$mean_diff
    summary_df$rank <- rank(-summary_df$benefit,ties.method="min")
    summary_df$spillover_rank <- rank(-summary_df$spillover_benefit,ties.method="min",na.last="keep")
    summary_df <- summary_df[order(-summary_df$benefit), ]
    rownames(summary_df) <- NULL
    result <- list(
      summary             = summary_df,
      ranked_targets      = summary_df$node,
      effects             = eff,
      baseline_scores     = NULL,
      intervention_scores = setNames(vector("list", n_vars), vars),
      trajectory          = NULL,   # deterministic; computed on demand when plotting
      delta               = delta,
      delta_prob_table    = dpt,
      paired              = NA,      # not applicable: no simulation
      method              = "exact",
      population_averaged = pt$population_averaged,
      base_mean           = pt$base,
      diagnostics         = pt$diagnostics,
      mc_diagnostics      = NULL,   # Monte Carlo only; see $mc_diagnostics below
      node_changes        = node_changes,
      direct_effects      = stats::setNames(decomposition$direct_contribution, vars),
      spillover_effects   = stats::setNames(decomposition$spillover_mean_diff, vars),
      ranked_spillover_targets = summary_df$node[order(summary_df$spillover_rank)],
      node_baseline       = pt$node_baseline,
      exact_conditions    = data.frame(condition = c("original", vars),
                                       mean_activation = c(pt$base, pt$base + eff)),
      kim_fit             = kim_fit,
      params = list(intervention_type = intervention_type, method = "exact",
                    n_sim = NA_integer_, n_steps = NA_integer_, burn_in = NA_integer_,
                    intervention_sd = intervention_sd, delta = delta, paired = NA,
                    dose_scale = dose_scale, target_dp = target_dp,
                    calibrate_on = calibrate_on, calibration = calib,
                    seed = seed))
    class(result) <- "kim_nira_result"
    if (verbose) {
      if (isTRUE(trajectory))
        message("  Note: exact paths are deterministic; kimPlotTrajectory() computes ",
                "the stepwise path on demand, so none is stored here.")
      message("  Done. Top-ranked target (exact; largest estimated benefit): [", result$ranked_targets[1], "]",
              "  (effect = ", eff[result$ranked_targets[1]], ")")
    }
    return(result)
  }

  # ---- MC engine --------------------------------------------------------------
  set.seed(seed)
  U_list <- if(paired) .kim_pregenerate_U(n_sim, n_vars, n_steps, seed = seed + 1000L) else NULL

  # Each independent trajectory samples one empirical person baseline.
  # This estimates the same finite-person mixture without N x M trajectories.
  person_indices <- if(use_pop) sample.int(nrow(shifts),n_sim,replace=TRUE) else rep(1L,n_sim)
  h_k <- if(use_pop) sweep(shifts[person_indices,,drop=FALSE],2,h,"+") else
    matrix(h,n_sim,n_vars,byrow=TRUE)
  # One transition matrix per trajectory when persons have their own slopes; a
  # shared J keeps the original single-matrix path unchanged.
  J_k <- if (is.null(kim_fit$person_J)) J else {
    arr <- kim_fit$person_J[person_indices, , , drop = FALSE]
    array(arr, dim = c(n_sim, n_vars, n_vars), dimnames = list(NULL, vars, vars))
  }
  # Initial states: an independent Bernoulli draw at the baseline probabilities.
  # The chain therefore has to reach the stationary law inside `warmup` at every p,
  # and warm-up adequacy is not assumed -- `$mc_diagnostics$rhat` compares chains
  # launched from dispersed states and `$relaxation` compares the retained window
  # with the local mixing time. Starting each chain from the stationary
  # distribution instead would remove the warm-up problem only where pi_0 can be
  # solved, which is exactly where the exact engine is already available; the
  # diagnostics cover the whole range instead of a subset of it.
  init      <- matrix(as.integer(matrix(runif(n_sim*n_vars),n_sim,n_vars)<stats::plogis(h_k)),n_sim,n_vars)
  colnames(init) <- vars
  init_kind <- "independent Bernoulli at plogis(h_k)"

  # One pass per condition records the whole path AND the retained-window mean, so
  # the stored trajectory and `effects` come from the same chains. `window_mean`
  # reproduces what kim_trajectory() would have returned for the same window, hence
  # mean_t(Delta_j(t)) over the window equals effects[j] by construction.
  keep_traj <- isTRUE(trajectory)
  do_diag   <- isTRUE(diagnose)
  n_traj    <- if(keep_traj || do_diag) min(n_sim, 200L) else 0L
  run_chain <- function(target, dl) .kim_trajectory_stepwise(
    init_states = init, J = J_k, h = h_k, n_steps = n_steps,
    target_node = target, delta = dl, U_list = U_list,
    n_traces = n_traj, burn_in = burn_in)

  base_traj <- run_chain(NULL, 0)
  int_traj  <- lapply(seq_len(n_vars), function(j) run_chain(j, delta[j]))
  names(int_traj) <- vars

  retained_rows <- seq.int(burn_in + 2L, n_steps + 1L)
  node_changes <- t(vapply(seq_len(n_vars), function(j) {
    colMeans(int_traj[[j]]$node_probs[retained_rows, , drop = FALSE] -
               base_traj$node_probs[retained_rows, , drop = FALSE])
  }, numeric(n_vars)))
  dimnames(node_changes) <- list(vars, vars)
  decomposition <- .kim_node_effect_decomposition(node_changes, intervention_type)

  baseline_scores     <- base_traj$window_mean
  intervention_scores <- setNames(lapply(int_traj, `[[`, "window_mean"), vars)

  summary_df <- do.call(rbind, lapply(seq_len(n_vars), function(j) {
    sc          <- intervention_scores[[j]]
    diff_scores <- sc - baseline_scores
    md          <- mean(diff_scores)
    data.frame(
      node       = vars[j],
      mean_diff  = md,
      delta_p    = dpt$delta_p[j],
      stringsAsFactors = FALSE
    )
  }))

  summary_df$benefit <- if(intervention_type=="alleviating") -summary_df$mean_diff else summary_df$mean_diff
  summary_df$rank <- rank(-summary_df$benefit,ties.method="min")
  summary_df <- cbind(summary_df, decomposition)
  summary_df$spillover_rank <- rank(-summary_df$spillover_benefit,ties.method="min",na.last="keep")
  summary_df <- summary_df[order(-summary_df$benefit), ]
  rownames(summary_df) <- NULL

  # ---- diagnostics: is the Monte Carlo estimate trustworthy? ------------------
  # Three distinct questions, deliberately not conflated: R-hat asks whether the
  # warm-up was long enough, ESS whether the retained window holds enough
  # independent information, and mc_se how precise the reported effect is. Only
  # the last is a standard error; warm-up bias is invisible to it.
  diag <- NULL
  if (do_diag) {
    relax <- .kim_relaxation(kim_fit)
    mc_se <- stats::setNames(vapply(seq_len(n_vars), function(j)
      .kim_mc_se(intervention_scores[[j]] - baseline_scores), numeric(1)), vars)

    # ESS is per chain and cannot exceed `sampling_steps`, so `tau = sampling_steps
    # / ess` is the integrated autocorrelation time. Trajectories that sit in one
    # state for the whole window have no defined ESS and are excluded; a large
    # `ess_degenerate` share means the score barely moved and tau is not
    # interpretable -- which happens for near-absorbing chains, where the estimate
    # is precise anyway. ESS is therefore reported, not gated on.
    ess <- NA_real_
    tau <- NA_real_
    ess_degenerate <- NA_real_
    if (!is.null(base_traj$traces) && ncol(base_traj$traces) > 0L) {
      win <- base_traj$traces[seq.int(burn_in + 2L, n_steps + 1L), , drop = FALSE]
      ev  <- apply(win, 2, .kim_ess_geyer)
      ess_degenerate <- mean(is.na(ev))
      if (any(!is.na(ev))) {
        ess <- mean(ev, na.rm = TRUE)
        tau <- sampling_steps / ess
      }
    }

    # Multi-start agreement. Each pilot chain gets its OWN innovation stream: with
    # a shared stream the chains synchronise after a few steps and then agree by
    # construction, which makes R-hat blind exactly where it is needed. Independent
    # noise keeps the between-chain variance meaningful.
    #
    # This runs at every p. The pilots are ordinary trajectories and draw nothing
    # from the dense engine, so the check used to be restricted to small p only
    # because the stationary initialisation was, and that is gone.
    # The SAME persons are assigned in every start. Otherwise between-start
    # differences would mix warm-up bias with person-sampling noise and R-hat
    # would be diluted exactly where it is needed.
    pilot <- min(n_sim, 500L)
    set.seed(seed + 7919L)
    starts <- list(rep(0L, n_vars), rep(1L, n_vars),
                   as.integer(runif(pilot * n_vars) < 0.5))
    pid <- if (use_pop) sample.int(nrow(shifts), pilot, replace = TRUE) else rep(1L, pilot)
    hh  <- if (use_pop) sweep(shifts[pid,,drop=FALSE], 2, h, "+") else
      matrix(h, pilot, n_vars, byrow = TRUE)
    # The pilot must use the same person-specific transitions as the main
    # chains, otherwise R-hat compares two different models.
    JJ  <- if (is.null(kim_fit$person_J)) J else
      array(kim_fit$person_J[pid, , , drop = FALSE],
            dim = c(pilot, n_vars, n_vars), dimnames = list(NULL, vars, vars))
    run_pilot <- function(start_states, s) {
      st <- if (is.matrix(start_states)) start_states else
        matrix(start_states, pilot, n_vars, byrow = TRUE)
      UU  <- if (paired) .kim_pregenerate_U(pilot, n_vars, n_steps,
                                            seed = s + 1000L) else NULL
      .kim_trajectory_stepwise(st, JJ, hh, n_steps, target_node = NULL, delta = 0,
                               U_list = UU, n_traces = 0L,
                               burn_in = burn_in)$window_mean
    }
    rhat <- .kim_rhat(lapply(seq_along(starts),
                            function(i) run_pilot(starts[[i]], seed + 7919L + 101L * i)))

    flags <- character(0)
    if (!is.finite(relax)) {
      flags <- c(flags, "weighted-Jacobian spectral radius >= 1; no finite relaxation time")
    } else if (sampling_steps < 5 * relax) {
      flags <- c(flags, sprintf("retained window %d steps < 5 x relaxation time %.1f",
                                sampling_steps, relax))
    }
    if (!is.na(rhat) && rhat > 1.05)
      flags <- c(flags, sprintf("multi-start R-hat %.3f > 1.05", rhat))
    if (is.finite(ess) && ess * 5 < sampling_steps)
      flags <- c(flags, sprintf(paste0("effective sample size %.1f per chain spans only ",
                                       "%.1f x its own window of %d steps; the window ",
                                       "average may not have settled"),
                                ess, tau, sampling_steps))
    if (length(flags))
      warning("Monte Carlo diagnostics flagged: ", paste(flags, collapse="; "),
              ". Warm-up bias is not covered by the Monte Carlo standard error; ",
              "raise `warmup`/`sampling_steps`.", call. = FALSE)

    diag <- list(mc_se = mc_se, ess = ess, tau = tau, ess_degenerate = ess_degenerate,
                 rhat = rhat,
                 rho_star = kim_fit$spectral_radius_star, relaxation = relax,
                 warmup = warmup, sampling_steps = sampling_steps,
                 initialization = init_kind, flags = flags,
                 note = paste("rhat detects gross warm-up failure;",
                              "relaxation compares the retained window with the",
                              "local mixing time; mc_se covers Monte Carlo noise only."))
  }

  # Shape the recorded path the way kimPlotTrajectory() consumes it.
  traj <- if(keep_traj) {
    keep_names <- c("node_probs","summary","q025","q975","traces","window_mean",
                    "burn_in","n_retained")
    list(baseline       = base_traj[keep_names],
         interventions  = lapply(int_traj, function(z) z[keep_names]),
         vars           = vars,
         n_steps        = n_steps,
         burn_in        = burn_in,
         delta          = delta,
         method         = "mc",
         initialization = init_kind)
  } else NULL

  result <- list(
    summary             = summary_df,
    ranked_targets      = summary_df$node,
    effects             = stats::setNames(summary_df$mean_diff, summary_df$node)[vars],
    node_changes        = node_changes,
    direct_effects      = stats::setNames(decomposition$direct_contribution, vars),
    spillover_effects   = stats::setNames(decomposition$spillover_mean_diff, vars),
    ranked_spillover_targets = summary_df$node[order(summary_df$spillover_rank)],
    baseline_scores     = baseline_scores,
    intervention_scores = intervention_scores,
    trajectory          = traj,
    mc_diagnostics      = diag,   # exact keeps its stationary-solver info in $diagnostics

    init_states         = init,
    person_indices      = person_indices,
    delta_prob_table    = dpt,
    delta               = delta,
    paired              = paired,
    base_mean           = mean(baseline_scores),
    method              = "mc",
    population_averaged = use_pop,
    scope               = if (use_pop) "population" else "typical",
    kim_fit             = kim_fit,
    params              = list(
      intervention_type  = intervention_type,
      method             = "mc",
      n_sim              = n_sim,
      n_steps            = n_steps,
      burn_in            = burn_in,
      warmup             = warmup,
      sampling_steps     = sampling_steps,
      scope              = if (use_pop) "population" else "typical",
      intervention_sd    = intervention_sd,
      delta              = delta,
      dose_scale         = dose_scale,
      target_dp          = target_dp,
      calibrate_on       = calibrate_on,
      calibration        = calib,
      paired             = paired,
      trajectory         = keep_traj,
      diagnose           = do_diag,
      initialization     = init_kind,
      seed               = seed
    )
  )
  class(result) <- "kim_nira_result"

  if (verbose) {
    message("  Done. Top-ranked target (largest estimated benefit): [", result$ranked_targets[1], "]",
            "  (effect = ", result$summary$mean_diff[1], ")")
  }
  result
}


# -----------------------------------------------------------------------------
# S3 print method
# -----------------------------------------------------------------------------

#' @export
print.kim_nira_result <- function(x, ...) {
  p <- x$params
  cat("\n-- KIM NodeIdentifyR Algorithm (NIRA) Result --\n")
  cat("Method       : ", p$method, "\n")
  cat("Intervention : ", p$intervention_type, "\n")
  cat("Nodes        : ", paste(x$kim_fit$vars, collapse = ", "), "\n")
  cat("Effect       : ",
      if (!is.null(x$population_averaged) && isTRUE(x$population_averaged))
        "population-average (marginal)" else "typical-person", "\n")
  if (identical(p$method, "mc")) {
    cat("n_sim =", p$n_sim, " | n_steps =", p$n_steps,
        " | burn_in =", p$burn_in, "\n")
    cat("Paired design:", if (isTRUE(p$paired)) "yes" else "no", "\n")
  }
  if (length(unique(x$delta)) == 1L) cat("delta        =", round(x$delta[1], 4), "logits\n\n") else {
    cat("delta (per target, logits):\n"); print(round(x$delta, 4)); cat("\n") }
  cat("Probability change table (conditional, at the typical intercept):\n")
  print(x$delta_prob_table[, c("node", "p_original", "p_intervened", "delta_p")],
        row.names = FALSE)
  cat("\nTarget ranking (sorted by directional benefit):\n")
  cols <- intersect(c("rank", "node", "mean_diff", "direct_contribution",
                      "spillover_contribution", "spillover_mean_diff",
                      "spillover_rank", "delta_p"), names(x$summary))
  print(x$summary[, cols, drop = FALSE], row.names = FALSE)
  cat("\nTop-ranked target (largest estimated benefit; see kimStabilityBootstrap() for whether it is distinguishable from the others): [", x$ranked_targets[1], "]\n")
  cat("\nNote: mean_diff = change in mean activation rate (post-burn_in steps,\n")
  cat("      averaged over persons when population-average). delta_p = change in\n")
  cat("      the CONDITIONAL activation probability at the all-zero preceding\n")
  cat("      state; it is not a marginal or stationary quantity. Sampling\n")
  cat("      uncertainty is reported by kimStabilityBootstrap() (SE/CI/percenttop_1).\n")
  cat("\n")
  invisible(x)
}
