# =============================================================================
# kimLandscape.R -- multistability, stability landscape and mixing of the KIM
#
# Why this file exists. With finite parameters every entry of the KIM transition
# kernel is strictly positive, so the stationary distribution pi is unique
# (Perron-Frobenius). A fitted KIM can therefore never be bistable in the strict
# sense of two absorbing equilibria. What it can show is METASTABILITY: pi puts
# its mass in two (or more) separated regions of the state space, the chain
# lingers in one region for a long time before crossing to the other, and the
# second eigenvalue of P approaches one. All three signatures are computed
# exactly here from the same finite-state kernel the exact NIRA engine uses:
#
#   * the stationary distribution of the total activation K = sum_i s_i and its
#     potential U(K) = -log P(K) (a stability landscape in the sense of
#     Cui et al., 2026, but obtained from the solved stationary law instead of
#     simulated frequencies);
#   * the modes of P(K) and the barrier between adjacent modes;
#   * the relaxation time -1/log|lambda_2(P)|.
#
# A useful reference point: with J = 0 the nodes are independent and K follows
# a Poisson-binomial law, which is log-concave and therefore unimodal. Any
# bimodality of P(K) within a person is thus produced by the interaction
# structure. Bimodality of the POPULATION mixture, by contrast, can also arise
# from persons with different baselines, so it is reported separately.
#
# The mean-field map m -> sigma(b + J m) is a deterministic skeleton of the
# chain. Its stable fixed points indicate where the stochastic chain is drawn;
# several stable fixed points are a mean-field signature of multistability that
# remains computable when p is too large for the exact engine. It is an
# approximation (it ignores fluctuations), not a property of the fitted chain.
# =============================================================================


# ---- mean-field skeleton ----------------------------------------------------

.kim_mf_map <- function(J, b, m) stats::plogis(b + as.numeric(J %*% m))

# One fixed point by damped iteration (same fixed points as the undamped map).
.kim_mf_fixed_point <- function(J, b, start = stats::plogis(b), maxit = 5000L,
                                tol = 1e-10) {
  m <- as.numeric(start)
  for (it in seq_len(maxit)) {
    mn <- 0.5 * m + 0.5 * .kim_mf_map(J, b, m)
    if (max(abs(mn - m)) < tol) return(mn)
    m <- mn
  }
  m
}

# Stable fixed points reached from dispersed starts, refined by Newton steps.
.kim_mf_fixed_points <- function(J, b, n_starts = 50L, seed = 1L) {
  p <- length(b)
  old <- if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) get(".Random.seed", .GlobalEnv) else NULL
  on.exit(if (is.null(old)) { if (exists(".Random.seed", .GlobalEnv, inherits = FALSE))
    rm(".Random.seed", envir = .GlobalEnv) } else assign(".Random.seed", old, .GlobalEnv), add = TRUE)
  set.seed(seed)
  starts <- c(list(rep(0, p), rep(1, p), stats::plogis(b)),
              lapply(seq_len(n_starts), function(i) stats::runif(p)))
  fps <- list()
  for (st in starts) {
    m <- .kim_mf_fixed_point(J, b, st)
    for (it in 1:20) {                          # Newton refinement on m - F(m)
      Fm <- .kim_mf_map(J, b, m); g <- m - Fm
      if (max(abs(g)) < 1e-13) break
      Jac <- diag(p) - diag(Fm * (1 - Fm), p) %*% J
      step <- tryCatch(solve(Jac, g), error = function(e) NULL)
      if (is.null(step)) break
      m <- pmin(pmax(m - step, 0), 1)
    }
    if (max(abs(m - .kim_mf_map(J, b, m))) > 1e-8) next
    if (!any(vapply(fps, function(f) max(abs(f - m)) < 1e-6, logical(1))))
      fps[[length(fps) + 1L]] <- m
  }
  if (!length(fps)) return(data.frame())
  out <- do.call(rbind, lapply(fps, function(m) {
    d <- m * (1 - m)
    rho <- max(abs(eigen(diag(d, p) %*% J, only.values = TRUE)$values))
    data.frame(mean_activation = mean(m), rho_local = rho, stable = rho < 1,
               t(stats::setNames(m, names(b))), check.names = FALSE)
  }))
  out <- out[out$stable, , drop = FALSE]
  out[order(out$mean_activation), , drop = FALSE]
}


# ---- modes of a discrete distribution ---------------------------------------

# Local maxima of q (indexed 0..p) after merging pairs separated by a barrier
# smaller than `min_barrier` (natural-log units of probability). The barrier
# between adjacent modes a < b is log(min(q_a, q_b) / min(q[a..b])).
.kim_modes <- function(q, min_barrier = 0.1) {
  n <- length(q)
  is_max <- vapply(seq_len(n), function(i)
    (i == 1L || q[i] > q[i - 1L]) && (i == n || q[i] >= q[i + 1L]), logical(1))
  modes <- which(is_max)
  repeat {
    if (length(modes) < 2L) break
    bar <- vapply(seq_len(length(modes) - 1L), function(a) {
      lo <- min(q[modes[a]:modes[a + 1L]])
      log(min(q[modes[a]], q[modes[a + 1L]]) / lo)
    }, numeric(1))
    if (min(bar) >= min_barrier) break
    a <- which.min(bar)
    drop <- if (q[modes[a]] < q[modes[a + 1L]]) a else a + 1L
    modes <- modes[-drop]
  }
  troughs <- barriers <- numeric(0)
  if (length(modes) >= 2L) for (a in seq_len(length(modes) - 1L)) {
    seg <- modes[a]:modes[a + 1L]
    tr <- seg[which.min(q[seg])]
    troughs <- c(troughs, tr - 1L)
    barriers <- c(barriers, log(min(q[modes[a]], q[modes[a + 1L]]) / q[tr]))
  }
  list(modes = modes - 1L, troughs = troughs, barriers = barriers,
       n_modes = length(modes))
}

# Stationary summaries of one (J, b): the count distribution and, optionally,
# the spectral gap.
.kim_landscape_one <- function(J, b, spectral = FALSE) {
  tr <- .kim_transition(J, b); st <- .kim_stationary(tr$P)
  K <- rowSums(tr$states); p <- ncol(tr$states)
  q <- as.numeric(tapply(st$pi, factor(K, levels = 0:p), sum))
  out <- list(pi = st$pi, count = q, mean = sum(q * (0:p)) / p,
              node = as.numeric(st$pi %*% tr$states))
  if (spectral) out$spectral <- .kim_spectral(tr$P)
  out
}

# Local extrema of a vector, including one-sided boundary extrema. These match
# Isinglandr's local_min_index()/local_max_index(), i.e. the selection rule used
# by Cui et al. (2026).
.kim_local_min_index <- function(v, include_margin = TRUE) {
  bounded_v <- if (include_margin) c(Inf, v, Inf) else c(-Inf, v, -Inf)
  diff1 <- diff(bounded_v)[-1]
  diff2 <- rev(diff(rev(bounded_v))[-1])
  which(diff1 >= 0 & diff2 >= 0)
}

.kim_local_max_index <- function(v, include_margin = TRUE) {
  bounded_v <- if (include_margin) c(-Inf, v, -Inf) else c(Inf, v, Inf)
  diff1 <- diff(bounded_v)[-1]
  diff2 <- rev(diff(rev(bounded_v))[-1])
  which(diff1 <= 0 & diff2 <= 0)
}

# Cui et al. (2026) phase stability from q = P(K), K = 0..p. The first phase is
# K < threshold and the second K >= threshold. Within a phase, stability is
# U(max) - U(min) with U = -log q; the minimum farther from the other phase and
# the facing maximum are used. Indices are returned both as 1-based rows of the
# full K = 0..p table and as 0-based K positions.
.kim_phase_stability_q <- function(q, threshold) {
  p <- length(q) - 1L
  U <- -log(q)
  one <- seq_len(threshold)
  two <- (threshold + 1L):(p + 1L)

  phase <- function(u, first) {
    mins <- .kim_local_min_index(u)
    if (!length(mins))
      return(list(stability = NA_real_, imin = NA_integer_, imax = NA_integer_,
                  status = "no_local_minimum"))
    imin <- if (first) mins[1L] else mins[length(mins)]
    maxs <- .kim_local_max_index(u)
    maxs <- if (first) maxs[maxs > imin] else maxs[maxs < imin]
    if (!is.finite(u[imin]))
      return(list(stability = NA_real_, imin = imin, imax = NA_integer_,
                  status = "nonfinite_potential"))
    if (!length(maxs)) {
      if (any(!is.finite(u)))
        return(list(stability = NA_real_, imin = imin, imax = NA_integer_,
                    status = "nonfinite_potential"))
      return(list(stability = 0, imin = imin, imax = imin,
                  status = "no_facing_maximum"))
    }
    imax <- if (first) maxs[1L] else maxs[length(maxs)]
    if (!is.finite(u[imax]))
      return(list(stability = NA_real_, imin = imin, imax = imax,
                  status = "nonfinite_potential"))
    list(stability = u[imax] - u[imin], imin = imin, imax = imax, status = "ok")
  }

  a <- phase(U[one], TRUE)
  b <- phase(U[two], FALSE)
  full <- function(z, offset) if (is.na(z)) NA_integer_ else z + offset
  kk <- function(z) if (is.na(z)) NA_integer_ else z - 1L
  status <- if (identical(a$status, "ok") && identical(b$status, "ok")) "ok" else
    paste0("first:", a$status, "; second:", b$status)
  list(stability1 = a$stability,
       stability2 = b$stability,
       stability_diff = a$stability - b$stability,
       effective_minindex1 = full(a$imin, 0L),
       effective_maxindex1 = full(a$imax, 0L),
       effective_minindex2 = full(b$imin, threshold),
       effective_maxindex2 = full(b$imax, threshold),
       K_min1 = kk(full(a$imin, 0L)),
       K_max1 = kk(full(a$imax, 0L)),
       K_min2 = kk(full(b$imin, threshold)),
       K_max2 = kk(full(b$imax, threshold)),
       status = status)
}

# Hamming-local maxima of a state distribution: s is a mode when no state that
# differs from it in one node has higher probability.
.kim_state_modes <- function(pi, p, vars, min_mass = 0.01) {
  S <- .kim_state_matrix(p); idx <- seq_along(pi) - 1L
  is_mode <- vapply(seq_along(pi), function(r) {
    nb <- bitwXor(idx[r], 2L^(0:(p - 1L))) + 1L
    all(pi[r] >= pi[nb])
  }, logical(1))
  r <- which(is_mode & pi >= min_mass)
  if (!length(r)) return(data.frame())
  data.frame(state = vapply(r, function(i) paste(vars[S[i, ] == 1L], collapse = "+"), ""),
             n_active = rowSums(S[r, , drop = FALSE]), probability = pi[r],
             stringsAsFactors = FALSE)[order(-pi[r]), ]
}


#' Metastability and stability landscape of a fitted KIM (exact)
#'
#' Computes, from the exact stationary distribution, whether the fitted
#' dynamics concentrate in separated low- and high-activation regions, how
#' strongly, and how a sustained intervention on each node reshapes that
#' landscape.
#'
#' A fitted KIM with finite parameters has a unique, strictly positive
#' stationary distribution, so it cannot be bistable in the strict sense of two
#' absorbing equilibria. Bistability in the clinically relevant sense appears as
#' metastability: a bimodal stationary distribution of the total activation
#' K = sum_i s_i, a barrier between the modes, and a second eigenvalue of P
#' close to one. All are exact here.
#'
#' @param x A `"dynaNIRA_kim"` fit or a `"kim_nira_result"`. With a result, its
#'   fit, dose and intervention direction are reused.
#' @param delta Optional dose (scalar or one per node). When given, each target in
#'   `targets` is intervened on and the landscape recomputed.
#' @param targets Node names to intervene on (default: all nodes when a dose is
#'   available, none otherwise).
#' @param min_barrier Smallest barrier, in natural-log units, for two local
#'   maxima of P(K) to count as separate modes (default 0.1, i.e. the trough is
#'   at least about one tenth below the smaller peak).
#' @param spectral `"typical"` (default) computes |lambda_2| and the relaxation
#'   time for the typical person, `"all"` for every person pattern, `"none"` skips.
#' @param mean_field Also locate the stable fixed points of the mean-field map.
#' @param n_starts Random starts for the mean-field search.
#' @param verbose Print a summary.
#' @return A `"kim_bistability"` list: `count_distribution` (long table of P(K) and
#'   U(K) = -log P(K) for the baseline and each intervention), `modes`,
#'   `persons` (per-person modes and barrier heights), `population_vs_within`
#'   (whether population bimodality is reproduced within persons),
#'   `state_modes`, `interventions`, `mean_field`, `spectral`, and settings.
#' @references Cui, J., Lunansky, G., Lichtwarck-Aschoff, A., Mendoza, N., &
#'   Hasselman, F. (2026). Quantifying the stability landscapes of psychological
#'   networks. Behavior Research Methods, 58, 68.
#' @export
kimBistability <- function(x, delta = NULL, targets = NULL,
                           min_barrier = 0.1,
                           spectral = c("typical", "all", "none"),
                           mean_field = TRUE, n_starts = 50L, verbose = TRUE) {
  spectral <- match.arg(spectral)
  if (inherits(x, "kim_nira_result")) {
    fit <- x$kim_fit
    if (is.null(delta)) delta <- x$delta
  } else if (inherits(x, "dynaNIRA_kim")) fit <- x else
    stop("x must be a fitKIM() fit or a runKIMNIRA() result.")
  vars <- fit$vars; p <- length(vars)
  .kim_check_exact_p(p)
  if (!is.null(delta)) delta <- .kim_expand_delta(delta, vars)
  if (is.null(targets)) targets <- if (is.null(delta)) character(0) else vars
  if (length(targets) && is.null(delta)) stop("Supply a dose (delta) to intervene.")
  if (!all(targets %in% vars)) stop("Unknown target node(s).")

  pat <- .kim_patterns(fit); bs <- pat$baselines; w <- pat$weights
  n_pat <- length(pat$index)

  # ---- baseline: distributions ------------------------------------------------
  base <- lapply(pat$index, function(k)
    .kim_landscape_one(.kim_J_for(fit, k), bs[k, ]))
  q_pop <- Reduce(`+`, Map(function(z, wk) wk * z$count, base, w))
  pi_pop <- Reduce(`+`, Map(function(z, wk) wk * z$pi, base, w))
  md <- .kim_modes(q_pop, min_barrier)

  if (spectral == "all") {
    base <- lapply(seq_len(n_pat), function(kk) {
      k <- pat$index[kk]
      .kim_landscape_one(.kim_J_for(fit, k), bs[k, ], spectral = TRUE)
    })
  }
  person_modes <- lapply(base, function(z) .kim_modes(z$count, min_barrier))
  persons <- data.frame(pattern = seq_len(n_pat), weight = w,
    mean_activation = vapply(base, `[[`, numeric(1), "mean"),
    n_modes = vapply(person_modes, `[[`, integer(1), "n_modes"),
    modes = vapply(person_modes, function(m) paste(m$modes, collapse = ","), ""),
    max_barrier = vapply(person_modes, function(m) if (length(m$barriers)) max(m$barriers) else 0, 0))
  if (spectral == "all") {
    persons$lambda2 <- vapply(base, function(z) z$spectral$lambda2, 0)
    persons$relaxation <- vapply(base, function(z) z$spectral$relaxation, 0)
  }
  within_share <- sum(w[persons$n_modes >= 2L])
  pop_bimodal <- md$n_modes >= 2L
  origin <- if (!pop_bimodal) "unimodal" else if (within_share >= 0.5)
    "within-person (dynamics)" else if (within_share > 0)
    "mixed: partly within-person, partly between-person heterogeneity" else
    "between-person heterogeneity (no person is bimodal)"

  # ---- interventions -----------------------------------------------------------
  count_long <- data.frame(condition = "baseline", K = 0:p, probability = q_pop,
                           stringsAsFactors = FALSE)
  iv <- NULL
  base_summary <- function(lst) {
    list(mean_activation = sum(w * vapply(lst, `[[`, 0, "mean")))
  }
  b0 <- base_summary(base)
  if (length(targets)) {
    rows <- list()
    for (tg in targets) {
      j <- match(tg, vars)
      lst <- lapply(pat$index, function(k) {
        bb <- bs[k, ]; bb[j] <- bb[j] + delta[j]
        .kim_landscape_one(.kim_J_for(fit, k), bb)
      })
      q <- Reduce(`+`, Map(function(z, wk) wk * z$count, lst, w))
      mj <- .kim_modes(q, min_barrier); sj <- base_summary(lst)
      count_long <- rbind(count_long, data.frame(condition = tg, K = 0:p,
                                                 probability = q, stringsAsFactors = FALSE))
      row <- data.frame(target = tg, delta = delta[j],
        mean_activation = sj$mean_activation,
        effect = sj$mean_activation - b0$mean_activation,
        n_modes = mj$n_modes,
        max_barrier = if (length(mj$barriers)) max(mj$barriers) else 0,
        share_persons_bimodal = sum(w[vapply(lst, function(z)
          .kim_modes(z$count, min_barrier)$n_modes, 0L) >= 2L]),
        stringsAsFactors = FALSE)
      rows[[tg]] <- row
    }
    iv <- do.call(rbind, rows); rownames(iv) <- NULL
    iv <- iv[order(-abs(iv$effect)), , drop = FALSE]
  }
  count_long$U <- -log(count_long$probability)

  # ---- spectral (typical person) and mean field --------------------------------
  spec <- NULL
  if (spectral != "none") {
    sp <- .kim_spectral(.kim_transition(fit$J, fit$h)$P)
    spec <- data.frame(who = "typical person", lambda2 = sp$lambda2,
                       relaxation = sp$relaxation)
    if (spectral == "all")
      spec <- rbind(spec, data.frame(who = c("persons: median", "persons: max"),
        lambda2 = c(stats::median(persons$lambda2), max(persons$lambda2)),
        relaxation = c(stats::median(persons$relaxation), max(persons$relaxation))))
  }
  mf <- if (mean_field) .kim_mf_fixed_points(fit$J, fit$h, n_starts) else NULL

  out <- list(
    count_distribution = count_long,
    modes = data.frame(K = md$modes, probability = q_pop[md$modes + 1L]),
    troughs = md$troughs, barriers = md$barriers,
    population_bimodal = pop_bimodal,
    persons = persons,
    population_vs_within = list(share_persons_bimodal = within_share, origin = origin),
    state_modes = .kim_state_modes(pi_pop, p, vars),
    interventions = iv, mean_field = mf, spectral = spec,
    settings = list(delta = delta, targets = targets,
                    min_barrier = min_barrier,
                    population_averaged = .kim_uses_population(fit)),
    kim_fit = fit)
  class(out) <- "kim_bistability"
  if (verbose) print(out)
  invisible(out)
}

#' @export
print.kim_bistability <- function(x, ...) {
  p <- length(x$kim_fit$vars)
  cat("\n-- KIM stationary landscape (exact) --\n")
  cat("P(K), K = number of active nodes (0..", p, "):\n", sep = "")
  b <- x$count_distribution[x$count_distribution$condition == "baseline", ]
  print(round(stats::setNames(b$probability, b$K), 4))
  cat("Modes of P(K): ", paste(x$modes$K, collapse = ", "),
      if (length(x$barriers)) paste0("  | barrier(s) = ", paste(round(x$barriers, 3), collapse = ", ")), "\n", sep = "")
  cat("Population distribution: ", if (x$population_bimodal) "multimodal" else "unimodal",
      "; persons that are themselves multimodal: ",
      round(100 * x$population_vs_within$share_persons_bimodal, 1), "%  -> ",
      x$population_vs_within$origin, "\n", sep = "")
  if (!is.null(x$spectral)) { cat("Spectral:\n"); print(x$spectral, row.names = FALSE) }
  if (!is.null(x$mean_field) && nrow(x$mean_field))
    cat("Stable mean-field fixed points (typical person): ", nrow(x$mean_field),
        " at mean activation ", paste(round(x$mean_field$mean_activation, 3), collapse = ", "),
        "\n", sep = "")
  if (!is.null(x$interventions)) {
    cat("\nIntervention effects on the landscape:\n")
    print(x$interventions, row.names = FALSE, digits = 4)
  }
  invisible(x)
}


#' Phase stability of a KIM stability landscape (Cui et al., 2026)
#'
#' Splits each condition's potential landscape U(K) = -log P(K) at a threshold
#' number of active nodes and reports the two phase stabilities and their
#' difference, following Cui et al. (2026) and the selection rule of
#' \code{Isinglandr::calculate_stability()}.
#'
#' @param land A `"kim_bistability"` object from [kimBistability()]. The
#'   baseline and every intervention condition stored in its
#'   `count_distribution` are available.
#' @param conditions Conditions to score. Default: all conditions present
#'   (`"baseline"` plus the intervention targets).
#' @param threshold Integer K in `1..p`: the first phase is `K < threshold`, the
#'   second `K >= threshold` (Cui et al.'s `split_value`). Default:
#'   `ceiling(p / 2)`. The value is used only for the phase split; it does not
#'   change the landscape or the intervention effects stored in `land`.
#' @return A data frame with one row per condition: `condition`, `threshold`,
#'   `stability1` (low-activation phase), `stability2` (high-activation phase),
#'   `stability_diff` (`stability1 - stability2`, positive means the
#'   low-activation phase is more stable), the selected `K_min*`/`K_max*`
#'   positions, and `status`.
#' @details Within each phase stability is `U(max) - U(min)`. The first phase
#'   uses its leftmost local minimum and the first local maximum to its right;
#'   the second phase uses its rightmost local minimum and the last local
#'   maximum to its left. Boundary extrema count. If no facing maximum exists
#'   the stability is 0; a missing local minimum or a non-finite potential gives
#'   `NA`, with the reason in `status`.
#' @references Cui, J., Lunansky, G., Lichtwarck-Aschoff, A., Mendoza, N., &
#'   Hasselman, F. (2026). Quantifying the stability landscapes of psychological
#'   networks. Behavior Research Methods, 58, 68.
#' @seealso [kimBistability()], [kimPlotLandscape()]
#' @export
kimPhaseStability <- function(land, conditions = NULL, threshold = NULL) {
  if (!inherits(land, "kim_bistability"))
    stop("`land` must come from kimBistability().", call. = FALSE)
  d <- land$count_distribution
  available <- unique(d$condition)
  if (is.null(conditions)) conditions <- available
  if (!is.character(conditions) || !length(conditions) || anyNA(conditions))
    stop("`conditions` must be a non-empty character vector.", call. = FALSE)
  miss <- setdiff(conditions, available)
  if (length(miss))
    stop("Unknown condition(s): ", paste(miss, collapse = ", "),
         ". Available: ", paste(available, collapse = ", "), ".", call. = FALSE)
  p <- length(land$kim_fit$vars)
  if (is.null(threshold)) threshold <- as.integer(ceiling(p / 2))
  if (length(threshold) != 1L || is.na(threshold) || !is.finite(threshold) ||
      threshold != round(threshold) || threshold < 1 || threshold > p)
    stop("threshold must be an integer in 1..p.", call. = FALSE)
  threshold <- as.integer(threshold)
  rows <- lapply(conditions, function(cond) {
    dd <- d[d$condition == cond, , drop = FALSE]
    dd <- dd[order(dd$K), , drop = FALSE]
    if (nrow(dd) != p + 1L || !identical(as.integer(dd$K), 0:p))
      stop("count_distribution for condition '", cond, "' is incomplete.",
           call. = FALSE)
    r <- .kim_phase_stability_q(dd$probability, threshold)
    data.frame(condition = cond, threshold = threshold,
               stability1 = r$stability1, stability2 = r$stability2,
               stability_diff = r$stability_diff,
               K_min1 = r$K_min1, K_max1 = r$K_max1,
               K_min2 = r$K_min2, K_max2 = r$K_max2,
               status = r$status, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}


#' Plot the stationary stability landscape U(K) = -log P(K)
#'
#' @param x A `"kim_bistability"` object.
#' @param conditions Conditions to draw: `"baseline"` plus target names. Default:
#'   baseline and the three targets with the largest effects.
#' @param per_state Divide P(K) by choose(p, K) before taking -log, i.e. show the
#'   landscape per microstate. This removes the combinatorial term that makes
#'   intermediate K common simply because many states have that count.
#' @param type `"potential"` (default) or `"probability"`.
#' @return A ggplot object.
#' @export
kimPlotLandscape <- function(x, conditions = NULL, per_state = FALSE,
                             type = c("potential", "probability")) {
  if (!inherits(x, "kim_bistability")) stop("x must come from kimBistability().")
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required.")
  type <- match.arg(type)
  d <- x$count_distribution; p <- length(x$kim_fit$vars)
  if (is.null(conditions)) {
    top <- if (is.null(x$interventions)) character(0) else
      utils::head(x$interventions$target, 3L)
    conditions <- c("baseline", top)
  }
  d <- d[d$condition %in% conditions, , drop = FALSE]
  d$condition <- factor(d$condition, levels = conditions)
  pr <- if (per_state) d$probability / choose(p, d$K) else d$probability
  d$y <- if (type == "potential") -log(pr) else pr
  ylab <- if (type == "potential") {
    if (per_state) "U(K) = -log[P(K) / C(p, K)]" else "U(K) = -log P(K)"
  } else if (per_state) "P(K) / C(p, K)" else "P(K)"
  cols <- c("grey35", grDevices::hcl.colors(max(1L, length(conditions) - 1L), "Dark 3"))
  gg <- ggplot2::ggplot(d, ggplot2::aes(x = .data$K, y = .data$y, colour = .data$condition,
                                        linetype = .data$condition)) +
    ggplot2::geom_line(linewidth = 0.9) + ggplot2::geom_point(size = 2) +
    ggplot2::scale_colour_manual(values = stats::setNames(cols[seq_along(conditions)], conditions)) +
    ggplot2::scale_linetype_manual(values = stats::setNames(c("dashed",
      rep("solid", length(conditions) - 1L)), conditions)) +
    ggplot2::scale_x_continuous(breaks = 0:p) +
    ggplot2::labs(x = "K = number of active nodes", y = ylab, colour = NULL, linetype = NULL) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "bottom", panel.grid.minor = ggplot2::element_blank())
  gg
}


#' Exact mixing (relaxation) time of a fitted KIM
#'
#' Reports the second-largest eigenvalue modulus |lambda_2| of each person's
#' transition kernel and the relaxation time -1/log|lambda_2|, in transitions
#' (one transition = one observation interval of the data the model was fitted
#' on). This replaces the linearised weighted-Jacobian heuristic stored in the
#' fit, which is reported alongside for comparison. The relaxation time is the
#' e-folding time of the slowest-decaying deviation from stationarity; a
#' Monte Carlo warm-up of several relaxation times removes the start-up bias.
#'
#' @param fit A `"dynaNIRA_kim"` fit.
#' @param persons `"typical"` (the typical-person h and J only) or `"all"`
#'   (every distinct person pattern; costs one dense eigen-decomposition each).
#' @return A data frame (and a `summary` attribute).
#' @export
kimMixing <- function(fit, persons = c("typical", "all")) {
  persons <- match.arg(persons)
  if (!inherits(fit, "dynaNIRA_kim")) stop("fit must come from fitKIM().")
  .kim_check_exact_p(length(fit$vars))
  rows <- list()
  sp <- .kim_spectral(.kim_transition(fit$J, fit$h)$P)
  rows[[1]] <- data.frame(who = "typical", weight = NA_real_, lambda2 = sp$lambda2,
                          relaxation = sp$relaxation)
  if (persons == "all") {
    pat <- .kim_patterns(fit)
    for (kk in seq_along(pat$index)) {
      k <- pat$index[kk]
      sp <- .kim_spectral(.kim_transition(.kim_J_for(fit, k), pat$baselines[k, ])$P)
      rows[[kk + 1L]] <- data.frame(who = paste0("pattern ", kk), weight = pat$weights[kk],
                                    lambda2 = sp$lambda2, relaxation = sp$relaxation)
    }
  }
  out <- do.call(rbind, rows)
  attr(out, "summary") <- list(
    exact_max_relaxation = max(out$relaxation),
    heuristic_rho_star = fit$spectral_radius_star,
    heuristic_rho_star_max = fit$spectral_radius_star_max,
    heuristic_relaxation = fit$mixing_time,
    unit = "transitions (one observation interval)")
  out
}
