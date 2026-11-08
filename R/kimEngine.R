# =============================================================================
# kimEngine.R -- internal synchronous finite-state engine
#
# Convention: J[i, j] is the effect of node j at t on node i at t+1 (column ->
# row). States are enumerated as the rows of .kim_state_matrix(p).
# =============================================================================

.kim_state_cache <- new.env(parent = emptyenv())

.kim_state_matrix <- function(p) {
  if (length(p) != 1 || p < 1 || p != as.integer(p)) stop("Invalid state dimension.")
  key <- as.character(p)
  if (is.null(.kim_state_cache[[key]]))
    .kim_state_cache[[key]] <- vapply(seq_len(p),
      function(j) as.integer((0:(2^p - 1)) %/% 2^(j - 1) %% 2), integer(2^p))
  .kim_state_cache[[key]]
}

# Soft size limit for the dense engine. The kernel is a dense 2^p x 2^p matrix,
# i.e. 8 * 4^p bytes: 8 MB at p = 10, 128 MB at p = 12, 2 GB at p = 14. Beyond
# the limit an allocation usually drives the machine into swap before R reports
# an error, so the engine refuses up front. The limit is an option, not a
# constant: options(KimNIRA.max_exact_p = 14) raises it deliberately.
.kim_check_exact_p <- function(p) {
  lim <- getOption("KimNIRA.max_exact_p", 12L)
  if (p > lim)
    stop("The exact engine builds a dense 2^p x 2^p kernel (", format(8 * 4^p / 2^20,
         big.mark = ","), " MB at p = ", p, "). This exceeds the limit p <= ", lim,
         ". Use method = 'mc', or raise options(KimNIRA.max_exact_p = ", p,
         ") if the machine has the memory.", call. = FALSE)
  invisible(TRUE)
}

# Recycle a scalar dose to every target, or validate a per-target dose vector.
.kim_expand_delta <- function(delta, vars) {
  p <- length(vars)
  if (length(delta) == 1L) delta <- rep(delta, p)
  if (length(delta) != p || any(!is.finite(delta)))
    stop("delta must be one finite value or one finite value per node.", call. = FALSE)
  if (!is.null(names(delta)) && all(nzchar(names(delta)))) {
    if (!setequal(names(delta), vars)) stop("Named delta must name every node exactly once.")
    delta <- delta[vars]
  }
  stats::setNames(as.numeric(delta), vars)
}

.kim_transition <- function(J, b) {
  p <- nrow(J)
  if (!is.matrix(J) || ncol(J) != p || length(b) != p || any(!is.finite(c(J, b))))
    stop("Invalid parameters.")
  .kim_check_exact_p(p)
  S <- .kim_state_matrix(p); eta <- sweep(S %*% t(J), 2, b, "+")
  # log sigma(eta) and log(1 - sigma(eta)) without overflow.
  lp <- -pmax(-eta, 0) - log1p(exp(-abs(eta)))
  lq <- -pmax(eta, 0) - log1p(exp(-abs(eta)))
  P <- exp(lp %*% t(S) + lq %*% t(1 - S))
  err <- max(abs(rowSums(P) - 1))
  if (!is.finite(err) || err > 1e-10) stop("Transition row normalization failed.")
  list(P = P, states = S, probabilities = stats::plogis(eta), row_error = err)
}

# Stationary distribution by direct solution of the constrained linear system.
#
# Existence and uniqueness. With finite J and b every conditional probability
# lies strictly inside (0, 1), so every entry of P is strictly positive. A
# positive stochastic matrix is irreducible and aperiodic, and by Perron-Frobenius
# its stationary distribution exists, is unique, and is strictly positive. This
# holds for every finite parameter value; it is a property of the model, not of
# the data.
#
# Solution. pi P = pi is (P^T - I) pi^T = 0, whose null space is one-dimensional.
# The last row is replaced by the normalisation 1^T pi^T = 1. The all-ones row
# cannot lie in the row space of (P^T - I) because every vector there is
# orthogonal to pi while 1^T pi = 1, so the replaced system is nonsingular.
#
# Accuracy. The solve has no iteration and no convergence tolerance, so there is
# no Monte Carlo or truncation error. It is not error-free: the forward error is
# bounded by roughly kappa(A) * machine epsilon, and kappa(A) grows as the
# spectral gap 1 - |lambda_2| shrinks, i.e. in strongly metastable chains whose
# logits are very large. The residual max|pi P - pi| reported here is a backward
# error; it can be small while pi is inaccurate in an ill-conditioned system.
# `condition = TRUE` additionally reports the reciprocal condition number.
.kim_stationary <- function(P, tol = 1e-12, condition = FALSE) {
  ns <- nrow(P)
  A <- t(P) - diag(ns); A[ns, ] <- 1
  v <- tryCatch(as.numeric(solve(A, c(rep(0, ns - 1), 1))),
                error = function(e) NULL, warning = function(w) NULL)
  if (is.null(v) || length(v) != ns || !all(is.finite(v)))
    stop("Stationary solve failed: the constrained system is numerically singular. ",
         "With finite parameters P is strictly positive and pi is unique, so this ",
         "signals extreme logits (conditional probabilities rounding to 0 or 1) ",
         "rather than a reducible chain. Inspect very large |J| or |h|.", call. = FALSE)
  if (any(v < -tol))
    stop("Stationary solver returned a negative probability (minimum ",
         format(min(v), digits = 3), "); the system is too ill-conditioned for a ",
         "reliable answer.", call. = FALSE)
  v <- pmax(v, 0); v <- v / sum(v)
  res <- max(abs(as.numeric(v %*% P) - v))
  if (!is.finite(res) || res >= tol)
    stop("Stationary solve failed its residual check (max|pi P - pi| = ",
         format(res, digits = 3), ").", call. = FALSE)
  out <- list(pi = v, residual = res, iterations = 1L, solver = "linear")
  if (condition) out$rcond <- rcond(A)
  out
}

.kim_stationary_rate <- function(J, b, delta = 0, target = NULL) {
  if (!is.null(target)) b[target] <- b[target] + delta
  tr <- .kim_transition(J, b); st <- .kim_stationary(tr$P)
  as.numeric(st$pi %*% rowMeans(tr$states))
}

.kim_baselines <- function(fit) {
  shifts <- fit$person_shifts
  if (is.null(shifts) || !nrow(shifts)) return(matrix(fit$h, nrow = 1))
  sweep(shifts, 2, fit$h, "+")
}

# The transition matrix of one person. Shared J unless the fit carries
# person-specific slopes, in which case person k propagates through their own.
.kim_J_for <- function(fit, k) {
  pJ <- fit$person_J
  if (is.null(pJ)) return(fit$J)
  matrix(pJ[k, , , drop = TRUE], nrow = nrow(fit$J), ncol = ncol(fit$J),
         dimnames = dimnames(fit$J))
}

# One definition of "does this fit describe a population of person-specific
# models?". Both engines call this, so they cannot disagree about the estimand.
.kim_uses_population <- function(fit) {
  sh <- fit$person_shifts
  (!is.null(sh) && nrow(sh) >= 1L && any(rowSums(abs(sh)) > 0)) || !is.null(fit$person_J)
}

# Distinct (baseline, J) person patterns and their population weights. Persons
# sharing both are one pattern, so the dense solve is done once per pattern.
.kim_patterns <- function(fit) {
  bs <- .kim_baselines(fit)
  bkey <- apply(bs, 1, function(z) paste(format(z, digits = 17), collapse = "|"))
  jkey <- if (is.null(fit$person_J)) rep("shared", nrow(bs)) else
    apply(fit$person_J, 1, function(M) paste(format(M, digits = 17), collapse = "|"))
  keys <- paste(bkey, jkey, sep = "@@")
  first <- !duplicated(keys); rep_idx <- which(first)
  list(index = rep_idx, weights = tabulate(match(keys, keys[first])) / nrow(bs),
       baselines = bs)
}

# Exact stationary intervention effects. `delta` is a scalar (one dose for every
# target) or a per-target vector; `J_override(j, Jk)` optionally edits the
# transition matrix used under the intervention on target j (used for the
# knock-out decomposition, never by runKIMNIRA itself).
.kim_exact_effects <- function(fit, delta, J_override = NULL) {
  p <- nrow(fit$J); vars <- fit$vars
  .kim_check_exact_p(p)
  delta <- .kim_expand_delta(delta, vars)
  pat <- .kim_patterns(fit); bs <- pat$baselines
  means <- matrix(0, length(pat$index), p + 1); node <- matrix(0, p + 1, p)
  diags <- list(); di <- 0L
  for (kk in seq_along(pat$index)) {
    k <- pat$index[kk]; Jk <- .kim_J_for(fit, k); b <- bs[k, ]
    for (j in 0:p) {
      bb <- b; JJ <- Jk
      if (j) { bb[j] <- bb[j] + delta[j]; if (!is.null(J_override)) JJ <- J_override(j, Jk) }
      tr <- .kim_transition(JJ, bb); st <- .kim_stationary(tr$P)
      nm <- as.numeric(st$pi %*% tr$states)
      means[kk, j + 1] <- mean(nm); node[j + 1, ] <- node[j + 1, ] + pat$weights[kk] * nm
      di <- di + 1L; diags[[di]] <- data.frame(person_pattern = kk, target = j,
        residual = st$residual, iterations = st$iterations, solver = st$solver)
    }
  }
  total <- as.numeric(pat$weights %*% means)
  nc <- sweep(node[-1, , drop = FALSE], 2, node[1, ], "-")
  dimnames(nc) <- list(vars, vars)
  list(base = total[1], effect = stats::setNames(total[-1] - total[1], vars),
       node_baseline = stats::setNames(node[1, ], vars), node_changes = nc,
       delta = delta, diagnostics = do.call(rbind, diags),
       population_averaged = .kim_uses_population(fit))
}

.kim_bc_int <- function(x, theta0, alpha = .05) {
  prop <- mean(x < theta0) + .5 * mean(x == theta0)
  z0 <- stats::qnorm(pmin(pmax(prop, 1 / (2 * length(x))), 1 - 1 / (2 * length(x))))
  stats::quantile(x, stats::pnorm(2 * z0 + stats::qnorm(c(alpha / 2, 1 - alpha / 2))), names = FALSE)
}

# Finite-time response B_j(m) from the baseline stationary distribution, and
# optionally the linear-response derivative d mu / d b_j at delta = 0.
.kim_response <- function(fit, delta = -.5, horizon = 30L, derivative = FALSE) {
  bs <- .kim_baselines(fit); p <- nrow(fit$J)
  delta <- .kim_expand_delta(delta, fit$vars)
  changes <- array(0, c(horizon + 1, p, p), dimnames = list(0:horizon, fit$vars, fit$vars))
  grad <- numeric(p)
  for (k in seq_len(nrow(bs))) {
    Jk <- .kim_J_for(fit, k)
    tr <- .kim_transition(Jk, bs[k, ]); pi <- .kim_stationary(tr$P)$pi
    base <- as.numeric(pi %*% tr$states)
    if (derivative) {
      f <- rowMeans(tr$states); mu <- sum(pi * f)
      z <- solve(diag(length(pi)) - tr$P + matrix(pi, length(pi), length(pi), byrow = TRUE), f - mu)
      for (j in seq_len(p)) {
        dP <- tr$P * (matrix(tr$states[, j], nrow = nrow(tr$P), ncol = ncol(tr$P), byrow = TRUE) -
                        tr$probabilities[, j])
        grad[j] <- grad[j] + as.numeric(pi %*% dP %*% z) / nrow(bs)
      }
    }
    for (j in seq_len(p)) {
      b <- bs[k, ]; b[j] <- b[j] + delta[j]; Pj <- .kim_transition(Jk, b)$P; v <- pi
      for (t in seq_len(horizon)) {
        v <- as.numeric(v %*% Pj)
        changes[t + 1, j, ] <- changes[t + 1, j, ] + (as.numeric(v %*% tr$states) - base) / nrow(bs)
      }
    }
  }
  list(node_changes = changes, mean_changes = apply(changes, c(1, 2), mean),
       derivative = stats::setNames(grad, fit$vars), initialization = "baseline stationary")
}

# Second-largest eigenvalue modulus of a stochastic matrix and the implied
# relaxation time -1 / log|lambda_2| (in transitions). Exact, unlike the
# weighted-Jacobian heuristic, but costs a dense eigen decomposition.
.kim_spectral <- function(P) {
  ev <- eigen(P, only.values = TRUE)$values
  mods <- sort(Mod(ev), decreasing = TRUE)
  l2 <- if (length(mods) > 1L) mods[2] else 0
  list(lambda2 = l2, gap = 1 - l2,
       relaxation = if (l2 <= 0) 0 else if (l2 >= 1) Inf else -1 / log(l2))
}
