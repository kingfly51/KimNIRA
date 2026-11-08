#' Check the conditional-independence (no contemporaneous dependence) assumption
#'
#' The KIM transition kernel assumes that, given s(t), the nodes at t+1 are
#' independent. This function adds the other nodes' contemporaneous states
#' s_j(t+1) to each node's lagged logistic regression and tests whether they
#' carry information beyond s(t). Contemporaneous states are centred within
#' person and the person means of both the lagged and the contemporaneous
#' predictors are included, so a positive result is not produced merely by
#' persons who are high on everything (between-person association).
#'
#' This is a diagnostic of the working assumption, not an estimator. Unpenalised
#' `glm` fits are used, so it is intended for the moderate p at which the exact
#' engine is used. Nonzero within-person contemporaneous coefficients mean that
#' part of the dependence between symptoms occurs within one observation
#' interval; the synchronous KIM then treats that dependence as absent, and the
#' spillover part of NIRA effects is typically understated when the
#' contemporaneous associations are positive.
#'
#' @param kim_fit A `"dynaNIRA_kim"` fit (its time rules are reused).
#' @param data The data used to fit it.
#' @param p_adjust Multiplicity adjustment across nodes (default `"holm"`).
#' @return A list with `tests` (one row per node: likelihood-ratio statistic,
#'   df, p-value and adjusted p-value for the block of within-person
#'   contemporaneous terms) and `coefficients` (p x p matrix; entry [i, j] is
#'   the within-person log-odds association of s_j(t+1) with s_i(t+1) given
#'   s(t); the diagonal is NA).
#' @export
kimContemporaneousCheck <- function(kim_fit, data, p_adjust = "holm") {
  if (!inherits(kim_fit, "dynaNIRA_kim") || is.null(kim_fit$settings))
    stop("kim_fit must be a fitKIM object with saved settings.")
  spec <- kim_fit$settings
  args <- spec[intersect(names(spec), names(formals(.kim_time_pairs)))]
  args$data <- data
  pr <- do.call(.kim_time_pairs, args)
  keep <- pr$person %in% kim_fit$person_ids
  X <- pr$X[[1L]][keep, , drop = FALSE]; Y <- pr$Y[keep, , drop = FALSE]
  person <- pr$person[keep]; vars <- spec$vars; p <- length(vars)
  storage.mode(X) <- "double"; storage.mode(Y) <- "double"
  pmean <- function(M) { m <- rowsum(M, person) / as.numeric(table(person)[rownames(rowsum(M, person))])
    m[person, , drop = FALSE] }
  Xm <- pmean(X); Ym <- pmean(Y); Yw <- Y - Ym
  colnames(X) <- paste0("lag_", vars); colnames(Xm) <- paste0("lagmean_", vars)
  colnames(Yw) <- paste0("now_", vars); colnames(Ym) <- paste0("nowmean_", vars)
  coefs <- matrix(NA_real_, p, p, dimnames = list(vars, vars))
  tests <- do.call(rbind, lapply(seq_len(p), function(i) {
    y <- Y[, i]
    if (stats::var(y) < 1e-10)
      return(data.frame(node = vars[i], lr = NA_real_, df = p - 1L, p_value = NA_real_))
    base <- data.frame(y = y, X, Xm, Ym[, -i, drop = FALSE])
    full <- data.frame(base, Yw[, -i, drop = FALSE])
    f0 <- suppressWarnings(stats::glm(y ~ ., family = stats::binomial(), data = base))
    f1 <- suppressWarnings(stats::glm(y ~ ., family = stats::binomial(), data = full))
    cf <- stats::coef(f1)[paste0("now_", vars[-i])]
    coefs[i, -i] <<- cf
    lr <- max(0, f0$deviance - f1$deviance)
    data.frame(node = vars[i], lr = lr, df = p - 1L,
               p_value = stats::pchisq(lr, p - 1L, lower.tail = FALSE))
  }))
  tests$p_adjusted <- stats::p.adjust(tests$p_value, method = p_adjust)
  list(tests = tests, coefficients = coefs,
       note = paste("Within-person contemporaneous associations given s(t); person means of",
                    "lagged and contemporaneous states are controlled. A significant block",
                    "indicates dependence within one observation interval, which the",
                    "synchronous KIM does not represent."))
}
