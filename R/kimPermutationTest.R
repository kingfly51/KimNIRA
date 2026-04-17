# -----------------------------------------------------------------------------
# kimPermutationTest() -- Single-node paired permutation test
# -----------------------------------------------------------------------------

#' Paired Permutation Test for KIM NodeIdentifyR Algorithm (NIRA) Intervention Effects
#'
#' @description
#' Tests whether a single node's intervention significantly shifts
#' the network's projected sum score, using a paired sign-flip
#' permutation design.
#'
#' Unlike the independent two-sample permutation in the original
#' \code{NIRApost::permutationTest}, this function exploits the paired
#' structure of KIM: baseline and intervention scores share the same
#' initial states (and, for \code{model = "kim"}, the same random number
#' stream), so the paired difference
#' \eqn{d_k = \text{intervention}_k - \text{baseline}_k} is a
#' within-subject quantity.
#'
#' Null hypothesis: the population mean of paired differences is zero.
#' Permutation: randomly flip the sign of each \eqn{d_k}.
#' p-value: \eqn{(k + 1) / (B + 1)} (plus-one, two-sided).
#'
#' This function is called internally by
#' \code{\link{kimPermutationNIRAtest}} for every node, but can also be
#' used stand-alone.
#'
#' @param baseline Numeric vector of baseline projected sum scores.
#' @param intervention Numeric vector of intervention projected sum scores.
#'   Same length as \code{baseline} (paired).
#' @param nPerm Positive integer, number of permutations. Default \code{5000}.
#'
#' @return A named \code{list} with 17 elements (matching NIRApost interface):
#' \describe{
#'   \item{\code{mean_original}}{Mean of baseline scores.}
#'   \item{\code{mean_other}}{Mean of intervention scores.}
#'   \item{\code{sd_original}, \code{sd_other}}{Standard deviations.}
#'   \item{\code{se_original}, \code{se_other}}{Standard errors.}
#'   \item{\code{ci_original}, \code{ci_other}}{95\% CI half-widths.}
#'   \item{\code{ciLower_original}, \code{ciUpper_original}}{95\% CI bounds.}
#'   \item{\code{ciLower_other}, \code{ciUpper_other}}{95\% CI bounds.}
#'   \item{\code{pooled_std_dev}}{SD of paired differences.}
#'   \item{\code{cohens_d}}{Paired Cohen's d = mean(d) / sd(d).}
#'   \item{\code{p_value}}{Two-sided permutation p-value.}
#'   \item{\code{n_permutations}}{Number of permutations.}
#'   \item{\code{permutation_distribution}}{Numeric vector of permuted
#'     |mean(d)| values.}
#' }
#'
#' @seealso \code{\link{kimPermutationNIRAtest}}, \code{\link{kimPlotNIRA}}
#' @importFrom stats sd
#' @export
kimPermutationTest <- function(baseline, intervention, nPerm = 5000L) {

  if (!is.numeric(baseline) || !is.numeric(intervention))
    stop("Both 'baseline' and 'intervention' must be numeric vectors.")
  if (length(baseline) != length(intervention))
    stop("'baseline' and 'intervention' must have the same length (paired).")
  n <- length(baseline)
  if (n < 3L) stop("At least 3 paired observations required.")
  nPerm <- as.integer(nPerm)
  if (nPerm < 100L) stop("nPerm must be >= 100.")

  mo <- mean(baseline);   mi <- mean(intervention)
  so <- stats::sd(baseline); si <- stats::sd(intervention)
  seo <- so / sqrt(n);    sei <- si / sqrt(n)

  pair_diffs <- intervention - baseline
  obs_md     <- mean(pair_diffs)
  obs_stat   <- abs(obs_md)
  sd_diff    <- stats::sd(pair_diffs)
  d_paired   <- if (sd_diff > 1e-10) obs_md / sd_diff else 0

  # Vectorised sign-flip
  sign_matrix <- matrix(
    sample(c(-1L, 1L), n * nPerm, replace = TRUE),
    nrow = nPerm, ncol = n
  )
  perm_means <- as.numeric(sign_matrix %*% pair_diffs) / n
  perm_stats <- abs(perm_means)

  extreme_count <- sum(perm_stats >= obs_stat)

  list(
    mean_original            = mo,
    mean_other               = mi,
    sd_original              = so,
    sd_other                 = si,
    se_original              = seo,
    se_other                 = sei,
    ci_original              = 1.96 * seo,
    ci_other                 = 1.96 * sei,
    ciLower_original         = mo - 1.96 * seo,
    ciUpper_original         = mo + 1.96 * seo,
    ciLower_other            = mi - 1.96 * sei,
    ciUpper_other            = mi + 1.96 * sei,
    pooled_std_dev           = sd_diff,
    cohens_d                 = d_paired,
    p_value                  = (extreme_count + 1L) / (nPerm + 1L),
    n_permutations           = nPerm,
    permutation_distribution = perm_stats
  )
}
