# =============================================================================
# kimSelectionStats.R -- score a fitted network against a known one
#
# Used by the simulation studies and by anyone validating the estimator on data
# with a known generating structure. Everything is computed on the off-diagonal
# cells by default, because the diagonal (self-persistence) is almost never the
# target of edge selection and is detected in essentially every fit.
# =============================================================================

#' Score edge selection against a known transition matrix
#'
#' @param fit A `"dynaNIRA_kim"` fit from [fitKIM()] (or any object with `J` and
#'   `vars`).
#' @param truth The true `p x p` matrix J, with the same `dimnames` as the fit's
#'   variables.
#' @param off_diagonal_only Ignore the diagonal. Default `TRUE`.
#'
#' @return A one-row data frame with
#'   `TP`, `FP`, `FN`, `TN` (cell counts),
#'   `sensitivity` (share of true edges found),
#'   `specificity` (share of true zeros excluded),
#'   `precision` (share of selected cells that are true),
#'   `F1`, `sign_accuracy` (among correctly detected true edges, the share whose
#'   estimated sign matches), and `n_selected`.
#'
#' @details A cell counts as "selected" when `fit$J` is non-zero. A fit produced
#'   with `screen = "none"` has no zero cells, so its specificity is 0 and its
#'   sensitivity 1 by construction: it never claims sparsity. Comparing it with a
#'   screened fit answers "what does screening buy?", not "which fit is better
#'   specified".
#' @export
kimSelectionStats <- function(fit, truth, off_diagonal_only = TRUE) {
  if (is.null(fit$J) || is.null(fit$vars))
    stop("`fit` must carry `J` and `vars`.")
  Jhat <- fit$J
  if (!identical(dim(truth), dim(Jhat)))
    stop("`truth` must be the same p x p shape as fit$J.")
  vars <- fit$vars
  if (!is.null(dimnames(truth)) && !is.null(dimnames(Jhat)) &&
      !identical(dimnames(truth)[[1]], dimnames(Jhat)[[1]]))
    truth <- truth[rownames(Jhat), colnames(Jhat), drop = FALSE]

  S  <- (Jhat != 0)
  Tt <- (truth != 0)
  off <- !diag(TRUE, nrow(S))
  if (off_diagonal_only) { S <- S & off; Tt <- Tt & off }

  TP <- sum(S & Tt); FP <- sum(S & !Tt); FN <- sum(!S & Tt); TN <- sum(!S & !Tt)
  det <- S & Tt
  data.frame(
    TP = TP, FP = FP, FN = FN, TN = TN,
    sensitivity = TP / max(TP + FN, 1),
    specificity = TN / max(TN + FP, 1),
    precision   = TP / max(TP + FP, 1),
    F1          = 2 * TP / max(2 * TP + FP + FN, 1),
    sign_accuracy = if (any(det)) mean(sign(Jhat[det]) == sign(truth[det])) else NA_real_,
    n_selected  = sum(S),
    row.names   = NULL)
}
