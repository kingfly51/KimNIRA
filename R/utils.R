# =============================================================================
# utils.R -- Package-level utilities
# =============================================================================

#' Null-coalescing operator
#'
#' Returns \code{a} if it is not \code{NULL}, otherwise returns \code{b}.
#' @keywords internal
#' @noRd
`%||%` <- function(a, b) if (!is.null(a)) a else b

#' @importFrom stats binomial coef complete.cases cor median na.omit predict runif sd setNames var
#' @noRd
NULL

utils::globalVariables(".data")
