# =============================================================================
# utils.R -- Package-level utilities
# =============================================================================

#' Null-coalescing operator
#'
#' Returns \code{a} if it is not \code{NULL}, otherwise returns \code{b}.
#' @keywords internal
`%||%` <- function(a, b) if (!is.null(a)) a else b
