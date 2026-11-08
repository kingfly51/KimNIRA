# =============================================================================
# kimPlotAllNetworks.R  ----  Side-by-side network plots (temporal + between-subject)
#
# Adapted for the "dynaNIRA_kim" object returned by fitKIM().
# Dependencies: qgraph (required), kimPlotNetwork() (same package)
# =============================================================================


#' Plot Temporal and Between-Subject Networks Side by Side
#'
#' @description
#' Draws the temporal network and the between-subject network on the same
#' canvas in two side-by-side panels.
#' Node positions are fixed (shared layout) to facilitate visual comparison
#' across networks.
#'
#' Notes:
#' \itemize{
#'   \item If the between-subject network could not be estimated
#'     (\code{between_network = NULL}), the right panel is left blank with
#'     an informative message.
#'   \item Node positions are determined by the spring layout of the first
#'     (temporal) network and reused for subsequent networks.
#' }
#'
#' @param kim_fit Object returned by fitKIM().
#' @param layout Character string, initial layout algorithm. Default \code{"spring"}.
#' @param nonsig Character, zero-edge handling: \code{"hide"} (default),
#'   \code{"dashed"}, \code{"show"}.
#' @param threshold Numeric. Edges with absolute weight below this value are
#'   treated as zero. Default \code{0}.
#' @param edge.labels Logical, whether to display edge weight values. Default \code{FALSE}.
#' @param labels Character vector or NULL, node labels. Uses vars when NULL.
#' @param node.size Numeric, node size. Default \code{8}.
#' @param title_temporal Character, title for the temporal network panel.
#' @param title_between Character, title for the between-subject network panel.
#' @param ... Additional arguments passed to \code{qgraph::qgraph()}.
#'
#' @return Invisibly returns a named list of qgraph objects.
#' @importFrom graphics plot.new
#' @export
kimPlotAllNetworks <- function(kim_fit,
                               layout                = "spring",
                               nonsig                = c("hide", "dashed", "show"),
                               threshold             = 0,
                               edge.labels           = FALSE,
                               labels                = NULL,
                               node.size             = 8,
                               title_temporal        = "Temporal",
                               title_between         = "Between-subjects",
                               ...) {

  # -- Input validation --------------------------------------------------------
  if (!inherits(kim_fit, "dynaNIRA_kim"))
    stop("`kim_fit` must be produced by fitKIM().")
  if (!requireNamespace("qgraph", quietly = TRUE))
    stop("Package 'qgraph' is required: install.packages('qgraph')")

  nonsig <- match.arg(nonsig)

  node_labels <- if (!is.null(labels)) labels else kim_fit$vars

  # -- Set up canvas: 2 columns ------------------------------------------------
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)

  graphics::layout(matrix(1:2, nrow = 1))

  results      <- list()
  fixed_layout <- NULL

  # -- Temporal network --------------------------------------------------------
  g <- kimPlotNetwork(
    kim_fit     = kim_fit,
    type        = "temporal",
    layout      = layout,
    nonsig      = nonsig,
    threshold   = threshold,
    edge.labels = edge.labels,
    title       = title_temporal,
    labels      = node_labels,
    node.size   = node.size,
    ...
  )
  results$temporal <- g
  fixed_layout <- g$layout

  # -- Between-subject network -------------------------------------------------
  if (is.null(kim_fit$between_network)) {
    plot.new()
    graphics::text(0.5, 0.5,
                   "Between-subjects network\nnot available",
                   cex = 1.1, col = "grey40")
    results$between <- NULL
  } else {
    g <- kimPlotNetwork(
      kim_fit     = kim_fit,
      type        = "between",
      layout      = fixed_layout,
      nonsig      = nonsig,
      threshold   = threshold,
      edge.labels = edge.labels,
      title       = title_between,
      labels      = node_labels,
      node.size   = node.size,
      ...
    )
    results$between <- g
  }

  graphics::layout(1)   # Restore single-panel layout

  invisible(results)
}
