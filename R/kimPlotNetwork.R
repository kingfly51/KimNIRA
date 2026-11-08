# =============================================================================
# kimPlotNetwork.R  ----  Plot a single KIM network
#
# Adapted for the "dynaNIRA_kim" object returned by fitKIM().
#
# Key fields in the fitKIM return structure:
#   J / temporal_network        Temporal interaction matrix (directed;
#                               J[i,j] = influence of j on i)
#   h / thresholds              Intercept vector
#   between_network             Between-subject network (EBICglasso partial
#                               correlation matrix)
#   vars                        Node name vector
#
# Dependency: qgraph (required)
# =============================================================================


#' Plot a Single KIM Network
#'
#' @description
#' Uses qgraph to plot one of the two networks produced by fitKIM:
#' \itemize{
#'   \item \code{"temporal"}: Temporal network (J matrix), directed graph,
#'     arrows point from j to i.
#'   \item \code{"between"}: Between-subject network (EBICglasso), undirected.
#' }
#'
#' @param kim_fit Object returned by fitKIM() (class contains "dynaNIRA_kim").
#' @param type Character, network type: \code{"temporal"} or \code{"between"}.
#'   Default \code{"temporal"}.
#' @param layout Character or matrix, node layout algorithm. Default \code{"spring"}.
#' @param nonsig Character, handling of zero/near-zero edge weights:
#'   \code{"hide"} (remove; recommended), \code{"dashed"} (dashed lines),
#'   \code{"show"} (display). Default \code{"hide"}.
#' @param threshold Numeric. Edges with absolute weight below this value are
#'   treated as zero. Default \code{0}.
#' @param edge.labels Logical, whether to display edge weight values. Default \code{FALSE}.
#' @param edge.label.cex Numeric, font size for edge labels. Default \code{0.6}.
#' @param title Character or NULL, plot title. Auto-generated when NULL.
#' @param title.cex Numeric, title font size. Default \code{1.2}.
#' @param labels Character vector or NULL, node labels. Uses kim_fit$vars when NULL.
#' @param label.cex Numeric, node label font size. Default \code{1.2}.
#' @param node.size Numeric, node size. Default \code{8}.
#' @param node.color Character, node fill colour. Default \code{"#BEB8DC"}.
#' @param ... Additional arguments passed to \code{qgraph::qgraph()}.
#'
#' @return Invisibly returns the qgraph object (includes \code{$layout}
#'   coordinates that can be extracted and passed to the next network plot).
#' @export
kimPlotNetwork <- function(kim_fit,
                           type           = c("temporal", "between"),
                           layout         = "spring",
                           nonsig         = c("hide", "dashed", "show"),
                           threshold      = 0,
                           edge.labels    = FALSE,
                           edge.label.cex = 0.6,
                           title          = NULL,
                           title.cex      = 1.2,
                           labels         = NULL,
                           label.cex      = 1.2,
                           node.size      = 8,
                           node.color     = "#BEB8DC",
                           ...) {

  # -- Input validation --------------------------------------------------------
  if (!inherits(kim_fit, "dynaNIRA_kim"))
    stop("`kim_fit` must be produced by fitKIM() (class must contain 'dynaNIRA_kim').")
  if (!requireNamespace("qgraph", quietly = TRUE))
    stop("Package 'qgraph' is required: install.packages('qgraph')")

  type   <- match.arg(type)
  nonsig <- match.arg(nonsig)

  vars   <- kim_fit$vars
  n_vars <- length(vars)

  # Node labels
  node_labels <- if (!is.null(labels)) labels else vars

  # -- Select matrix -----------------------------------------------------------
  if (type == "temporal") {

    # J[i,j] = influence of node j at t on node i at t+1
    # qgraph convention: input[i,j] = edge weight from node i to node j
    # So we pass t(J): t(J)[j,i] = J[i,j] -> draws arrow j -> i
    J <- kim_fit$J
    if (is.null(J)) J <- kim_fit$temporal_network
    if (is.null(J)) stop("Cannot find temporal matrix (J / temporal_network) in kim_fit.")
    NET      <- t(J)
    directed <- TRUE
    auto_title <- "Temporal Network"

  } else {  # "between"

    NET <- kim_fit$between_network
    if (is.null(NET)) {
      stop("Between-subject network is NULL. ",
           "Check whether fitKIM() successfully estimated it (requires 'qgraph').")
    }
    directed   <- FALSE
    auto_title <- "Between-subjects Network"
  }

  rownames(NET) <- vars
  colnames(NET) <- vars

  # -- Threshold truncation ----------------------------------------------------
  if (threshold > 0) {
    NET[abs(NET) < threshold] <- 0
  }

  # -- Handle nonsig -----------------------------------------------------------
  if (nonsig == "dashed") {
    lty <- ifelse(NET == 0, 2, 1)
  } else {
    lty <- 1
  }

  # -- Plot title --------------------------------------------------------------
  plot_title <- if (!is.null(title)) title else auto_title

  # -- Call qgraph -------------------------------------------------------------
  g <- qgraph::qgraph(
    input          = NET,
    directed       = directed,
    layout         = layout,
    posCol         = "#2171B5",
    negCol         = "#CB181D",
    fade           = TRUE,
    lty            = lty,
    labels         = node_labels,
    label.cex      = label.cex,
    vsize          = node.size,
    edge.labels    = edge.labels,
    edge.label.cex = edge.label.cex,
    title          = plot_title,
    title.cex      = title.cex,
    mar            = c(3, 3, 5, 3),
    color          = node.color,
    ...
  )

  invisible(g)
}
