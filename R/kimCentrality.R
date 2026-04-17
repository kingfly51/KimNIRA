# =============================================================================
# kimCentrality.R
# Compute and plot centrality indices for KIM networks
#
# Handles the fitKIM convention: J[i,j] = j -> i (column points to row)
# Internally transposes before passing to qgraph (which uses i -> j).
# =============================================================================

#' Compute and Plot Centrality Indices for KIM Networks
#'
#' @description
#' Computes centrality indices for both the temporal (directed, asymmetric)
#' network and the between-subject (undirected, partial correlation) network
#' from a fitted KIM object. Correctly handles the fitKIM convention where
#' \code{J[i,j]} represents the effect of node j (at time t) on node i
#' (at time t+1) -- i.e., column drives row -- which is the transpose of
#' qgraph's default convention.
#'
#' \strong{Temporal network indices (directed):}
#' \itemize{
#'   \item \strong{OutStrength}: sum of absolute outgoing edge weights from a
#'     node (how strongly a node drives others at the next time point).
#'     Equals \code{colSums(abs(J))} in the fitKIM convention.
#'     After transposing J for qgraph, this corresponds to
#'     \code{qgraph::centrality()$InDegree}.
#'   \item \strong{InStrength}: sum of absolute incoming edge weights to a
#'     node (how strongly a node is driven by others).
#'     Equals \code{rowSums(abs(J))} in the fitKIM convention.
#'     After transposing J for qgraph, this corresponds to
#'     \code{qgraph::centrality()$OutDegree}.
#'   \item \strong{OutEI (Out-Expected Influence)}: signed sum of outgoing
#'     edges (\code{colSums(J)}); positive = net activating driver.
#'   \item \strong{InEI (In-Expected Influence)}: signed sum of incoming
#'     edges (\code{rowSums(J)}).
#' }
#'
#' \strong{Between-subject network indices (undirected):}
#' Strength, Closeness, Betweenness, and Expected Influence computed via
#' \code{qgraph::centrality()}.
#'
#' @param fit A \code{"dynaNIRA_kim"} object returned by \code{\link{fitKIM}}.
#' @param plot Logical. Whether to draw centrality plots via
#'   \code{qgraph::centralityPlot}. Default \code{TRUE}.
#' @param scale Character. Scaling method for \code{qgraph::centralityPlot}.
#'   One of \code{"z-scores"} (default), \code{"raw"}, \code{"raw0"},
#'   \code{"relative"}.
#' @param order_by_temporal Character. Temporal-network index used to sort
#'   nodes in the plot. Default \code{"OutStrength"}.
#' @param order_by_between Character. Between-network index used to sort nodes
#'   in the plot. Default \code{"Strength"}.
#' @param include_temporal Character vector of temporal-network indices to
#'   display. Default:
#'   \code{c("InStrength","OutStrength","InExpectedInfluence","OutExpectedInfluence")}.
#' @param include_between Character vector of between-network indices to
#'   display. Default:
#'   \code{c("Strength","Closeness","Betweenness","ExpectedInfluence")}.
#' @param signed Logical. Passed to \code{qgraph::centralityPlot}. Whether
#'   to compute signed (Expected Influence) or unsigned (Strength) versions.
#'   Default \code{TRUE}.
#' @param use_ggplot Logical. If \code{TRUE}, additionally return custom
#'   \pkg{ggplot2} versions of the centrality plots. Default \code{FALSE}.
#' @param verbose Logical. Print progress messages. Default \code{TRUE}.
#'
#' @return Invisibly returns a named list:
#' \describe{
#'   \item{\code{temporal}}{Data frame of temporal-network centrality indices
#'     (InStrength, OutStrength, InEI, OutEI) for each node.}
#'   \item{\code{between}}{Data frame of between-network centrality indices
#'     (Strength, Closeness, Betweenness, EI) for each node. \code{NULL} if
#'     the between-subject network was not estimated.}
#'   \item{\code{temporal_raw}}{Full output of \code{qgraph::centrality()}
#'     applied to the transposed J matrix.}
#'   \item{\code{between_raw}}{Full output of \code{qgraph::centrality()}
#'     applied to the between-subject matrix. \code{NULL} if unavailable.}
#'   \item{\code{plots}}{A list with elements \code{temporal} and
#'     \code{between} containing \pkg{ggplot2} objects. Only populated when
#'     \code{use_ggplot = TRUE}.}
#' }
#'
#' @section Direction convention:
#' \code{fitKIM} stores the temporal network such that \code{J[i,j]} is the
#' directed effect of node j (at time t) on node i (at time t+1) -- column
#' is the source, row is the target.  \pkg{qgraph} assumes the opposite
#' (\code{M[i,j]} = i drives j), so this function passes \code{t(J)} to
#' \code{qgraph::centrality()}.
#'
#' After this transposition the mapping between qgraph labels and fitKIM
#' concepts is:
#' \tabular{ll}{
#'   \strong{qgraph label} \tab \strong{fitKIM meaning} \cr
#'   \code{InDegree}  \tab InStrength  = \code{rowSums(|J|)} (node is driven by others) \cr
#'   \code{OutDegree} \tab OutStrength = \code{colSums(|J|)} (node drives others) \cr
#' }
#' The returned data frame uses the \strong{fitKIM-sense} column names
#' (OutStrength, InStrength) and is verified against manual \code{colSums}/
#' \code{rowSums} of the original J matrix.
#'
#' @examples
#' \dontrun{
#' fit <- fitKIM(data = ema_long,
#'               vars = c("Dep","Loneliness","RIS","Stress","SA","Anx"),
#'               id = "BH", time = "day_num")
#'
#' # Default: qgraph plots + data frames
#' cent <- kimCentrality(fit)
#' cent$temporal    # data frame
#' cent$between
#'
#' # With custom ggplot2 output
#' cent <- kimCentrality(fit, use_ggplot = TRUE)
#' cent$plots$temporal
#' cent$plots$between
#' }
#'
#' @seealso \code{\link{fitKIM}}, \code{\link{kimPlotNetwork}},
#'   \code{\link{kimPlotAllNetworks}}
#' @importFrom stats reorder
#' @export
kimCentrality <- function(fit,
                          plot              = TRUE,
                          scale             = c("z-scores", "raw", "raw0", "relative"),
                          order_by_temporal = "OutStrength",
                          order_by_between  = "Strength",
                          include_temporal  = c("InStrength", "OutStrength",
                                                "InExpectedInfluence",
                                                "OutExpectedInfluence"),
                          include_between   = c("Strength", "Closeness",
                                                "Betweenness",
                                                "ExpectedInfluence"),
                          signed            = TRUE,
                          use_ggplot        = FALSE,
                          verbose           = TRUE) {

  # -- 0. Input checks ---------------------------------------------------------
  if (!requireNamespace("qgraph", quietly = TRUE))
    stop("Package 'qgraph' is required: install.packages('qgraph').")

  if (!inherits(fit, "dynaNIRA_kim"))
    warning("Input does not inherit from 'dynaNIRA_kim'; ",
            "proceeding but verify this is a fitKIM() result.")

  if (is.null(fit$temporal_network))
    stop("fit$temporal_network is NULL; cannot compute temporal centrality.")

  scale <- match.arg(scale)

  # -- 1. Temporal network ------------------------------------------------------
  # fitKIM convention: J[i,j] = effect of j (at t) on i (at t+1)
  #   => column j is the SOURCE, row i is the TARGET
  #   => OutStrength of node j = colSums(abs(J))[j]  (how much j drives others)
  #   => InStrength  of node i = rowSums(abs(J))[i]  (how much i is driven)
  #
  # After transposing to J_qg = t(J): J_qg[j, i] = J[i, j]
  #   => row j is the SOURCE, col i is the TARGET  (qgraph row->col convention)
  #
  # qgraph::centrality on J_qg:
  #   InDegree[i]  = sum of edges COMING IN  to i = colSums(abs(J_qg))[i]
  #                = rowSums(abs(J))[i]  = fitKIM InStrength[i]   (node is driven)
  #   OutDegree[i] = sum of edges GOING OUT from i = rowSums(abs(J_qg))[i]
  #                = colSums(abs(J))[i]  = fitKIM OutStrength[i]  (node drives others)
  J    <- fit$temporal_network
  vars <- if (!is.null(fit$vars)) fit$vars else rownames(J)
  if (is.null(rownames(J))) dimnames(J) <- list(vars, vars)

  J_qg <- t(J)    # J_qg[j, i] = J[i, j]; row = source, col = target

  cent_temp <- qgraph::centrality(J_qg)

  temporal_df <- data.frame(
    Node        = vars,
    # qgraph InDegree  = colSums(abs(J_qg)) = rowSums(abs(J)) = fitKIM InStrength
    InStrength  = cent_temp$InDegree,    # rowSums(abs(J)): node is driven by others
    # qgraph OutDegree = rowSums(abs(J_qg)) = colSums(abs(J)) = fitKIM OutStrength
    OutStrength = cent_temp$OutDegree,   # colSums(abs(J)): node drives others
    InEI        = rowSums(J),            # signed In-Expected Influence
    OutEI       = colSums(J),            # signed Out-Expected Influence
    stringsAsFactors = FALSE
  )

  # Direction sanity check
  stopifnot(
    isTRUE(all.equal(unname(temporal_df$InStrength),  unname(rowSums(abs(J))))),
    isTRUE(all.equal(unname(temporal_df$OutStrength), unname(colSums(abs(J)))))
  )

  if (verbose) {
    message("Temporal network centrality computed for ",
            nrow(temporal_df), " node(s).")
    message("  Convention: J[i,j] = j -> i (fitKIM).")
    message("  OutStrength = colSums(|J|): how strongly a node drives others.")
    message("  InStrength  = rowSums(|J|): how strongly a node is driven.")
  }

  # -- 2. Between-subject network -----------------------------------------------
  between_df  <- NULL
  cent_btw    <- NULL
  between_mat <- fit$between_network

  if (is.null(between_mat)) {
    if (verbose)
      message("fit$between_network is NULL; skipping between-network centrality.")
  } else {

    # Replace NAs (cells involving nodes with zero between-person variance,
    # which were excluded from EBICglasso and back-filled with NA rather than 0)
    # with 0 for centrality computation, but report which nodes were affected.
    na_cells <- sum(is.na(between_mat))
    if (na_cells > 0L) {
      dropped <- fit$between_dropped_nodes
      if (!is.null(dropped) && length(dropped) > 0L && verbose)
        message("Between network: ", length(dropped),
                " node(s) excluded from EBICglasso (zero between-person ",
                "variance): ", paste(dropped, collapse = ", "),
                ". Their ", na_cells, " NA cell(s) set to 0 for centrality.")
      between_mat[is.na(between_mat)] <- 0
    }

    btw_vars <- if (!is.null(rownames(between_mat))) rownames(between_mat) else vars
    cent_btw <- qgraph::centrality(between_mat)

    between_df <- data.frame(
      Node        = btw_vars,
      Strength    = cent_btw$OutDegree,   # undirected: InDegree = OutDegree
      Closeness   = cent_btw$Closeness,
      Betweenness = cent_btw$Betweenness,
      EI          = colSums(between_mat), # signed Expected Influence
      stringsAsFactors = FALSE
    )

    if (verbose)
      message("Between-subject network centrality computed for ",
              nrow(between_df), " node(s).")
  }

  # -- 3. qgraph centrality plots -----------------------------------------------
  if (plot) {
    tryCatch(
      print(qgraph::centralityPlot(
        list(Temporal = J_qg),
        include = include_temporal,
        scale   = scale,
        orderBy = order_by_temporal,
        signed  = signed
      )),
      error = function(e)
        warning("Temporal centralityPlot failed: ", e$message)
    )

    if (!is.null(between_mat)) {
      tryCatch(
        print(qgraph::centralityPlot(
          list(Between = between_mat),
          include = include_between,
          scale   = scale,
          orderBy = order_by_between,
          signed  = signed
        )),
        error = function(e)
          warning("Between centralityPlot failed: ", e$message)
      )
    }
  }

  # -- 4. Optional ggplot2 versions ---------------------------------------------
  gg_plots <- NULL
  if (use_ggplot) {
    need_pkgs    <- c("ggplot2", "dplyr", "tidyr")
    missing_pkgs <- need_pkgs[
      !vapply(need_pkgs, requireNamespace, logical(1), quietly = TRUE)
    ]
    if (length(missing_pkgs) > 0L) {
      warning("Custom ggplot2 output requires: ",
              paste(missing_pkgs, collapse = ", "), ". Skipping use_ggplot.")
    } else {
      gg_plots <- list(
        temporal = .kim_gg_temporal(temporal_df),
        between  = if (!is.null(between_df)) .kim_gg_between(between_df) else NULL
      )
      if (plot) {
        print(gg_plots$temporal)
        if (!is.null(gg_plots$between)) print(gg_plots$between)
      }
    }
  }

  # -- 5. Return ----------------------------------------------------------------
  out <- list(
    temporal     = temporal_df,
    between      = between_df,
    temporal_raw = cent_temp,
    between_raw  = cent_btw,
    plots        = gg_plots
  )
  invisible(out)
}


# -----------------------------------------------------------------------------
# Internal helpers: custom ggplot2 centrality panels
# -----------------------------------------------------------------------------

#' @keywords internal
.kim_gg_temporal <- function(df) {
  long <- tidyr::pivot_longer(
    df[, c("Node", "InStrength", "OutStrength")],
    cols      = c("InStrength", "OutStrength"),
    names_to  = "Index",
    values_to = "Value"
  )

  ggplot2::ggplot(
    long,
    ggplot2::aes(
      x    = stats::reorder(.data$Node, .data$Value),
      y    = .data$Value,
      fill = .data$Index
    )
  ) +
    ggplot2::geom_col(
      position = ggplot2::position_dodge(width = 0.7),
      width    = 0.6
    ) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title    = "Temporal Network Centrality",
      subtitle = "OutStrength = node drives others  |  InStrength = node is driven",
      x        = NULL,
      y        = "Centrality (sum of |edge weights|)",
      fill     = NULL
    ) +
    ggplot2::scale_fill_manual(
      values = c("InStrength" = "#2C7FB8", "OutStrength" = "#E76F51")
    ) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      legend.position    = "top",
      panel.grid.major.y = ggplot2::element_blank()
    )
}

#' @keywords internal
.kim_gg_between <- function(df) {
  # z-score each index column before pivoting
  num_cols <- setdiff(names(df), "Node")
  z_df     <- df
  z_df[num_cols] <- lapply(df[num_cols], function(x) as.numeric(base::scale(x)))

  long <- tidyr::pivot_longer(
    z_df,
    cols      = -"Node",
    names_to  = "Index",
    values_to = "Value"
  )

  ggplot2::ggplot(
    long,
    ggplot2::aes(
      x = .data$Value,
      y = stats::reorder(.data$Node, .data$Value)
    )
  ) +
    ggplot2::geom_line(
      ggplot2::aes(group = .data$Node),
      colour  = "grey75",
      linewidth = 0.4
    ) +
    ggplot2::geom_point(size = 2.5, colour = "#2A9D8F") +
    ggplot2::facet_wrap(~ Index, nrow = 1) +
    ggplot2::labs(
      title = "Between-Subject Network Centrality (z-scores)",
      x     = NULL,
      y     = NULL
    ) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank()
    )
}
