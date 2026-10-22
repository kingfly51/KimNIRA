#' Plot NIRA intervention effects
#'
#' @description
#' Plots the projected mean activation for the baseline and for each intervention
#' target.
#'
#' \itemize{
#'   \item \strong{Exact results} are shown as exact point estimates with no
#'     interval. There is no numerical uncertainty to display, so no significance
#'     stars are drawn.
#'   \item \strong{Monte Carlo results} show a normal mean interval based on the
#'     retained trajectory scores. A star marks an effect whose paired-difference
#'     test versus no change is significant; stars are Benjamini-Hochberg
#'     adjusted across nodes. This reflects Monte Carlo precision, not sampling
#'     uncertainty. No equality mark is ever added for comparisons with the
#'     top-ranked target.
#' }
#'
#' The function accepts only a \code{runKIMNIRA()} result. No bootstrap object
#' and no \code{stars} argument are needed.
#'
#' @param statresult A \code{"kim_nira_result"} from \code{\link{runKIMNIRA}}.
#' @param perturbation_type \code{"alleviating"} or \code{"aggravating"};
#'   inferred from the result when omitted.
#' @param node_labels Optional named character vector of display labels,
#'   including an entry for \code{"baseline"}.
#' @param level Confidence level for the Monte Carlo interval. Ignored for
#'   exact results.
#' @param effect Plot the total NIRA effect (\code{"total"}) or the mean change
#'   among non-target nodes (\code{"spillover"}). Spillover plots are point
#'   estimates only.
#' @param p_adjust Multiple-comparison adjustment for the Monte Carlo
#'   significance stars, passed to \code{\link[stats]{p.adjust}}. Default
#'   \code{"BH"}. Use \code{"none"} for unadjusted p-values.
#'
#' @return A \pkg{ggplot2} object.
#' @export
kimPlotNIRA <- function(statresult,
                        perturbation_type = NULL,
                        node_labels = NULL,
                        level = 0.95,
                        effect = c("total", "spillover"),
                        p_adjust = "BH") {

  if (!inherits(statresult, "kim_nira_result"))
    stop("`statresult` must be produced by runKIMNIRA().", call. = FALSE)

  effect <- match.arg(effect)
  if (length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1)
    stop("`level` must be between 0 and 1.", call. = FALSE)
  p_adjust <- match.arg(p_adjust, stats::p.adjust.methods)

  is_mc <- identical(statresult$method, "mc")
  direction <- statresult$params$intervention_type
  if (is.null(perturbation_type)) perturbation_type <- direction
  perturbation_type <- match.arg(perturbation_type,
                                 c("alleviating", "aggravating"))

  vars <- statresult$kim_fit$vars
  p    <- length(vars)

  # ------------------------------------------------------------------ data ---
  if (effect == "spillover") {
    spill <- statresult$summary$spillover_mean_diff[
      match(vars, statresult$summary$node)]
    if (length(spill) != p || any(!is.finite(spill)))
      stop("This result does not contain finite spillover effects.",
           call. = FALSE)
    pd <- data.frame(
      node    = c("baseline", vars),
      mean    = c(0, spill),
      ciLower = NA_real_,
      ciUpper = NA_real_,
      star    = "",
      stringsAsFactors = FALSE
    )
  } else {
    pd <- data.frame(
      node    = c("baseline", vars),
      mean    = c(statresult$base_mean,
                  statresult$base_mean + statresult$effects[vars]),
      ciLower = NA_real_,
      ciUpper = NA_real_,
      star    = "",
      stringsAsFactors = FALSE
    )

    if (is_mc) {
      # Visual interval for the projected activation level.
      scores <- c(list(statresult$baseline_scores),
                  statresult$intervention_scores[vars])
      bounds <- vapply(scores, function(x) {
        if (!is.numeric(x) || length(x) < 2L || any(!is.finite(x)))
          stop("At least two finite trajectory scores are required per condition.",
               call. = FALSE)
        mean(x) + c(-1, 1) * stats::qnorm((1 + level) / 2) *
          stats::sd(x) / sqrt(length(x))
      }, numeric(2))
      pd$ciLower <- bounds[1, ]
      pd$ciUpper <- bounds[2, ]

      # Significance of the intervention effect versus no change, from the
      # paired baseline/intervention trajectory scores.
      pv <- vapply(seq_len(p), function(j) {
        d <- statresult$intervention_scores[[j]] - statresult$baseline_scores
        if (length(d) < 2L || any(!is.finite(d))) return(NA_real_)
        se <- stats::sd(d) / sqrt(length(d))
        if (!is.finite(se) || se <= 0) return(NA_real_)
        2 * stats::pnorm(-abs(mean(d) / se))
      }, numeric(1))
      pv_adj <- stats::p.adjust(pv, method = p_adjust)
      # 显著性星号：p < .05 一颗，p < .01 两颗，p < .001 三颗（校正后）。
      pd$star[-1] <- .kim_star_label(pv_adj)
    }
  }

  # ----------------------------------------------------------------- order ---
  ord <- order(pd$mean[pd$node != "baseline"],
               decreasing = perturbation_type == "aggravating")
  pd <- rbind(
    pd[pd$node == "baseline", , drop = FALSE],
    pd[pd$node != "baseline", , drop = FALSE][ord, , drop = FALSE]
  )
  pd$newid <- seq_len(nrow(pd)) - 1L

  labels <- pd$node
  if (!is.null(node_labels)) {
    if (is.null(names(node_labels)))
      stop("`node_labels` must be a named vector.", call. = FALSE)
    hit <- match(pd$node, names(node_labels))
    ok  <- !is.na(hit)
    labels[ok] <- unname(node_labels[hit[ok]])
  }

  pd$star_y <- ifelse(
    is.na(pd$ciUpper), pd$mean,
    ifelse(pd$mean >= 0, pd$ciUpper, pd$ciLower)
  )

  # ------------------------------------------------------------------ plot ---
  gg <- ggplot2::ggplot(
    pd, ggplot2::aes(x = .data$newid, y = .data$mean, group = 1)
  ) +
    ggplot2::geom_line(color = "steelblue", linewidth = 1) +
    ggplot2::geom_point(size = 3, color = "tomato")

  if (is_mc) {
    gg <- gg + ggplot2::geom_errorbar(
      ggplot2::aes(ymin = .data$ciLower, ymax = .data$ciUpper),
      width = 0.2, color = "gray40"
    )
  }

  labelled <- nzchar(pd$star)
  if (any(labelled)) {
    up <- labelled & !is.na(pd$mean) & pd$mean >= 0
    dn <- labelled & !is.na(pd$mean) & pd$mean <  0
    if (any(up)) gg <- gg + ggplot2::geom_text(
      data = pd[up, , drop = FALSE],
      ggplot2::aes(x = .data$newid, y = .data$star_y, label = .data$star),
      vjust = -0.7, size = 4
    )
    if (any(dn)) gg <- gg + ggplot2::geom_text(
      data = pd[dn, , drop = FALSE],
      ggplot2::aes(x = .data$newid, y = .data$star_y, label = .data$star),
      vjust = 1.7, size = 4
    )
  }

  gg <- gg +
    ggplot2::scale_x_continuous(breaks = pd$newid, labels = labels) +
    ggplot2::labs(
      x = "Intervention Target",
      y = if (effect == "spillover") {
        "Mean Change in Non-target Activation"
      } else {
        "Mean Projected Activation Rate"
      },
      title = ""
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.title.x       = ggplot2::element_text(face = "bold", size = 12),
      axis.title.y       = ggplot2::element_text(face = "bold", size = 12),
      axis.text.x        = ggplot2::element_text(angle = 45, hjust = 1,
                                                face = "bold", size = 10),
      axis.text.y        = ggplot2::element_text(face = "bold", size = 10),
      panel.grid.major.x = ggplot2::element_blank()
    )

  gg
}
