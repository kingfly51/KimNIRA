# =============================================================================
# kimPlotTrajectory.R  ----  KIM NodeIdentifyR Algorithm (NIRA) Longitudinal Intervention Trajectory Plots
#
# Visualises the longitudinal propagation of a simulated intervention across
# multiple time steps. Supports:
#   (A) Summary trajectory  (baseline vs intervention mean activation rate)
#   (B) Node-level trajectories  (per-node activation probability over time)
#   (C) Difference trajectory  (intervention - baseline at each step)
#   (D) Multi-target panel  (side-by-side comparison of different targets)
#   (E) Spaghetti plot  (individual-level trajectories showing variability)
#
# Dependency: ggplot2 (required)
# Input: "kim_nira_result" object returned by runKIMNIRA()
#
# Changes vs. previous version:
#   - All plot subtypes: title, subtitle, and caption removed entirely;
#     only axis titles (x / y) are shown.
#   - .plot_diff: y-axis label direction now reflects intervention_type
#     ("alleviating" -> negative = alleviating; "aggravating" -> positive = aggravating).
#   - %||% defined at file level so this file can be sourced standalone.
# =============================================================================




# -----------------------------------------------------------------------------
# Internal helper: stepwise trajectory simulation (records each step)
# -----------------------------------------------------------------------------

#' Stepwise Trajectory Simulation (Internal)
#'
#' Similar to kim_trajectory() but records node activation rates and
#' summary scores at every step.
#'
#' @keywords internal
.kim_trajectory_stepwise <- function(init_states, J, h, n_steps,
                                     target_node = NULL, delta = 0,
                                     U_list = NULL, n_traces = 0L) {
  n_sim  <- nrow(init_states)
  n_vars <- ncol(init_states)
  vars   <- colnames(init_states)

  h_eff <- h
  if (!is.null(target_node) && delta != 0)
    h_eff[target_node] <- h[target_node] + delta

  # Pre-allocate
  node_probs     <- matrix(NA_real_, nrow = n_steps + 1L, ncol = n_vars)
  colnames(node_probs) <- vars
  summary_scores <- numeric(n_steps + 1L)
  q025           <- numeric(n_steps + 1L)
  q975           <- numeric(n_steps + 1L)

  # Individual traces (optional)
  trace_mat      <- NULL
  n_trace_actual <- 0L
  if (n_traces > 0L) {
    n_trace_actual <- min(as.integer(n_traces), n_sim)
    trace_mat <- matrix(NA_real_, nrow = n_steps + 1L, ncol = n_trace_actual)
  }

  # Step 0
  s <- init_states
  storage.mode(s) <- "double"
  row_scores           <- rowSums(s) / n_vars
  node_probs[1L, ]     <- colMeans(s)
  summary_scores[1L]   <- mean(row_scores)
  q025[1L]             <- stats::quantile(row_scores, 0.025, names = FALSE)
  q975[1L]             <- stats::quantile(row_scores, 0.975, names = FALSE)
  if (n_trace_actual > 0L)
    trace_mat[1L, ] <- row_scores[seq_len(n_trace_actual)]

  # Step-wise simulation
  for (step in seq_len(n_steps)) {
    U_step <- if (!is.null(U_list)) U_list[[step]] else NULL
    s          <- kim_step(s, J, h_eff, U_step = U_step)
    row_scores <- rowSums(s) / n_vars
    node_probs[step + 1L, ]   <- colMeans(s)
    summary_scores[step + 1L] <- mean(row_scores)
    q025[step + 1L]           <- stats::quantile(row_scores, 0.025, names = FALSE)
    q975[step + 1L]           <- stats::quantile(row_scores, 0.975, names = FALSE)
    if (n_trace_actual > 0L)
      trace_mat[step + 1L, ] <- row_scores[seq_len(n_trace_actual)]
  }

  list(node_probs = node_probs,
       summary    = summary_scores,
       q025       = q025,
       q975       = q975,
       traces     = trace_mat)
}


# -----------------------------------------------------------------------------
# Internal helper: run paired stepwise simulation
# -----------------------------------------------------------------------------

#' @keywords internal
.kim_run_stepwise <- function(nira_result, target_nodes, n_traces = 0L) {

  kim_fit <- nira_result$kim_fit
  p       <- nira_result$params

  J       <- kim_fit$J
  h       <- kim_fit$h
  vars    <- kim_fit$vars
  n_vars  <- length(vars)
  n_sim   <- p$n_sim
  n_steps <- p$n_steps
  seed    <- p$seed
  delta   <- nira_result$delta

  # Sample shared initial states
  init_states <- kim_init_states(
    h     = h,
    n_sim = n_sim,
    vars  = vars,
    seed  = seed
  )

  # Pre-generate shared random numbers (paired design)
  U_list <- .kim_pregenerate_U(
    n_sim   = n_sim,
    n_vars  = n_vars,
    n_steps = n_steps,
    seed    = seed + 1000L
  )

  # Baseline simulation
  baseline <- .kim_trajectory_stepwise(
    init_states = init_states, J = J, h = h,
    n_steps     = n_steps,
    target_node = NULL, delta = 0,
    U_list      = U_list, n_traces = n_traces
  )

  # Intervention simulations
  interventions <- list()
  for (j in target_nodes) {
    interventions[[vars[j]]] <- .kim_trajectory_stepwise(
      init_states = init_states, J = J, h = h,
      n_steps     = n_steps,
      target_node = j, delta = delta,
      U_list      = U_list, n_traces = n_traces
    )
  }

  list(baseline      = baseline,
       interventions = interventions,
       vars          = vars,
       n_steps       = n_steps,
       delta         = delta)
}


# -----------------------------------------------------------------------------
# Internal helper: build plotting data frames
# -----------------------------------------------------------------------------

#' @keywords internal
.kim_build_plot_data <- function(stepwise_result, target_node_name) {

  vars    <- stepwise_result$vars
  n_steps <- stepwise_result$n_steps
  bl      <- stepwise_result$baseline
  iv      <- stepwise_result$interventions[[target_node_name]]
  steps   <- 0:n_steps

  if (is.null(iv))
    stop("Node '", target_node_name, "' not found in intervention results. ",
         "Available nodes: ", paste(names(stepwise_result$interventions), collapse = ", "),
         ". Check the target argument or ensure the node was included in the simulation.")

  # Summary trajectory
  df_summary <- data.frame(
    step      = rep(steps, 2L),
    score     = c(bl$summary, iv$summary),
    condition = rep(c("Baseline", "Intervention"), each = n_steps + 1L),
    stringsAsFactors = FALSE
  )

  # Node-level trajectories
  df_node_list <- lapply(seq_along(vars), function(i) {
    data.frame(
      step      = rep(steps, 2L),
      prob      = c(bl$node_probs[, i], iv$node_probs[, i]),
      node      = vars[i],
      condition = rep(c("Baseline", "Intervention"), each = n_steps + 1L),
      stringsAsFactors = FALSE
    )
  })
  df_node <- do.call(rbind, df_node_list)

  # Difference trajectories
  df_diff_list <- lapply(seq_along(vars), function(i) {
    data.frame(
      step = steps,
      diff = iv$node_probs[, i] - bl$node_probs[, i],
      node = vars[i],
      stringsAsFactors = FALSE
    )
  })
  df_diff <- do.call(rbind, df_diff_list)

  # Summary score difference
  df_summary_diff <- data.frame(
    step = steps,
    diff = iv$summary - bl$summary,
    stringsAsFactors = FALSE
  )

  list(summary      = df_summary,
       node         = df_node,
       diff         = df_diff,
       summary_diff = df_summary_diff)
}


# -----------------------------------------------------------------------------
# Main plotting function
# -----------------------------------------------------------------------------

#' Plot KIM NodeIdentifyR Algorithm (NIRA) Longitudinal Intervention Trajectories
#'
#' @export
kimPlotTrajectory <- function(nira_result,
                              type             = c("summary", "node", "diff",
                                                   "panel", "spaghetti"),
                              target           = NULL,
                              top_n            = 3L,
                              n_traces         = 50L,
                              ribbon           = TRUE,
                              ci_band          = FALSE,
                              highlight_target = TRUE,
                              node_labels      = NULL,
                              colors           = NULL,
                              base_size        = 12,
                              line_width       = 1,
                              verbose          = FALSE) {

  if (!inherits(nira_result, "kim_nira_result"))
    stop("`nira_result` must be produced by runKIMNIRA() (class = 'kim_nira_result').")
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required: install.packages('ggplot2')")

  type <- match.arg(type)

  vars   <- nira_result$kim_fit$vars
  n_vars <- length(vars)
  p      <- nira_result$params
  ranked <- nira_result$ranked_targets

  # -- Resolve target ---------------------------------------------------------
  if (is.null(target)) {
    target_name <- ranked[1L]
  } else if (is.character(target)) {
    if (!target %in% vars)
      stop("target '", target, "' is not in the node list: ",
           paste(vars, collapse = ", "))
    target_name <- target
  } else if (is.numeric(target)) {
    idx <- as.integer(target)
    if (idx < 1L || idx > n_vars)
      stop("target index ", idx, " is out of range [1, ", n_vars, "].")
    target_name <- vars[idx]
  } else {
    stop("target must be a node name (character) or node index (integer).")
  }

  # -- Colour palette ----------------------------------------------------------
  default_colors <- c(
    "#E64B35", "#4DBBD5", "#00A087", "#3C5488",
    "#F39B7F", "#8491B4", "#91D1C2", "#DC0000",
    "#7E6148", "#B09C85", "#E377C2", "#BCBD22"
  )
  if (is.null(colors)) {
    if (n_vars <= length(default_colors)) {
      colors <- default_colors[seq_len(n_vars)]
    } else {
      colors <- c(default_colors,
                  grDevices::hcl.colors(n_vars - length(default_colors),
                                        palette = "Dark 3"))
    }
  } else if (length(colors) < n_vars) {
    warning("Length of colors (", length(colors), ") < number of nodes (",
            n_vars, "); padding with default colours.")
    colors <- c(colors, default_colors)[seq_len(n_vars)]
  }
  names(colors) <- vars

  # -- Node label mapping ------------------------------------------------------
  label_map <- stats::setNames(vars, vars)
  if (!is.null(node_labels)) {
    matched <- intersect(names(node_labels), vars)
    label_map[matched] <- node_labels[matched]
  }

  # -- Hard cap on n_traces ----------------------------------------------------
  if (type == "spaghetti") {
    n_traces <- as.integer(n_traces)
    if (n_traces > 200L) {
      warning("n_traces exceeds the upper limit of 200; truncating.")
      n_traces <- 200L
    }
  }

  # -- Determine which nodes to simulate --------------------------------------
  if (type == "panel") {
    top_n <- min(as.integer(top_n), n_vars)
    target_indices <- match(ranked[seq_len(top_n)], vars)
    need_traces <- 0L
  } else if (type == "spaghetti") {
    target_indices <- match(target_name, vars)
    need_traces <- n_traces
  } else {
    target_indices <- match(target_name, vars)
    need_traces <- 0L
  }

  if (verbose)
    message("Running stepwise simulation (", length(target_indices),
            " intervention condition(s))...")

  sw <- .kim_run_stepwise(nira_result,
                          target_nodes = target_indices,
                          n_traces     = need_traces)

  if (verbose) message("Simulation complete; generating plot...")

  # -- Retrieve intervention_type for direction-aware labels ------------------
  intervention_type <- nira_result$params$intervention_type %||% "alleviating"

  # -- Dispatch to plot subtype ------------------------------------------------
  gg <- switch(type,
               summary   = .plot_summary(sw, target_name, label_map, colors,
                                         ribbon, ci_band, line_width, base_size),
               node      = .plot_node(sw, target_name, label_map, colors,
                                      line_width, base_size),
               diff      = .plot_diff(sw, target_name, label_map, colors,
                                      highlight_target, line_width, base_size,
                                      intervention_type),
               panel     = .plot_panel(sw, ranked, top_n, label_map, colors,
                                       ribbon, line_width, base_size),
               spaghetti = .plot_spaghetti(sw, target_name, label_map, colors,
                                           line_width, base_size)
  )

  if (verbose) message("Plot complete.")
  invisible(gg)
}


# -----------------------------------------------------------------------------
# Shared theme  ---- no title / subtitle / caption; clean axis titles only
# -----------------------------------------------------------------------------

#' @keywords internal
.kim_base_theme <- function(base_size = 12) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      # Remove title, subtitle, caption entirely
      plot.title       = ggplot2::element_blank(),
      plot.subtitle    = ggplot2::element_blank(),
      plot.caption     = ggplot2::element_blank(),
      # Grid
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(color = "grey92", linewidth = 0.4),
      # Facet strip
      strip.text       = ggplot2::element_text(face = "bold", size = base_size - 1),
      # Legend
      legend.position  = "bottom",
      legend.title     = ggplot2::element_blank(),
      legend.text      = ggplot2::element_text(size = base_size - 2),
      # Margins & axes
      plot.margin      = ggplot2::margin(15, 15, 10, 15),
      axis.title       = ggplot2::element_text(size = base_size - 1, color = "grey30"),
      axis.text        = ggplot2::element_text(size = base_size - 2, color = "grey40")
    )
}


# -----------------------------------------------------------------------------
# type = "summary"
# -----------------------------------------------------------------------------

#' @keywords internal
.plot_summary <- function(sw, target_name, label_map, colors,
                          ribbon, ci_band, line_width, base_size) {

  bl         <- sw$baseline
  iv         <- sw$interventions[[target_name]]
  target_col <- unname(colors[target_name])
  steps      <- 0:sw$n_steps

  df_wide <- data.frame(
    step         = steps,
    baseline     = bl$summary,
    intervention = iv$summary,
    bl_q025      = bl$q025,
    bl_q975      = bl$q975,
    iv_q025      = iv$q025,
    iv_q975      = iv$q975,
    stringsAsFactors = FALSE
  )

  gg <- ggplot2::ggplot(df_wide, ggplot2::aes(x = .data$step))

  if (ribbon) {
    gg <- gg + ggplot2::geom_ribbon(
      ggplot2::aes(ymin = pmin(.data$baseline, .data$intervention),
                   ymax = pmax(.data$baseline, .data$intervention)),
      fill = target_col, alpha = 0.12
    )
  }

  if (ci_band) {
    gg <- gg +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = .data$bl_q025, ymax = .data$bl_q975),
        fill = "grey50", alpha = 0.10
      ) +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = .data$iv_q025, ymax = .data$iv_q975),
        fill = target_col, alpha = 0.08
      )
  }

  gg <- gg +
    ggplot2::geom_line(ggplot2::aes(y = .data$baseline),
                       color = "grey50", linewidth = line_width,
                       linetype = "dashed") +
    ggplot2::geom_line(ggplot2::aes(y = .data$intervention),
                       color = target_col, linewidth = line_width * 1.3) +
    ggplot2::annotate(
      "point", x = sw$n_steps,
      y = df_wide$intervention[nrow(df_wide)],
      color = target_col, size = 3
    ) +
    ggplot2::annotate(
      "point", x = sw$n_steps,
      y = df_wide$baseline[nrow(df_wide)],
      color = "grey50", size = 2.5, shape = 1
    ) +
    ggplot2::scale_x_continuous(
      name   = "Time Step",
      breaks = seq(0, sw$n_steps, by = max(1L, sw$n_steps %/% 7L)),
      expand = ggplot2::expansion(mult = c(0.02, 0.12))
    ) +
    ggplot2::scale_y_continuous(
      name   = "Mean Activation Rate",
      labels = function(x) paste0(round(x * 100, 1), "%")
    ) +
    .kim_base_theme(base_size)

  gg
}


# -----------------------------------------------------------------------------
# type = "node"
# -----------------------------------------------------------------------------

#' @keywords internal
.plot_node <- function(sw, target_name, label_map, colors,
                       line_width, base_size) {

  pd <- .kim_build_plot_data(sw, target_name)
  df <- pd$node

  df$node_label <- label_map[df$node]
  node_order    <- c(target_name, setdiff(sw$vars, target_name))
  df$node_label <- factor(df$node_label, levels = label_map[node_order])

  lvls      <- levels(df$node_label)
  lvls[1L]  <- paste0(lvls[1L], " [\u2605 target]")
  levels(df$node_label) <- lvls

  gg <- ggplot2::ggplot(df, ggplot2::aes(x = .data$step, y = .data$prob,
                                         color    = .data$condition,
                                         linetype = .data$condition)) +
    ggplot2::geom_line(linewidth = line_width) +
    ggplot2::facet_wrap(~ node_label, scales = "free_y", ncol = 3L) +
    ggplot2::scale_color_manual(
      values = c("Baseline" = "grey55", "Intervention" = "#E64B35")
    ) +
    ggplot2::scale_linetype_manual(
      values = c("Baseline" = "dashed", "Intervention" = "solid")
    ) +
    ggplot2::scale_x_continuous(
      name   = "Time Step",
      breaks = seq(0, sw$n_steps, by = max(1L, sw$n_steps %/% 5L))
    ) +
    ggplot2::scale_y_continuous(
      name   = "P(s = 1)",
      labels = function(x) paste0(round(x * 100, 1), "%")
    ) +
    .kim_base_theme(base_size) +
    ggplot2::theme(legend.position = "bottom")

  gg
}


# -----------------------------------------------------------------------------
# type = "diff"
# -----------------------------------------------------------------------------

#' @keywords internal
.plot_diff <- function(sw, target_name, label_map, colors,
                       highlight_target, line_width, base_size,
                       intervention_type = "alleviating") {

  pd <- .kim_build_plot_data(sw, target_name)
  df <- pd$diff
  df$node_label <- label_map[df$node]
  df$is_target  <- df$node == target_name

  # Y-axis label reflects actual intervention direction
  y_label <- if (intervention_type == "alleviating") {
    "\u0394P (Intervention \u2212 Baseline) | Negative = alleviating"
  } else {
    "\u0394P (Intervention \u2212 Baseline) | Positive = aggravating"
  }

  gg <- ggplot2::ggplot(df, ggplot2::aes(x = .data$step, y = .data$diff,
                                         color = .data$node_label,
                                         group = .data$node)) +
    ggplot2::geom_hline(yintercept = 0, color = "grey70",
                        linetype = "dotted", linewidth = 0.5)

  if (highlight_target) {
    df_other  <- df[!df$is_target, ]
    df_target <- df[df$is_target, ]
    gg <- gg +
      ggplot2::geom_line(data      = df_other,
                         linewidth = line_width * 0.6, alpha = 0.5) +
      ggplot2::geom_line(data      = df_target,
                         linewidth = line_width * 1.8)
  } else {
    gg <- gg +
      ggplot2::geom_line(linewidth = line_width)
  }

  gg <- gg +
    ggplot2::scale_color_manual(
      values = stats::setNames(unname(colors), label_map[names(colors)])
    ) +
    ggplot2::scale_x_continuous(
      name   = "Time Step",
      breaks = seq(0, sw$n_steps, by = max(1L, sw$n_steps %/% 7L))
    ) +
    ggplot2::scale_y_continuous(
      name   = y_label,
      labels = function(x) paste0(round(x * 100, 1), "%")
    ) +
    .kim_base_theme(base_size)

  gg
}


# -----------------------------------------------------------------------------
# type = "panel"
# -----------------------------------------------------------------------------

#' @keywords internal
.plot_panel <- function(sw, ranked, top_n, label_map, colors,
                        ribbon, line_width, base_size) {

  top_n       <- min(as.integer(top_n), length(ranked))
  top_targets <- ranked[seq_len(top_n)]

  missing <- setdiff(top_targets, names(sw$interventions))
  if (length(missing) > 0L)
    stop("The following targets are not in the simulation results: ",
         paste(missing, collapse = ", "), ".")

  steps <- 0:sw$n_steps

  df_list <- lapply(seq_along(top_targets), function(k) {
    tn <- top_targets[k]
    iv <- sw$interventions[[tn]]
    data.frame(
      step      = rep(steps, 2L),
      score     = c(sw$baseline$summary, iv$summary),
      condition = rep(c("Baseline", "Intervention"), each = sw$n_steps + 1L),
      target    = paste0("#", k, " ", label_map[tn]),
      stringsAsFactors = FALSE
    )
  })
  df        <- do.call(rbind, df_list)
  df$target <- factor(df$target, levels = unique(df$target))

  gg <- ggplot2::ggplot(df, ggplot2::aes(x = .data$step, y = .data$score))

  if (ribbon) {
    df_ribbon_list <- lapply(seq_along(top_targets), function(k) {
      tn <- top_targets[k]
      iv <- sw$interventions[[tn]]
      data.frame(
        step   = steps,
        ylo    = pmin(sw$baseline$summary, iv$summary),
        yhi    = pmax(sw$baseline$summary, iv$summary),
        target = paste0("#", k, " ", label_map[tn]),
        stringsAsFactors = FALSE
      )
    })
    df_ribbon        <- do.call(rbind, df_ribbon_list)
    df_ribbon$target <- factor(df_ribbon$target, levels = levels(df$target))

    gg <- gg +
      ggplot2::geom_ribbon(
        data        = df_ribbon,
        mapping     = ggplot2::aes(x    = .data$step,
                                   ymin = .data$ylo,
                                   ymax = .data$yhi),
        inherit.aes = FALSE,
        fill        = "#E64B35",
        alpha       = 0.10
      )
  }

  gg <- gg +
    ggplot2::geom_line(ggplot2::aes(color    = .data$condition,
                                    linetype = .data$condition),
                       linewidth = line_width) +
    ggplot2::facet_wrap(~ target, ncol = min(top_n, 3L), scales = "free_y") +
    ggplot2::scale_color_manual(
      values = c("Baseline" = "grey55", "Intervention" = "#E64B35")
    ) +
    ggplot2::scale_linetype_manual(
      values = c("Baseline" = "dashed", "Intervention" = "solid")
    ) +
    ggplot2::scale_x_continuous(
      name   = "Time Step",
      breaks = seq(0, sw$n_steps, by = max(1L, sw$n_steps %/% 5L))
    ) +
    ggplot2::scale_y_continuous(
      name   = "Mean Activation Rate",
      labels = function(x) paste0(round(x * 100, 1), "%")
    ) +
    .kim_base_theme(base_size)

  gg
}


# -----------------------------------------------------------------------------
# type = "spaghetti"
# -----------------------------------------------------------------------------

#' @keywords internal
.plot_spaghetti <- function(sw, target_name, label_map, colors,
                            line_width, base_size) {

  target_col <- unname(colors[target_name])
  steps      <- 0:sw$n_steps

  bl_traces <- sw$baseline$traces
  iv_traces <- sw$interventions[[target_name]]$traces

  if (is.null(bl_traces) || is.null(iv_traces))
    stop("Spaghetti plot requires individual trace data (n_traces > 0).")

  n_actual <- ncol(bl_traces)

  df_traces_list <- lapply(seq_len(n_actual), function(k) {
    data.frame(
      step      = rep(steps, 2L),
      score     = c(bl_traces[, k], iv_traces[, k]),
      condition = rep(c("Baseline", "Intervention"), each = length(steps)),
      id        = k,
      stringsAsFactors = FALSE
    )
  })
  df_traces <- do.call(rbind, df_traces_list)

  df_mean <- data.frame(
    step      = rep(steps, 2L),
    score     = c(sw$baseline$summary,
                  sw$interventions[[target_name]]$summary),
    condition = rep(c("Baseline", "Intervention"), each = length(steps)),
    stringsAsFactors = FALSE
  )

  gg <- ggplot2::ggplot() +
    ggplot2::geom_line(
      data    = df_traces,
      mapping = ggplot2::aes(x     = .data$step,
                             y     = .data$score,
                             group = interaction(.data$id, .data$condition),
                             color = .data$condition),
      alpha = 0.08, linewidth = 0.3
    ) +
    ggplot2::geom_line(
      data    = df_mean,
      mapping = ggplot2::aes(x        = .data$step,
                             y        = .data$score,
                             color    = .data$condition,
                             linetype = .data$condition),
      linewidth = line_width * 1.5
    ) +
    ggplot2::scale_color_manual(
      values = c("Baseline" = "grey45", "Intervention" = target_col)
    ) +
    ggplot2::scale_linetype_manual(
      values = c("Baseline" = "dashed", "Intervention" = "solid")
    ) +
    ggplot2::scale_x_continuous(
      name   = "Time Step",
      breaks = seq(0, sw$n_steps, by = max(1L, sw$n_steps %/% 7L))
    ) +
    ggplot2::scale_y_continuous(
      name   = "Mean Activation Rate",
      labels = function(x) paste0(round(x * 100, 1), "%")
    ) +
    .kim_base_theme(base_size)

  gg
}
