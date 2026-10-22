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
#' Records node activation rates, summary scores and per-trajectory scores at
#' every step. `window_mean` reproduces exactly what kim_trajectory() returns
#' for the same chain and window, without a second simulation -- this is what
#' keeps the stored trajectory and the reported effects on the same run.
#'
#' @param burn_in Leading steps excluded from `window_mean` (still recorded).
#' @keywords internal
.kim_trajectory_stepwise <- function(init_states, J, h, n_steps,
                                     target_node = NULL, delta = 0,
                                     U_list = NULL, n_traces = 0L,
                                     burn_in = 0L) {
  n_sim  <- nrow(init_states)
  n_vars <- ncol(init_states)
  vars   <- colnames(init_states)

  burn_in <- max(0L, as.integer(burn_in))
  if (burn_in >= n_steps) stop("burn_in must be < n_steps.")
  n_retained <- n_steps - burn_in

  h_eff <- h
  if (!is.null(target_node) && delta != 0) {
    if (is.matrix(h_eff)) h_eff[, target_node] <- h_eff[, target_node] + delta
    else h_eff[target_node] <- h_eff[target_node] + delta
  }

  # Pre-allocate
  node_probs     <- matrix(NA_real_, nrow = n_steps + 1L, ncol = n_vars)
  colnames(node_probs) <- vars
  summary_scores <- numeric(n_steps + 1L)
  q025           <- numeric(n_steps + 1L)
  q975           <- numeric(n_steps + 1L)
  cum_score      <- numeric(n_sim)        # retained-window accumulator

  # Individual traces (optional)
  trace_mat      <- NULL
  n_trace_actual <- 0L
  if (n_traces > 0L) {
    n_trace_actual <- min(as.integer(n_traces), n_sim)
    trace_mat <- matrix(NA_real_, nrow = n_steps + 1L, ncol = n_trace_actual)
  }

  # Step 0 (never enters the retained window: burn_in >= 0)
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
    if (step > burn_in) cum_score <- cum_score + row_scores
    if (n_trace_actual > 0L)
      trace_mat[step + 1L, ] <- row_scores[seq_len(n_trace_actual)]
  }

  list(node_probs  = node_probs,
       summary     = summary_scores,
       q025        = q025,
       q975        = q975,
       traces      = trace_mat,
       window_mean = cum_score / n_retained,
       burn_in     = burn_in,
       n_retained  = n_retained)
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

  list(summary = df_summary,
       node    = df_node,
       diff    = df_diff)
}


# -----------------------------------------------------------------------------
# Internal: EXACT per-step trajectory (deterministic distribution propagation)
# -----------------------------------------------------------------------------
.kim_exact_stepwise <- function(kim_fit, delta, target = NULL, n_steps = 20L) {
  h <- kim_fit$h; vars <- kim_fit$vars; p <- length(vars)
  ns  <- 2^p
  S   <- .kim_state_matrix(p)                  # ns x p (0/1)
  act <- rowSums(S)
  shifts <- kim_fit$person_shifts
  use_pop <- .kim_uses_population(kim_fit)
  blist <- if (use_pop && !is.null(shifts) && nrow(shifts) >= 1L) {
    lapply(seq_len(nrow(shifts)), function(k) h + shifts[k, ])
  } else list(h)

  node_probs <- matrix(0, n_steps + 1L, p, dimnames = list(NULL, vars))
  summary_scores <- numeric(n_steps + 1L)

  for (k in seq_along(blist)) {
    J  <- .kim_J_for(kim_fit, k)
    b  <- blist[[k]]
    b2 <- b; if (!is.null(target)) b2[target] <- b2[target] + delta
    # Use the same transition direction as runKIMNIRA.
    pr <- .kim_transition(J, b2)$probabilities
    # Baseline and intervention share the baseline stationary initial state.
    v0 <- .kim_stationary(.kim_transition(J, b)$P)$pi
    prod_dist <- function(pv) { w <- 1.0
      for (j in seq_len(p)) w <- c(w * (1 - pv[j]), w * pv[j]); w }
    # ns = 2^p. Below the threshold, materialising the full ns x ns kernel turns
    # each step into one matrix-vector product; above it that matrix costs more
    # memory than the loop costs time, so the steps are accumulated on the fly.
    # This is a speed/memory trade-off, not a cap on p -- both branches give the
    # same answer, and neither refuses to run.
    buildT <- ns <= 2048L
    if (buildT) {
      Tmat <- matrix(0, ns, ns)
      for (f in seq_len(ns)) Tmat[f, ] <- prod_dist(pr[f, ])
      step1 <- function(v) as.numeric(v %*% Tmat)
    } else {
      step1 <- function(v) { vn <- numeric(ns)
        for (f in seq_len(ns)) { if (v[f] == 0) next
          vn <- vn + v[f] * prod_dist(pr[f, ]) }
        vn }
    }
    v <- v0
    node_probs[1L, ]   <- node_probs[1L, ]   + as.numeric(v %*% S)
    summary_scores[1L] <- summary_scores[1L] + as.numeric(v %*% act) / p
    for (t in seq_len(n_steps)) {
      v <- step1(v)
      node_probs[t + 1L, ]   <- node_probs[t + 1L, ]   + as.numeric(v %*% S)
      summary_scores[t + 1L] <- summary_scores[t + 1L] + as.numeric(v %*% act) / p
    }
  }
  n <- length(blist)
  node_probs <- node_probs / n
  summary_scores <- summary_scores / n
  list(node_probs = node_probs, summary = summary_scores,
       q025 = summary_scores, q975 = summary_scores,   # no MC spread
       traces = NULL)
}

#' @keywords internal
.kim_run_stepwise_exact <- function(nira_result, target_nodes, n_steps = 20L) {
  kim_fit <- nira_result$kim_fit
  vars    <- kim_fit$vars
  delta   <- .kim_expand_delta(nira_result$delta, vars)
  bl <- .kim_exact_stepwise(kim_fit, 0, NULL, n_steps)
  iv <- lapply(target_nodes, function(j) .kim_exact_stepwise(kim_fit, delta[j], j, n_steps))
  names(iv) <- vars[target_nodes]
  list(baseline = bl, interventions = iv, vars = vars,
       n_steps = n_steps, delta = delta)
}


# -----------------------------------------------------------------------------
# Main plotting function
# -----------------------------------------------------------------------------

#' Plot KIM NodeIdentifyR Algorithm (NIRA) Longitudinal Intervention Trajectories
#'
#' @description
#' Draws the time course of a sustained perturbation.
#'
#' For a Monte Carlo result the path is read from `nira_result$trajectory`, which
#' [runKIMNIRA()] records from the very chains that produced `$effects`. Plotting
#' therefore describes the same run as the reported numbers instead of simulating
#' a second, different chain; the plotted per-step means averaged over the retained
#' window `(burn_in, n_steps]` reproduce `$effects`. A Monte Carlo result without a
#' stored trajectory is an error, not a silent re-simulation.
#'
#' An exact result stores no path: exact propagation is deterministic, so the path
#' is recomputed here and cannot drift away from `$effects`.
#'
#' Two quantities are shown and must not be conflated. `$effects` estimates
#' `Delta_j_stat`, the change in the stationary mean of g(S) under a sustained
#' perturbation of node j. This plot shows `Delta_j(t)`, the response path from the baseline stationary
#' state; averaging that path over the retained window `(burn_in, n_steps]`
#' reproduces `$effects`, but no single step is claimed to equal it.
#'
#' Responses start from the baseline stationary distribution under
#' `method = "exact"`, so the exact curve is the intervention response B_j(m):
#' the system is in its long-run baseline state when the sustained perturbation
#' begins. Only the exact curve may be read as "how quickly the effect emerges".
#' Under `method = "mc"` chains start from an independent Bernoulli draw at each
#' person's baseline, so the leading part of the curve is warm-up towards the
#' stationary law, not an intervention response; it is shown as a convergence
#' diagnostic. Monte Carlo ribbons describe trajectory dispersion, not sampling
#' uncertainty. One step is one fitted transition, i.e. one observation
#' interval of the data.
#'
#' @param nira_result Result returned by runKIMNIRA.
#' @param type Plot layout. `"summary"` shows baseline versus intervention mean
#'   trajectories; `"node"` shows per-node probabilities in free-scale facets;
#'   `"diff"` shows the per-node intervention-minus-baseline difference for one
#'   target; `"panel"` compares the top `top_n` targets on a shared scale;
#'   `"spaghetti"` draws individual Monte Carlo traces (MC only).
#' @param target Target node name or index; defaults to the first ranked target.
#' @param top_n Number of targets in a panel plot.
#' @param n_traces Number of individual Monte Carlo traces to retain. Capped at
#'   200, and at the number stored in `nira_result$trajectory`.
#' @param ribbon Shade the gap between baseline and intervention mean trajectories.
#' @param ci_band Show the pointwise central 95 percent simulation range, not a
#'   mean confidence interval or a person-bootstrap interval. MC only.
#' @param highlight_target Highlight the perturbed node.
#' @param node_labels Optional display labels for nodes.
#' @param colors Optional named vector of plotting colors.
#' @param base_size Base font size.
#' @param line_width Line width.
#' @param verbose Print progress messages.
#'
#' @return A ggplot object, returned invisibly.
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

  eng <- if (!is.null(nira_result$method)) nira_result$method else "mc"
  is_exact <- eng == "exact"
  # No silent switch of engine: a Monte Carlo result without its replicates is an
  # error, because an exact path would describe a different computation.
  if (!is_exact && is.null(nira_result$baseline_scores))
    stop("This Monte Carlo result carries no trajectory replicates. Re-run ",
         "runKIMNIRA(..., trajectory = TRUE), or use method = 'exact'.")

  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required: install.packages('ggplot2')")

  type <- match.arg(type)
  if (type == "spaghetti" && is_exact)
    stop("type = 'spaghetti' draws individual Monte-Carlo trajectories and ",
         "therefore requires runKIMNIRA(..., method = 'mc').")
  if (is_exact) {                        # no MC replicates: no ribbons/traces
    ribbon  <- FALSE; ci_band <- FALSE; n_traces <- 0L
  }

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
  # Keep the original palette; extend short palettes without missing colours.
  palette <- c(default_colors,
               grDevices::hcl.colors(max(n_vars, 30L), palette = "Dark 3"))
  supplied <- colors
  colors <- stats::setNames(rep(NA_character_, n_vars), vars)
  if (length(supplied)) {
    if (anyNA(supplied)) stop("colors cannot contain NA.")
    if (!is.null(names(supplied))) {
      if (anyDuplicated(names(supplied)) ||
          any(!names(supplied) %in% vars))
        stop("Named colors must use unique node names.")
      colors[names(supplied)] <- supplied
    } else {
      take <- seq_len(min(length(supplied), n_vars))
      colors[take] <- supplied[take]
    }
  }
  available <- unique(palette[!palette %in% colors[!is.na(colors)]])
  colors[is.na(colors)] <- available[seq_len(sum(is.na(colors)))]
  grDevices::col2rgb(colors) # Validate colour specifications.

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
  } else {
    target_indices <- match(target_name, vars)
  }

  if (verbose)
    message("Running stepwise simulation (", length(target_indices),
            " intervention condition(s))...")

  sw <- if (!is.null(nira_result$trajectory)) {
    # Same chains that produced the reported effects -- no second simulation.
    nira_result$trajectory
  } else if (is_exact) {
    # Exact paths are deterministic propagation, not simulation: recomputing them
    # returns the identical path and cannot drift away from the reported effects.
    .kim_run_stepwise_exact(nira_result,
                            target_nodes = target_indices,
                            n_steps      = 20L)
  } else {
    stop("This Monte Carlo result carries no stored trajectory, so plotting it ",
         "would mean simulating a second, different chain that need not reproduce ",
         "the reported effects. Re-run with runKIMNIRA(..., trajectory = TRUE) to ",
         "record the path of the same run.")
  }

  # Stored trajectories keep a fixed trace budget; honour the requested count.
  if (type == "spaghetti" && !is.null(sw$baseline$traces)) {
    nc <- min(n_traces, ncol(sw$baseline$traces))
    if (nc < n_traces)
      warning("Only ", nc, " stored trace(s) available; requested ", n_traces, ".")
    sw$baseline$traces <- sw$baseline$traces[, seq_len(nc), drop = FALSE]
    sw$interventions <- lapply(sw$interventions, function(z) {
      if (!is.null(z$traces)) z$traces <- z$traces[, seq_len(nc), drop = FALSE]
      z
    })
  }

  if (verbose) message("Generating plot from recorded trajectory...")

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

# One step is one fitted transition, i.e. one observation interval.
step_axis_label <- "Transitions since onset"

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
      name   = step_axis_label,
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

  # Align each initial baseline near the top, retaining absolute labels.
  # Equal percentage-point changes have equal heights among non-target nodes.
  # Reserve 90% of the range below baseline and 10% above it. Include all
  # upward responses/MC fluctuations; never clip data to force a flat baseline.
  baseline_rows <- df[df$condition == "Baseline" & df$step == 0, ]
  centres <- stats::setNames(baseline_rows$prob, baseline_rows$node)
  offsets <- df$prob - unname(centres[df$node])
  other_offsets <- offsets[df$node != target_name]
  span_for <- function(x) max(c(-x / .9, x / .1, 1e-4), na.rm = TRUE)
  shared_span <- span_for(other_offsets)
  target_span <- span_for(offsets[df$node == target_name])
  axis_bounds <- do.call(rbind, lapply(node_order, function(nm) {
    span <- if (nm == target_name) target_span else shared_span
    data.frame(step = 0, prob = centres[[nm]] + c(-.9, .1) * span,
               node_label = as.character(df$node_label[match(nm, df$node)]))
  }))
  axis_bounds$node_label <- factor(axis_bounds$node_label, levels = lvls)

  # 三条线：基线（灰虚线）、非靶点干预（红实线）、靶点干预（蓝实线）。
  df$line_type <- ifelse(df$condition == "Baseline", "Baseline",
                         ifelse(df$node == target_name, "Target", "Intervention"))
  df$line_type <- factor(df$line_type,
                         levels = c("Baseline", "Intervention", "Target"))

  gg <- ggplot2::ggplot(df, ggplot2::aes(x = .data$step, y = .data$prob,
                                         color    = .data$line_type,
                                         linetype = .data$line_type)) +
    ggplot2::geom_line(linewidth = line_width) +
    ggplot2::geom_blank(data = axis_bounds,
      mapping = ggplot2::aes(x = .data$step, y = .data$prob),
      inherit.aes = FALSE) +
    ggplot2::facet_wrap(~ node_label, scales = "free_y", ncol = 3L) +
    ggplot2::scale_color_manual(
      values = c("Baseline" = "grey55", "Intervention" = "#E64B35",
                 "Target" = "#3C5488"),
      name = NULL
    ) +
    ggplot2::scale_linetype_manual(
      values = c("Baseline" = "dashed", "Intervention" = "solid",
                 "Target" = "solid"),
      name = NULL
    ) +
    ggplot2::scale_x_continuous(
      name   = step_axis_label,
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
      name   = step_axis_label,
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
    ggplot2::facet_wrap(~ target, ncol = min(top_n, 3L), scales = "fixed") +
    ggplot2::scale_color_manual(
      values = c("Baseline" = "grey55", "Intervention" = "#E64B35")
    ) +
    ggplot2::scale_linetype_manual(
      values = c("Baseline" = "dashed", "Intervention" = "solid")
    ) +
    ggplot2::scale_x_continuous(
      name   = step_axis_label,
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
      name   = step_axis_label,
      breaks = seq(0, sw$n_steps, by = max(1L, sw$n_steps %/% 7L))
    ) +
    ggplot2::scale_y_continuous(
      name   = "Mean Activation Rate",
      labels = function(x) paste0(round(x * 100, 1), "%")
    ) +
    .kim_base_theme(base_size)

  gg
}
