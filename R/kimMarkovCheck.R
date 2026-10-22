# Contiguous blocks within person. Original positions are retained for embargoes.
.kim_order_folds <- function(person, n_folds) {
  out <- integer(length(person))
  for (pid in unique(person)) {
    ix <- which(person == pid)
    out[ix] <- pmin(n_folds, ceiling(seq_along(ix) * n_folds / length(ix)))
  }
  out
}

.kim_order_train <- function(pool, val, person, position, gap) {
  keep <- setdiff(pool, val)
  keep[!vapply(keep, function(i) any(person[val] == person[i] &
    abs(position[val] - position[i]) <= gap), logical(1))]
}

# Training and held-out designs: lagged states plus optional context columns.
.kim_order_design <- function(x, person, train, test, method, context) {
  list(train = cbind(x[train, , drop = FALSE], context[train, , drop = FALSE]),
       test = cbind(x[test, , drop = FALSE], context[test, , drop = FALSE]))
}

# Only the specific boundary message, confirmed by isSingular(), is benign
# for prediction. A singular fit can ALSO have optimizer/convergence failures.
.kim_mixed_diagnostics <- function(optinfo, singular, warnings = character(), messages = character()) {
  text <- function(x) unique(as.character(unlist(x, use.names = FALSE)))
  boundary <- function(x) singular & grepl("^boundary \\(singular\\) fit", x)
  conv <- text(optinfo$conv$lme4$messages)
  failures <- c(conv[!boundary(conv)], text(optinfo$warnings))
  for (part in c("opt", "lme4")) {
    codes <- if (part == "opt") optinfo$conv$opt else optinfo$conv$lme4$code
    codes <- unlist(codes, use.names = FALSE)
    if (length(codes) && any(is.na(codes) | codes != 0))
      failures <- c(failures, paste0(part, " convergence code: ", paste(codes, collapse = ", ")))
  }
  observed <- unique(c(warnings, messages))
  list(singular = singular, failures = unique(failures),
    notices = observed[!boundary(observed) & !observed %in% failures],
    boundary_messages = unique(c(conv[boundary(conv)], observed[boundary(observed)])))
}

.kim_valid_predictions <- function(prob, n) {
  if (!is.matrix(prob) || nrow(prob) != n || !ncol(prob) ||
      any(!is.finite(prob) | prob < 0 | prob > 1))
    stop("Invalid held-out predictions: require finite probabilities in [0, 1] for every validation row.")
  prob
}

.kim_order_predict <- function(x, y, person, train, test, method, context,
                               slope_structure) {
  if (length(train) < 5L || length(unique(y[train])) < 2L)
    stop("Insufficient training observations or no outcome variation.")
  dm <- .kim_order_design(x, person, train, test, method, context)
  colnames(dm$train) <- colnames(dm$test) <- paste0("x", seq_len(ncol(dm$train)))
  tr <- data.frame(y = y[train], dm$train, person = factor(person[train]))
  te <- data.frame(dm$test, person = factor(person[test], levels = levels(tr$person)))
  if (anyNA(te$person)) stop("No training observations for a validation person.")
  rhs <- colnames(dm$train)
  random <- if (method == "re") "(1 | person)" else paste0("(1 + ",
    paste(rhs[seq_len(ncol(x))], collapse = " + "),
    if (slope_structure == "full") " | person)" else " || person)")
  formula <- stats::as.formula(paste("y ~", paste(c(rhs, random), collapse = " + ")))
  fit_warnings <- fit_messages <- character()
  model <- tryCatch(withCallingHandlers(
    lme4::glmer(formula, data = tr, family = stats::binomial(),
      control = lme4::glmerControl(optimizer = "bobyqa")),
    warning = function(w) {
      fit_warnings <<- c(fit_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }, message = function(m) {
      fit_messages <<- c(fit_messages, trimws(conditionMessage(m)))
      invokeRestart("muffleMessage")
    }), error = function(e) stop(paste(unique(c(conditionMessage(e),
      fit_warnings, fit_messages)), collapse = " | "), call. = FALSE))
  dg <- .kim_mixed_diagnostics(model@optinfo, lme4::isSingular(model), fit_warnings, fit_messages)
  if (length(dg$failures)) stop(paste("Mixed model fit failed:",
    paste(unique(c(dg$failures, dg$notices, dg$boundary_messages)), collapse = " | ")), call. = FALSE)
  for (notice in dg$notices) warning(notice, call. = FALSE)
  pred <- .kim_valid_predictions(matrix(stats::predict(model, newdata = te,
    type = "response"), ncol = 1L), length(test))
  if (dg$singular) warning(structure(list(
    message = paste0("Singular random-effects fit; valid predictions retained. ",
      paste(dg$boundary_messages, collapse = " | ")),
    call = NULL), class = c("kim_singular_fit", "warning", "condition")))
  pred
}

.kim_order_loss <- function(y, prob) {
  prob <- pmin(pmax(prob, 1e-12), 1 - 1e-12)
  -(y * log(prob) + (1 - y) * log1p(-prob))
}

# Resample the same persons for every complete score column. A failed node
# never removes people/folds, and never suppresses another node's interval.
.kim_order_intervals <- function(gains, n_boot, conf_level) {
  interval <- matrix(NA_real_, ncol(gains), 2L)
  valid <- apply(gains, 2L, function(z) all(is.finite(z)))
  if (nrow(gains) >= 2L && any(valid)) {
    boot <- matrix(replicate(n_boot, colMeans(gains[
      sample.int(nrow(gains), nrow(gains), replace = TRUE), valid, drop = FALSE])),
      nrow = sum(valid), ncol = n_boot)
    interval[valid, ] <- t(apply(boot, 1L, stats::quantile,
      probs = c((1 - conf_level)/2, (1 + conf_level)/2), names = FALSE))
  }
  interval
}

#' Assess the incremental predictive value of a second lag
#'
#' Compares lag-1 with joint lag-1/lag-2 binary regressions on identical complete
#' three-occasion windows. This is a predictive sensitivity analysis, not a test
#' proving the Markov property, and its output cannot be used by runKIMNIRA.
#'
#' @param kim_fit A fitKIM object with saved settings and lags = 1.
#' @param data Original data, including original schedule indices and missing rows.
#' @param n_folds Number of outer within-person contiguous validation blocks.
#' @param cv_gap Embargo in original within-person row positions, at least two.
#' @param practical_gain Smallest meaningful absolute reduction in log loss
#'   (natural-log units per node/occasion). Default 0.005 is a user-adjustable
#'   working threshold, not a universal scientific cutoff.
#' @param n_boot Number of person-level resamples of held-out losses.
#' @param conf_level Confidence level for descriptive percentile intervals.
#' @param min_persons Minimum persons for an automated recommendation (default
#'   10, a reporting safeguard rather than a power guarantee).
#' @param context Optional column names for observed transition context, such as
#'   morning/evening. Response-occasion values enter both models as fixed
#'   effects. No context is inferred.
#' @param seed Seed for resampling; restored on exit.
#' @param verbose Print outer-fold progress.
#'
#' @details Time rules, min_obs, between_person, screen and slope_structure are
#'   inherited from kim_fit. Only fitted persons are analysed. Both models use the
#'   same eligible triplets and outer folds, and there is no penalty tuning: the
#'   same random-effect structure is used, with additional random lag-2 slopes
#'   for res. Outer validation estimates interpolation into held-out time blocks
#'   within known persons, not strictly prospective forecasting or
#'   generalisation to new persons.
#'
#'   Positive gain means lag-2 improves prediction. Persons receive equal weight
#'   after averaging their occasions; nodes receive equal weight in the overall
#'   score. Intervals resample persons' already-computed held-out losses without
#'   refitting. They are conditional, descriptive uncertainty summaries: shared
#'   training fits and model-selection uncertainty are not fully represented.
#'   Node intervals are exploratory and not multiplicity-adjusted.
#'
#'   lag2_improves means the overall interval lies above practical_gain;
#'   first_order_adequate means it lies below that threshold; otherwise the
#'   result is inconclusive. Inadequate person count, failed folds or warnings
#'   other than a confirmed singular boundary fit also force inconclusive.
#'   Singular random-effects fits with successful optimizer codes, no other
#'   reported convergence failures and valid predictions are retained and
#'   flagged. Singularity is not evidence that numerical reliability or the
#'   random-effects specification has been established. Recommendations and
#'   intervals with singular fits require particular caution.
#'   An adequate result concerns only this
#'   additive second-lag alternative, not all higher-order/nonlinear dependence.
#'   Excluding overnight transitions with two occasions per day leaves no
#'   eligible triplets; the function stops with a clear message.
#'
#' @return A kim_markov_check list with recommendation, overall and node score
#'   tables, person-level gains, out-of-fold predictions, pairing audit, outer
#'   split indices, fit issues and settings.
#'   `issue_details` distinguishes singular, warning and error events, preserving
#'   original messages. `fit_diagnostics` gives one row per node/model/outer fold,
#'   with status, training/validation counts and valid-prediction counts.
#'   Complete nodes receive intervals even when another node fails. Overall
#'   comparison requires all nodes and folds; no failed node/person is dropped.
#' @export
kimMarkovCheck <- function(kim_fit, data, n_folds = 5L,
    cv_gap = 2L, practical_gain = 0.005, n_boot = 1000L, conf_level = 0.95,
    min_persons = 10L, context = NULL, seed = 2026L, verbose = TRUE) {
  if (!inherits(kim_fit, "dynaNIRA_kim") || is.null(kim_fit$settings))
    stop("kim_fit must be a fitKIM object with saved settings.")
  spec <- kim_fit$settings
  if (!identical(as.integer(spec$lags), 1L)) stop("Fit a lag-1 main model before kimMarkovCheck.")
  for (value in list(n_folds, cv_gap, n_boot, min_persons))
    if (length(value) != 1L || !is.numeric(value) || !is.finite(value) || value != floor(value))
      stop("Fold, embargo, bootstrap and person counts must be finite integers.")
  if (n_folds < 3 || cv_gap < 2 || n_boot < 100 || min_persons < 2)
    stop("Require folds >= 3, cv_gap >= 2, n_boot >= 100 and min_persons >= 2.")
  if (length(practical_gain) != 1L || !is.finite(practical_gain) || practical_gain < 0 ||
      length(conf_level) != 1L || !is.finite(conf_level) || conf_level <= 0 || conf_level >= 1)
    stop("Invalid practical_gain or conf_level.")
  method <- spec$between_person
  if (!method %in% c("re", "res")) stop("Unsupported between_person method; refit with 're' or 'res'.")
  if (!requireNamespace("lme4", quietly = TRUE)) stop("lme4 is required.")
  args <- spec[intersect(names(spec), names(formals(.kim_time_pairs)))]
  args$data <- data; args$lags <- 1:2
  pairs <- do.call(.kim_time_pairs, args)
  eligible <- pairs$person %in% kim_fit$person_ids
  if (!sum(eligible)) stop("No eligible complete triplets. Three consecutive occasions are required; two beeps/day with overnight excluded cannot support this check.")
  person <- pairs$person[eligible]; position <- pairs$position[eligible]
  y <- pairs$Y[eligible, , drop = FALSE]
  x1 <- pairs$X[[1L]][eligible, , drop = FALSE]
  x2 <- cbind(x1, pairs$X[[2L]][eligible, , drop = FALSE])
  rawrows <- pairs$rows[eligible, , drop = FALSE]
  pairs$audit$used_in_diagnostic <- pairs$audit$response_row %in% rawrows[, 1L]
  if (any(table(person) < n_folds * 2L))
    stop("Too few eligible triplets per person for outer folds; reduce n_folds or revise the sampling rule.")
  ctx <- matrix(numeric(), nrow(y), 0L)
  if (length(context)) {
    if (!all(context %in% names(data))) stop("Missing context columns.")
    dd <- data[rawrows[, 1L], context, drop = FALSE]
    if (anyNA(dd)) stop("Missing context values on eligible response occasions.")
    names(dd) <- paste0("context", seq_along(context))
    ctx <- stats::model.matrix(~ ., dd)[, -1L, drop = FALSE]
    if (any(!is.finite(ctx))) stop("Non-finite context predictors.")
  }
  n <- nrow(y); p <- ncol(y); ids <- unique(person)
  folds <- .kim_order_folds(person, n_folds)
  pred1 <- pred2 <- matrix(NA_real_, n, p, dimnames = list(NULL, spec$vars))
  issues <- character(); issue_rows <- list(); fit_rows <- list()
  splits <- vector("list", n_folds)
  for (fold in seq_len(n_folds)) {
    if (verbose) message("Markov sensitivity: outer fold ", fold, "/", n_folds)
    val <- which(folds == fold)
    tr <- .kim_order_train(seq_len(n), val, person, position, cv_gap)
    splits[[fold]] <- list(train = tr, validation = val)
    for (i in seq_len(p)) for (order in 1:2) {
      x <- if (order == 1L) x1 else x2
      label <- paste0("fold ", fold, ", node ", spec$vars[i], ", order ", order)
      current_types <- current_messages <- character()
      record_issue <- function(type, message) {
        current_types <<- c(current_types, type)
        current_messages <<- c(current_messages, message)
        issues <<- c(issues, paste(label, message, sep = ": "))
        issue_rows[[length(issue_rows) + 1L]] <<- data.frame(fold = fold,
          node = spec$vars[i], order = order, type = type, message = message,
          stringsAsFactors = FALSE)
      }
      ans <- tryCatch(withCallingHandlers({
        .kim_order_predict(x, y[, i], person, tr, val, method, ctx,
          spec$slope_structure)
      }, warning = function(w) {
        record_issue(if (inherits(w, "kim_singular_fit")) "singular" else "warning", conditionMessage(w))
        invokeRestart("muffleWarning")
      }), error = function(e) {
        record_issue("error", conditionMessage(e))
        NULL
      })
      fit_rows[[length(fit_rows) + 1L]] <- data.frame(fold = fold, node = spec$vars[i],
        order = order, status = if (is.null(ans)) "failed" else
          if ("warning" %in% current_types) "warning" else
          if ("singular" %in% current_types) "singular" else "ok",
        singular = "singular" %in% current_types, n_train = length(tr),
        n_validation = length(val), n_predicted = if (is.null(ans)) 0L else nrow(ans),
        detail = paste(unique(current_messages), collapse = " | "), stringsAsFactors = FALSE)
      if (!is.null(ans)) {
        if (order == 1L) pred1[val, i] <- ans[, 1L] else pred2[val, i] <- ans[, 1L]
      }
    }
  }
  l1 <- .kim_order_loss(y, pred1); l2 <- .kim_order_loss(y, pred2)
  b1 <- (y - pred1)^2; b2 <- (y - pred2)^2
  average <- function(z) {
    out <- t(vapply(ids, function(pid) colMeans(z[person == pid, , drop = FALSE]), numeric(p)))
    dimnames(out) <- list(ids, spec$vars); out
  }
  a1 <- average(l1); a2 <- average(l2)
  gains <- cbind(overall = rowMeans(a1 - a2), a1 - a2)
  complete <- all(is.finite(c(pred1, pred2)))
  oldseed <- if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) get(".Random.seed", .GlobalEnv) else NULL
  on.exit(if (is.null(oldseed)) {
    if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
  } else assign(".Random.seed", oldseed, .GlobalEnv), add = TRUE)
  set.seed(seed)
  interval <- .kim_order_intervals(gains, n_boot, conf_level)
  issue_details <- if (length(issue_rows)) do.call(rbind, issue_rows) else
    data.frame(fold = integer(), node = character(), order = integer(),
      type = character(), message = character(), stringsAsFactors = FALSE)
  fit_diagnostics <- do.call(rbind, fit_rows)
  score <- data.frame(node = c("overall", spec$vars),
    log_loss_lag1 = c(mean(a1), colMeans(a1)), log_loss_lag12 = c(mean(a2), colMeans(a2)),
    gain = colMeans(gains), lower = interval[, 1L], upper = interval[, 2L],
    brier_lag1 = c(mean(average(b1)), colMeans(average(b1))),
    brier_lag12 = c(mean(average(b2)), colMeans(average(b2))), row.names = NULL)
  score$complete <- c(complete, apply(is.finite(pred1) & is.finite(pred2), 2L, all))
  count_status <- function(status) vapply(spec$vars, function(v)
    sum(fit_diagnostics$node == v & fit_diagnostics$status %in% status), integer(1))
  singular_counts <- vapply(spec$vars, function(v)
    sum(fit_diagnostics$node == v & fit_diagnostics$singular), integer(1))
  failed_counts <- count_status("failed")
  score$singular_fits <- c(sum(singular_counts), singular_counts)
  score$failed_fits <- c(sum(failed_counts), failed_counts)
  status <- "inconclusive"
  if (complete && !any(issue_details$type != "singular") && length(ids) >= min_persons) {
    if (score$lower[1L] > practical_gain) status <- "lag2_improves"
    else if (score$upper[1L] < practical_gain) status <- "first_order_adequate"
  }
  explanation <- switch(status,
    lag2_improves = "The joint lag-1/lag-2 model shows additional predictive value beyond the chosen practical threshold; reconsider the first-order approximation.",
    first_order_adequate = "The first-order model is provisionally adequate relative to this additive lag-2 alternative and the chosen threshold. This does not establish the Markov property.",
    inconclusive = "Evidence is insufficient for a recommendation: inspect the interval, person count and fit issues. This is not evidence that the first-order assumption holds.")
  if (any(fit_diagnostics$singular)) explanation <- paste(explanation,
    "Caution: singular random-effects fits were retained for prediction; review $fit_diagnostics and random-effects adequacy.")
  result <- list(recommendation = status, explanation = explanation,
    overall = score[1L, , drop = FALSE], nodes = score[-1L, , drop = FALSE],
    person_gains = gains, predictions = list(lag1 = pred1, lag12 = pred2, observed = y),
    windows = data.frame(person = person, position = position, response_row = rawrows[, 1L],
      lag1_row = rawrows[, 2L], lag2_row = rawrows[, 3L], fold = folds),
    pair_audit = pairs$audit, pair_summary = pairs$summary,
    pair_time_unit = pairs$time_unit,
    splits = splits, issues = unique(issues),
    issue_details = issue_details, fit_diagnostics = fit_diagnostics,
    n_persons = length(ids), n_triplets = n, complete = complete,
    settings = list(main_model = spec, n_folds = n_folds,
      cv_gap = cv_gap, practical_gain = practical_gain,
      n_boot = n_boot, conf_level = conf_level, min_persons = min_persons,
      context = context, seed = seed),
    uncertainty = "Descriptive person bootstrap of held-out losses, conditional on fitted CV models; no refitting or multiplicity adjustment.")
  class(result) <- "kim_markov_check"
  result
}

#' @export
print.kim_markov_check <- function(x, ...) {
  cat("KIM lag-order predictive sensitivity\n", x$explanation, "\n", sep = "")
  cat("Persons:", x$n_persons, " Complete triplets:", x$n_triplets, "\n")
  print(x$overall, row.names = FALSE)
  cat("Practical gain threshold:", x$settings$practical_gain, "\n", x$uncertainty, "\n")
  if (!is.null(x$fit_diagnostics)) {
    cat("Singular fits with retained predictions:", sum(x$fit_diagnostics$singular),
      " Failed fits:", sum(x$fit_diagnostics$status == "failed"), "\n")
    if (length(x$issues)) cat("See $issue_details and $fit_diagnostics for original messages and fold status.\n")
  } else if (length(x$issues)) cat("Fit issues:", length(x$issues), "(see $issues)\n")
  invisible(x)
}
