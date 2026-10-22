# =============================================================================
# fitKIM.R  ----  KIM estimator (Kinetic Ising Model for dynamic networks)
#
# Model:
#   For symptom i of person k at measurement t:
#     P(s_ik(t+1) = 1 | s_k(t)) = sigma( h_i + a_ki + sum_j J_ij * s_jk(t) )
#   where a_ki is the person baseline (deviation from the typical person).
#
# `between_person` chooses how a_ki is handled. Only two options are offered:
#   "re"  : random intercept a_ki ~ N(0, tau_i^2) (lme4::glmer). Default.
#           Partial pooling; the natural multilevel form.
#   "res" : random intercept AND random slopes, i.e. person k also gets their
#           own transition matrix J_k = J + u_k. Answers "what if people differ
#           in the temporal relations themselves, not only in their baseline
#           levels?". slope_structure = "diag" (independent slope variances) or
#           "full" (slope covariance matrix); the default "auto" picks "full"
#           below 6 nodes and "diag" from 6 nodes on, mirroring mlVAR's default
#           rule, because correlated slope blocks are routinely singular.
#
# Edge screening (`screen`)
#   A cell of J is kept only if it is distinguishable from zero; the kept cells
#   are then refitted WITHOUT a penalty so the reported weights carry no
#   shrinkage. This is the two-stage ("relaxed lasso") structure, and the
#   screening rule was chosen from a simulation study:
#     "fdr"          (default) Benjamini-Hochberg over every candidate cell of
#                    the whole p x p network, at level `fdr_q`. It matched
#                    Bonferroni wherever edges were detectable and kept usable
#                    power when they were not (F1 0.79 vs 0.54 at |J| in
#                    0.25-0.45).
#     "significance" Bonferroni at level `alpha`; higher specificity, lower
#                    power on weak edges.
#     "none"         no screening: every cell is a candidate (the old dense
#                    behaviour, and what `screen = "none"` reproduces).
#   Note that an information-criterion penalty path (rpql/glmnet style) is NOT
#   used: on the same simulations every such criterion put its minimum at the
#   least-penalised end of the path and gave specificity near 0.4, because the
#   likelihood gap between the sparse and dense model is far larger than
#   log(n) per parameter.
#
# What the fit carries, in addition to the population (J, h):
#   support       : p x p logical, the cells that survived screening
#   stage1_beta   : the unpenalised coefficients that were screened
#   stage1_se     : their standard errors
#   screen_threshold : the |z| cut-off (Bonferroni) or the BH level (FDR)
#   person_shifts : n_persons x p matrix with a_ki. Zero when screen leaves no
#                   cell for a node.
#   person_sd     : length-p vector of per-node SD of a_ki.
# The (J, h + person_shifts) combination defines the population generative
# model used for *population-average* (marginal) intervention simulation.
# For "res", person_J replaces J per person in the same way, and the mean of
# person_J over persons is exactly J (the slopes are centred like the intercepts).
# =============================================================================


# -----------------------------------------------------------------------------
# Internal helper: between-subject (cross-sectional, EBICglasso) network
# -----------------------------------------------------------------------------
.estimate_between_kim <- function(data, vars, id, gamma = 0.5, verbose = TRUE) {

  .null_result <- function(dropped = character(0))
    list(network = NULL, dropped_nodes = dropped)

  if (!requireNamespace("qgraph", quietly = TRUE)) {
    warning("Package 'qgraph' is required to estimate the between-subject network: ",
            "install.packages('qgraph'). Returning NULL.")
    return(.null_result())
  }

  person_ids <- unique(data[[id]])
  person_means <- do.call(rbind, lapply(person_ids, function(pid) {
    pdata <- data[data[[id]] == pid, vars, drop = FALSE]
    colMeans(pdata, na.rm = TRUE)
  }))
  rownames(person_means) <- person_ids
  colnames(person_means) <- vars
  n_persons <- nrow(person_means); n_vars <- length(vars)

  if (n_persons < n_vars + 1L) {
    warning("Number of persons (", n_persons, ") is less than variables + 1 (",
            n_vars + 1L, "). Returning NULL.")
    return(.null_result())
  }

  col_var      <- apply(person_means, 2, var, na.rm = TRUE)
  dropped_vars <- vars[col_var < 1e-10]
  valid_vars   <- vars[col_var >= 1e-10]
  if (length(dropped_vars) > 0L) {
    if (verbose) message("  Between-subject network: ",
                         length(dropped_vars), " node(s) dropped (zero between-person variance).")
    if (length(valid_vars) < 3L) return(.null_result(dropped_vars))
    person_means <- person_means[, valid_vars, drop = FALSE]
  }
  n_valid <- length(valid_vars)

  cor_mat <- cor(person_means, use = "pairwise.complete.obs")
  if (any(is.na(cor_mat))) { cor_mat[is.na(cor_mat)] <- 0
    warning("Correlation matrix contained NAs; set to 0.") }

  eig_vals <- eigen(cor_mat, only.values = TRUE)$values
  if (any(eig_vals <= 0)) {
    if (verbose) message("  Correlation matrix not positive definite; nearPD.")
    if (!requireNamespace("Matrix", quietly = TRUE)) {
      cor_mat <- cor_mat + diag(n_valid) * (abs(min(eig_vals)) + 0.001)
      D_inv   <- diag(1 / sqrt(diag(cor_mat))); cor_mat <- D_inv %*% cor_mat %*% D_inv
    } else {
      cor_mat <- as.matrix(Matrix::nearPD(cor_mat, corr = TRUE, keepDiag = TRUE)$mat)
    }
    dimnames(cor_mat) <- list(valid_vars, valid_vars)
  }

  between_net_small <- tryCatch(
    qgraph::EBICglasso(S = cor_mat, n = n_persons, gamma = gamma,
                       returnAllResults = FALSE),
    error = function(e) { warning("EBICglasso failed: ", e$message); NULL })
  if (is.null(between_net_small)) return(.null_result(dropped_vars))

  if (length(dropped_vars) > 0L) {
    between_net <- matrix(NA_real_, n_vars, n_vars, dimnames = list(vars, vars))
    between_net[valid_vars, valid_vars] <- between_net_small
  } else {
    between_net <- between_net_small
    dimnames(between_net) <- list(vars, vars)
  }
  list(network = between_net, dropped_nodes = dropped_vars)
}


# -----------------------------------------------------------------------------
# Internal helper: fit ONE node of the multilevel transition model
# -----------------------------------------------------------------------------
# `keep` is the logical p-vector saying which lagged predictors enter this node's
# regression. Returns NULL when the fit fails so the caller can fall back.
.kim_fit_node <- function(y, Xlag, person, pids, keep, between_person,
                          slope_structure, control) {
  if (!any(keep)) return(NULL)
  Xs  <- Xlag[, keep, drop = FALSE]
  nms <- sub("^lag_", "", colnames(Xs))
  gd  <- data.frame(y = y, Xs, person = factor(person))
  rhs <- colnames(Xs)
  form <- if (between_person == "re") {
    stats::as.formula(sprintf("y ~ %s + (1 | person)", paste(rhs, collapse = " + ")))
  } else {
    op <- if (slope_structure == "full") "|" else "||"
    stats::as.formula(sprintf("y ~ %s + (1 + %s %s person)",
                              paste(rhs, collapse = " + "),
                              paste(rhs, collapse = " + "), op))
  }
  fit <- tryCatch(lme4::glmer(form, data = gd, family = stats::binomial,
                              control = control),
                  error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  co <- summary(fit)$coefficients
  re <- lme4::ranef(fit)$person
  re <- re[match(pids, rownames(re)), , drop = FALSE]
  list(fit = fit, co = co, re = re, nms = nms,
       singular = if (between_person == "res") lme4::isSingular(fit) else FALSE)
}


# -----------------------------------------------------------------------------
# Main estimation function
# -----------------------------------------------------------------------------

#' Fit a Kinetic Ising Model to multi-person binary EMA data
#'
#' Estimates the temporal interaction matrix J and intercept vector h with
#' conditional transition likelihood (node-wise multilevel logistic regressions).
#' `between_person` selects how the between-person heterogeneity a_ki is handled:
#' `"re"` (random intercepts; default) or `"res"` (random intercepts AND random
#' slopes, so each person also gets their own transition matrix).
#'
#' Edges are screened before the reported weights are computed: a cell of J is
#' kept only if it is distinguishable from zero, and the kept cells are then
#' refitted without a penalty. See `screen`.
#'
#' @param data Long-format data frame, one row per EMA measurement.
#' @param vars Character vector of binary (0/1) symptom variable names.
#' @param id Character, subject ID column. Default `"id"`.
#' @param time Original global measurement-index column, OR an actual Date/POSIX
#'   time column for advanced elapsed-time filtering. Integer numeric indices
#'   automatically detect missing scheduled occasions (1, 2, 4 is not continuous).
#'   Omit when using dayvar + beepvar. If neither is specified, a column named
#'   `time` is used for backward compatibility. Never renumber after missingness.
#' @param lags Integer lag distance. Default `1L`; a value of two uses only
#'   the second lag, not jointly the first and second lags.
#' @param between_person One of `"re"` (default) or `"res"`.
#'   `"re"` gives each person a random intercept per node (partial pooling).
#'   `"res"` adds random slopes, so each person also has their own transition
#'   matrix; see `slope_structure`.
#' @param slope_structure Only used when `between_person = "res"`. `"auto"`
#'   (default) mirrors the default rule of `mlVAR` for multilevel temporal
#'   networks: `"full"` (the complete slope covariance matrix) for fewer than six
#'   nodes and `"diag"` (independent slope variances) from six nodes on.
#'   Correlated slope covariances add p(p+1)/2 parameters per node on top of
#'   the p variances and are routinely singular. Passing `"diag"` or `"full"`
#'   explicitly always wins over the rule.
#' @param slope_warn_singular Warn when a `"res"` fit is singular.
#' @param screen How edges are screened before the reported weights:
#'   `"fdr"` (default) keeps cells whose unpenalised coefficient is significant
#'   under Benjamini-Hochberg control at level `fdr_q`, applied over all `p^2`
#'   candidate cells of the network; `"significance"` uses a Bonferroni-corrected
#'   level `alpha` instead; `"none"` screens nothing, so every cell is a
#'   candidate and J is dense (the behaviour of earlier versions).
#' @param fdr_q Benjamini-Hochberg level for `screen = "fdr"`. Default `0.10`.
#' @param alpha Family-wise level for `screen = "significance"`. Default `0.05`.
#' @param min_obs Minimum lag pairs retained per subject. Default `10L`.
#' @param max_time_gap Maximum elapsed time between each adjacent recorded
#'   measurement in a lag window (NULL = no upper bound). For lags greater
#'   than one, checked per step, not on the total window duration.
#' @param min_time_gap Minimum elapsed time per adjacent measurement (NULL = no
#'   lower bound). Neither bound establishes a first-order Markov process.
#' @param gamma EBIC hyperparameter for the between-subject network.
#' @param seed Random seed (kept for reproducibility of any downstream sampling).
#' @param verbose Logical.
#' @param estimate_between Whether to additionally estimate a between-person network.
#' @param consec Compatibility/advanced option: a separate global integer
#'   schedule-index column, for example alongside a real timestamp.
#' @param dayvar,beepvar Optional day and within-day scheduled measurement columns.
#'   Together these replace time entirely: records are sorted by person, day
#'   and beep. Days must be integer study-day indices or Date; beeps start at
#'   one and preserve missing slots. Nonconsecutive within-day beeps are excluded.
#' @param allow_overnight Keep cross-day transitions when dayvar is supplied?
#'   Default FALSE.
#' @param beeps_per_day Planned number of measurements per day, only needed for
#'   overnight day/beep pairing without a global index. Never inferred from
#'   observed data.
#' @param time_unit Units for Date/POSIX time-gap bounds: auto, secs, mins, hours,
#'   days or weeks.
#'
#' @return A `"dynaNIRA_kim"` list with `J`, `h`, `support` (p x p logical, the
#'   cells that survived screening), `stage1_beta`/`stage1_se` (the unpenalised
#'   coefficients that were screened and their standard errors),
#'   `screen_threshold`, `person_shifts`, `person_sd`, `between_network`,
#'   spectral diagnostics, and call metadata.
#'
#'   For `between_person = "res"` it additionally carries `person_J`
#'   (`n_persons x p x p`, person k's transition matrix, whose mean over persons
#'   is `J`), `slope_sd_est` (`p x p`, model-estimated random-slope SD: the
#'   quantity to report for "do slopes vary across persons?"), `slope_sd`
#'   (`p x p`, SD of the fitted deviations behind `person_J`, which is smaller
#'   because conditional modes are shrunk), `slope_varcorr` (per-node `VarCorr`
#'   table, including slope correlations when `slope_structure = "full"`), and
#'   `slope_singular` (per-node boundary-fit flag).
#'
#'   `person_constant` (`n_persons x p`) marks person-node cells whose outcome
#'   never varies. Under `"re"`/`"res"` such a person's random intercept is
#'   determined by shrinkage towards the population distribution rather than by
#'   their own data (the conditional mode is finite but prior-driven).
#'
#'   `spectral_radius_star` is the weighted-Jacobian spectral radius at the
#'   mean-field fixed point of the typical person; `spectral_radius_star_max` is
#'   its maximum over persons. Both are local, linearised heuristics; the exact
#'   relaxation time of the finite chain is given by [kimMixing()].
#' @export
fitKIM <- function(data,
                   vars,
                   id               = "id",
                   time             = NULL,
                   lags             = 1L,
                   between_person   = c("re", "res"),
                   slope_structure  = c("auto", "diag", "full"),
                   slope_warn_singular = TRUE,
                   screen           = c("fdr", "significance", "none"),
                   fdr_q            = 0.10,
                   alpha            = 0.05,
                   min_obs          = 10L,
                   max_time_gap     = NULL,
                   gamma            = 0.5,
                   seed             = 42L,
                   verbose          = TRUE,
                   estimate_between = TRUE,
                   min_time_gap     = NULL,
                   consec           = NULL,
                   dayvar           = NULL,
                   beepvar          = NULL,
                   allow_overnight  = FALSE,
                   time_unit        = "auto",
                   beeps_per_day    = NULL) {

  between_person  <- match.arg(between_person)
  slope_structure <- match.arg(slope_structure)
  screen          <- match.arg(screen)
  if (length(lags) != 1L || !is.finite(lags) || lags < 1 ||
      lags != as.integer(lags) || min_obs < 2)
    stop("Invalid lags or min_obs.")
  if (fdr_q <= 0 || fdr_q >= 1) stop("`fdr_q` must be in (0, 1).")
  if (alpha <= 0 || alpha >= 1) stop("`alpha` must be in (0, 1).")

  if (!requireNamespace("lme4", quietly = TRUE))
    stop("between_person='", between_person,
         "' requires package 'lme4': install.packages('lme4')")

  if (!is.data.frame(data)) stop("`data` must be a data.frame.")
  miss <- setdiff(c(vars, id, time), names(data))
  if (length(miss) > 0) stop("Columns missing from data: ", paste(miss, collapse = ", "))
  non_binary <- vars[!vapply(vars, function(v)
    all(stats::na.omit(unique(data[[v]])) %in% c(0, 1)), logical(1))]
  if (length(non_binary) > 0) stop("Non-binary variables: ", paste(non_binary, collapse = ", "))
  n_vars <- length(vars)
  if (n_vars < 2 || anyDuplicated(vars)) stop("Use at least two distinct node variables.")

  # Random-slope structure, resolved from the node count. This mirrors the
  # default rule in mlVAR::mlVAR (temporal = "default"): correlated slopes for
  # fewer than six nodes, orthogonal slopes otherwise. The reason is
  # identification, not taste -- a correlated block adds p(p+1)/2 covariance
  # parameters per node, and mlVAR's own documentation marks correlated random
  # effects "no longer practical".
  slope_requested <- slope_structure
  if (slope_structure == "auto")
    slope_structure <- if (between_person == "res" && n_vars < 6L) "full" else "diag"
  if (slope_requested == "full" && between_person != "res")
    warning("`slope_structure` is only used when between_person = 'res'; ignored.",
            call. = FALSE)
  if (slope_requested == "auto" && between_person == "res" && verbose)
    message("  slope_structure = 'auto' resolved to '", slope_structure,
            "' for ", n_vars, " nodes", if (n_vars >= 6L) " (>= 6)." else ".")
  if (anyNA(data[c(id, time)])) stop("Missing person IDs or time indices.")
  if (!is.null(time) && anyDuplicated(data[c(id, time)]))
    stop("Repeated person-time indices: use a unique ordered measurement index.")

  # -- keep subjects with enough observations ---------------------------------
  obs_counts <- tapply(seq_len(nrow(data)), data[[id]], length)
  excluded   <- names(obs_counts)[obs_counts < (min_obs + lags)]
  if (length(excluded) > 0 && verbose)
    message("Excluding ", length(excluded), " subject(s) with fewer than ",
            min_obs + lags, " observations.")
  data_clean <- data[!data[[id]] %in% excluded, ]
  n_persons  <- length(unique(data_clean[[id]]))
  if (n_persons == 0L) stop("No subjects with sufficient observations.")

  # ==========================================================================
  # Section 1: build lag pairs (s(t) -> s(t+1)) within person
  # ==========================================================================
  pairs <- .kim_time_pairs(data, vars, id, time, lags, min_obs,
    max_time_gap, min_time_gap, consec, dayvar, beepvar, allow_overnight, time_unit,
    beeps_per_day)
  X_lag <- pairs$X[[1L]]; Y_t1 <- pairs$Y; person_id_pool <- pairs$person
  storage.mode(X_lag) <- "double"; storage.mode(Y_t1) <- "double"
  total_na_removed <- sum(pairs$audit$reason == "missing")
  total_gap_removed <- sum(pairs$audit$reason == "time_gap")
  excluded_post_clean <- setdiff(pairs$excluded, excluded)

  if (verbose) {
    message("  Valid lag pairs: ", nrow(X_lag))
    if (total_na_removed > 0L) message("  Lag pairs removed (missing): ", total_na_removed)
    if (!is.null(max_time_gap) && total_gap_removed > 0L)
      message("  Lag pairs removed (time gaps): ", total_gap_removed)
  }
  if (is.null(X_lag) || nrow(X_lag) < n_vars * 2L) stop("Insufficient lag pairs.")

  transition_types <- NULL
  if (!is.null(dayvar) && !is.null(beepvar)) {
    kept <- pairs$rows
    from_b <- data[[beepvar]][kept[, ncol(kept)]]; to_b <- data[[beepvar]][kept[, 1L]]
    cross <- data[[dayvar]][kept[, ncol(kept)]] != data[[dayvar]][kept[, 1L]]
    lab <- paste0(ifelse(cross, "overnight ", ""), "beep ", from_b, " -> ", to_b)
    tt <- table(lab)
    transition_types <- data.frame(transition = names(tt), n = as.integer(tt),
                                   stringsAsFactors = FALSE)
    if (verbose) {
      message("  Retained transition types: ",
              paste(sprintf("%s (n=%d)", transition_types$transition, transition_types$n),
                    collapse = "; "))
      if (nrow(transition_types) == 1L)
        message("  All retained transitions are of one type, so one model step is ",
                "exactly that interval; simulated multi-step paths repeat it.")
    }
  }

  colnames(X_lag) <- paste0("lag_", vars)
  pids_clean <- sort(unique(person_id_pool))
  n_p <- length(pids_clean)
  ctrl <- lme4::glmerControl(optimizer = "bobyqa")
  eps <- 0.01

  # ==========================================================================
  # Section 2: pass 1 -- unpenalised multilevel fit of every node
  # ==========================================================================
  # This fit serves two purposes: its coefficients are what the screen tests,
  # and when nothing is screened out it is also the reported estimate.
  beta1 <- matrix(NA_real_, n_vars, n_vars, dimnames = list(vars, vars))
  se1   <- beta1
  pass1 <- vector("list", n_vars)
  node_fit_failed <- character(0); node_zero_var <- character(0)
  person_constant <- matrix(FALSE, n_p, n_vars, dimnames = list(pids_clean, vars))
  all_keep <- rep(TRUE, n_vars)

  for (i in seq_len(n_vars)) {
    y_i <- Y_t1[, i]
    pc <- tapply(y_i, person_id_pool, function(z) length(unique(z)) == 1L)
    person_constant[, i] <- as.logical(pc[as.character(pids_clean)])
    if (var(y_i) < 1e-10) { node_zero_var <- c(node_zero_var, vars[i]); next }
    if (verbose) message("  [fit] node ", i, "/", n_vars, ": ", vars[i])
    f1 <- .kim_fit_node(y_i, X_lag, person_id_pool, pids_clean, all_keep,
                        between_person, slope_structure, ctrl)
    if (is.null(f1)) { node_fit_failed <- c(node_fit_failed, vars[i]); next }
    pass1[[i]] <- f1
    co <- f1$co
    beta1[i, f1$nms] <- co[paste0("lag_", f1$nms), 1L]
    se1[i, f1$nms]   <- co[paste0("lag_", f1$nms), 2L]
  }

  # ==========================================================================
  # Section 3: screen
  # ==========================================================================
  support <- matrix(TRUE, n_vars, n_vars, dimnames = list(vars, vars))
  screen_threshold <- NA_real_
  if (screen == "none") {
    if (verbose) message("  [screen] none: all ", n_vars^2, " cells are candidates.")
  } else if (screen == "significance") {
    z <- stats::qnorm(1 - alpha / (2 * n_vars^2))
    zz <- abs(beta1) / se1
    zz[!is.finite(zz)] <- 0
    support[] <- zz > z
    screen_threshold <- z
    if (verbose) message("  [screen] Bonferroni |z| > ", round(z, 3),
                         " (alpha = ", alpha, " over ", n_vars^2, " cells): ",
                         sum(support), " selected.")
  } else {   # "fdr"
    pv <- 2 * stats::pnorm(-abs(beta1 / se1))
    pv[!is.finite(pv)] <- 1
    ord <- order(pv); m <- length(pv)
    pass <- pv[ord] <= (seq_len(m) / m) * fdr_q
    support[] <- FALSE
    if (any(pass)) support[ord[seq_len(max(which(pass)))]] <- TRUE
    screen_threshold <- fdr_q
    if (verbose) message("  [screen] BH-FDR q = ", fdr_q, " over ", m,
                         " cells: ", sum(support), " selected.")
  }

  # ==========================================================================
  # Section 4: pass 2 -- refit the surviving cells WITHOUT a penalty
  # ==========================================================================
  J_mat  <- matrix(0, n_vars, n_vars, dimnames = list(vars, vars))
  h_vec  <- setNames(numeric(n_vars), vars)
  person_shifts <- matrix(0, n_p, n_vars, dimnames = list(pids_clean, vars))
  person_sd <- setNames(numeric(n_vars), vars)
  person_J <- if (between_person == "res")
    array(0, c(n_p, n_vars, n_vars), dimnames = list(pids_clean, vars, vars)) else NULL
  slope_sd  <- if (between_person == "res")
    matrix(NA_real_, n_vars, n_vars, dimnames = list(vars, vars)) else NULL
  slope_sd_est <- if (between_person == "res")
    matrix(NA_real_, n_vars, n_vars, dimnames = list(vars, vars)) else NULL
  slope_varcorr  <- setNames(vector("list", n_vars), vars)
  slope_singular <- setNames(logical(n_vars), vars)

  for (i in seq_len(n_vars)) {
    y_i <- Y_t1[, i]
    if (var(y_i) < 1e-10) {
      p_m <- mean(y_i); p_m <- max(min(p_m, 1 - eps), eps)
      h_vec[i] <- log(p_m / (1 - p_m)); next
    }
    keep <- support[i, ]
    if (!any(keep)) {                     # screened out entirely
      f1 <- pass1[[i]]
      if (is.null(f1)) next
      u1 <- f1$re[, 1L]
      h_vec[i] <- f1$co["(Intercept)", 1L] + mean(u1)
      person_shifts[, i] <- u1 - mean(u1)
      person_sd[i] <- stats::sd(u1)
      next
    }
    # reuse the pass-1 fit when the screen kept everything
    fit_obj <- if (all(keep) && !is.null(pass1[[i]])) pass1[[i]] else {
      if (verbose) message("  [refit] node ", i, "/", n_vars, ": ", vars[i])
      .kim_fit_node(y_i, X_lag, person_id_pool, pids_clean, keep,
                    between_person, slope_structure, ctrl)
    }
    if (is.null(fit_obj)) {
      if (!vars[i] %in% node_fit_failed) node_fit_failed <- c(node_fit_failed, vars[i])
      next
    }
    fe  <- lme4::fixef(fit_obj$fit)
    re  <- fit_obj$re
    nms <- fit_obj$nms
    h_vec[i] <- fe["(Intercept)"]
    J_vec <- stats::setNames(as.numeric(fe[paste0("lag_", nms)]), nms)
    J_vec[!is.finite(J_vec)] <- 0

    if (between_person == "res") {
      if (fit_obj$singular) {
        slope_singular[i] <- TRUE
        if (slope_warn_singular)
          warning("Node '", vars[i], "': random-slope fit is singular, so at least ",
                  "one slope variance or covariance sits on the boundary. These ",
                  "slopes are not identified by the data; see $slope_singular and ",
                  "consider slope_structure = 'diag'.", call. = FALSE)
      }
      u  <- re[, 1L]
      h_vec[i] <- h_vec[i] + mean(u)
      person_shifts[, i] <- u - mean(u)
      person_sd[i] <- stats::sd(u)
      vc <- as.data.frame(lme4::VarCorr(fit_obj$fit))
      slope_varcorr[[i]] <- vc
      want <- paste0("lag_", nms)
      if (all(want %in% colnames(re))) {
        dev <- as.matrix(re[, want, drop = FALSE]); colnames(dev) <- nms
        # Centre exactly as the random intercept is centred, so the reported J is
        # the mean of the person-specific transition matrices and person_J
        # deviates from it with mean zero (sum_k u_kij = 0 for every edge).
        J_vec <- J_vec + colMeans(dev)
        dev <- sweep(dev, 2, colMeans(dev), "-")
        person_J[, i, nms] <- sweep(dev, 2, J_vec, "+")
        slope_sd[i, nms]   <- apply(dev, 2, stats::sd)
        for (j in seq_along(nms)) {
          hit <- vc$sdcor[vc$var1 == want[j] & (is.na(vc$var2) | vc$var1 == vc$var2)]
          if (length(hit) == 1L) slope_sd_est[i, nms[j]] <- hit
        }
      } else {
        warning("Node '", vars[i], "': could not map random-slope columns; ",
                "person_J for this node falls back to the common J.", call. = FALSE)
      }
    } else {
      u <- re[, 1L]
      h_vec[i] <- h_vec[i] + mean(u)
      person_shifts[, i] <- u - mean(u)
      person_sd[i] <- stats::sd(u)
    }
    J_mat[i, nms] <- J_vec
  }

  if (verbose) {
    n_const <- colSums(person_constant)
    if (any(n_const > 0L))
      message("  Constant-outcome person-node cells (intercept determined by ",
              "shrinkage towards the population, not by own data): ",
              paste(sprintf("%s %d/%d", names(n_const), n_const, n_p)[n_const > 0L],
                    collapse = ", "), ". See $person_constant.")
    if (length(node_zero_var) > 0L)
      message("  Zero-variance nodes: ", paste(node_zero_var, collapse = ", "),
              " (intercept = marginal logit, edges = 0).")
    if (length(node_fit_failed) > 0L)
      message("  Nodes whose multilevel fit failed: ", paste(node_fit_failed, collapse = ", "))
    if (!is.null(person_J)) {
      message("  Person-specific transition matrices: ", n_p, " x ", n_vars, " x ",
              n_vars, ". Use $person_J.")
      message("  Random-slope SD by node (row) and predictor (column), ",
              slope_structure, " structure:")
      print(round(slope_sd_est, 3))
      if (any(slope_singular))
        message("  Singular node(s): ", paste(names(which(slope_singular)),
                collapse = ", "), ": slope variances/covariances on the boundary.")
    }
  }

  # ==========================================================================
  # Section 5: spectral diagnostics (weighted Jacobian)
  # ==========================================================================
  # Linearised (heuristic) contraction rate of the mean-field map
  # m -> sigma(b + J m), evaluated at its fixed point rather than at the
  # all-zero state, for the typical person and for every person. The exact
  # relaxation time of the finite chain is -1/log|lambda_2(P)| (kimMixing()).
  rho_J <- max(abs(eigen(J_mat, only.values = TRUE)$values))
  rho_at <- function(Jk, bk) {
    m <- .kim_mf_fixed_point(Jk, bk)
    d <- m * (1 - m)
    max(abs(eigen(diag(d, length(d)) %*% Jk, only.values = TRUE)$values))
  }
  rho_star <- rho_at(J_mat, h_vec)
  rho_persons <- vapply(seq_len(n_p), function(k) {
    Jk <- if (is.null(person_J)) J_mat else matrix(person_J[k, , ], n_vars, n_vars)
    rho_at(Jk, h_vec + person_shifts[k, ])
  }, numeric(1))
  rho_star_max <- max(c(rho_star, rho_persons))
  mixing_time_approx <- if (rho_star_max < 1) max(1, 1 / (1 - rho_star_max)) else {
    warning("Weighted-Jacobian spectral radius >= 1 for at least one person; ",
            "the linearised heuristic gives no finite relaxation time. Use kimMixing().")
    Inf
  }

  if (verbose) {
    message("  Spectral radius rho(J)            = ", round(rho_J, 4), " (reference)")
    message("  rho_star typical / max over persons = ", round(rho_star, 4), " / ",
            round(rho_star_max, 4), " (linearised heuristic)")
    message("  Heuristic relaxation time         = ",
            if (is.finite(mixing_time_approx)) round(mixing_time_approx, 2) else "Inf",
            " transitions (exact value: kimMixing())")
  }

  # ==========================================================================
  # Section 6: between-subject network + return
  # ==========================================================================
  if (verbose) message("Estimating between-subject network...")
  between_est <- if (estimate_between) .estimate_between_kim(
    data_clean[as.character(data_clean[[id]]) %in% pids_clean, , drop = FALSE],
    vars, id, gamma, verbose) else
    list(network = NULL, dropped_nodes = character())
  between_net           <- between_est$network
  between_dropped_nodes <- between_est$dropped_nodes

  result <- list(
    J                     = J_mat,
    h                     = h_vec,
    temporal_network      = J_mat,
    thresholds            = h_vec,
    support               = support,
    stage1_beta           = beta1,
    stage1_se             = se1,
    screen                = screen,
    screen_threshold      = screen_threshold,
    between_network       = between_net,
    between_dropped_nodes = between_dropped_nodes,
    spectral_radius       = rho_star,
    spectral_radius_J     = rho_J,
    spectral_radius_star  = rho_star,
    spectral_radius_star_max = rho_star_max,
    mixing_time           = mixing_time_approx,
    transition_types      = transition_types,
    n_lag_pairs           = nrow(X_lag),
    n_na_removed          = total_na_removed,
    n_gap_removed         = total_gap_removed,
    pair_audit            = pairs$audit,
    pair_summary          = pairs$summary,
    pair_time_unit        = pairs$time_unit,
    excluded_post_clean   = excluded_post_clean,
    vars                  = vars,
    n_persons             = n_p,
    excluded              = excluded,
    node_fit_failed       = node_fit_failed,
    node_zero_var         = node_zero_var,
    between_person        = between_person,
    person_shifts         = person_shifts,     # n_p x p ; a_ki (mean 0 per node)
    person_sd             = person_sd,         # per-node SD of a_ki
    person_J              = person_J,          # n_p x p x p ; NULL unless "res"
    slope_sd              = slope_sd,          # p x p SD of the FITTED (shrunk) deviations
    slope_sd_est          = slope_sd_est,      # p x p model-estimated slope SD (use this one)
    slope_varcorr         = slope_varcorr,     # per-node VarCorr table
    slope_singular        = slope_singular,    # per-node singular-fit flag
    person_constant       = person_constant,   # n_p x p ; outcome never varies
    person_ids            = pids_clean,
    model                 = "kim",
    estimation_method     = paste0("screened conditional transition likelihood (",
                                   between_person, ", screen = ", screen, ")"),
    max_time_gap          = max_time_gap,
    call                  = match.call()
  )
  class(result) <- "dynaNIRA_kim"
  result$settings <- list(vars = vars, id = id, time = time, lags = lags,
    between_person = between_person, screen = screen, fdr_q = fdr_q, alpha = alpha,
    min_obs = min_obs, max_time_gap = max_time_gap, gamma = gamma, seed = seed,
    estimate_between = estimate_between, slope_structure = slope_structure,
    slope_warn_singular = slope_warn_singular, min_time_gap = min_time_gap,
    consec = consec, dayvar = dayvar, beepvar = beepvar,
    allow_overnight = allow_overnight, time_unit = time_unit,
    beeps_per_day = beeps_per_day)
  if (verbose) message("fitKIM complete (between_person = ", between_person,
                       ", screen = ", screen, ").")
  result
}
