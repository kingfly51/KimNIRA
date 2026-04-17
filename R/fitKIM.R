# =============================================================================
# fitKIM.R  ----  KIM (Kinetic Ising Model for dynamic network analysis)
#
# Model:
#   P(s(t+1) | s(t)) = prod_i sigma( sum_j J_ij * s_j(t) + h_i )
#
#   where:
#     J[i,j] = temporal influence of node j at t on node i at t+1
#              (directed, asymmetric)
#     h[i]   = intercept (threshold) for node i
#
#   Nodes at the current time point are mutually independent (no contemporaneous
#   network). Activation of each node at t+1 is determined solely by all nodes
#   at the previous time point t.
#
#   A between-subject network is also estimated.
#
# Estimation:
#   Within-person (temporal) network:
#     Pseudo-likelihood (Besag 1975)
#     Predictors = s(t), estimating the temporal network J
#
#   Between-person network: EBICglasso (Epskamp & Fried, 2018)
#     Regularised partial correlations on the matrix of person-level node means.
#
# References:
#   Besag, J. (1975). Statistical analysis of non-lattice data.
#     The Statistician, 24(3), 179-195.
#   Campajola et al. (2022). Scientific Reports, 12, 19339.
#   Epskamp, S. & Fried, E. I. (2018). A tutorial on regularized partial
#     correlation networks. Psychological Methods, 23(4), 617-634.
# =============================================================================


# -----------------------------------------------------------------------------
# Internal helper: between-subject network estimation
# -----------------------------------------------------------------------------
#
#' Estimate Between-Subject Network
#'
#' @description
#' Steps:
#'   1. Compute each person's node means (activation rate across all time points)
#'   2. Compute correlation matrix of between-person means
#'   3. Estimate regularised partial correlation network via EBICglasso
#'
#' @param data Cleaned data frame
#' @param vars Symptom variable names
#' @param id Subject ID column name
#' @param gamma EBIC hyperparameter (0 = BIC, 0.5 = recommended default)
#' @param verbose Whether to print messages
#'
#' @return A named list with two elements:
#'   \describe{
#'     \item{\code{network}}{Full p x p partial correlation matrix. Cells
#'       involving nodes that could not be estimated (zero between-person
#'       variance) are set to \code{NA} -- not \code{0} -- to distinguish
#'       "not estimable" from "estimated as zero". Returns \code{NULL} if
#'       estimation fails entirely.}
#'     \item{\code{dropped_nodes}}{Character vector of node names that were
#'       excluded from EBICglasso because their between-person variance was
#'       effectively zero. Length-0 character vector when all nodes were used.}
#'   }
#' @keywords internal
.estimate_between_kim <- function(data, vars, id, gamma = 0.5, verbose = TRUE) {

  # Helper: build a full-size NA-filled result for clean early returns
  .null_result <- function(dropped = character(0))
    list(network = NULL, dropped_nodes = dropped)

  if (!requireNamespace("qgraph", quietly = TRUE)) {
    warning("Package 'qgraph' is required to estimate the between-subject network: ",
            "install.packages('qgraph'). Returning NULL.")
    return(.null_result())
  }

  # 1. Compute each person's mean on each node
  person_ids <- unique(data[[id]])
  person_means <- do.call(rbind, lapply(person_ids, function(pid) {
    pdata <- data[data[[id]] == pid, vars, drop = FALSE]
    colMeans(pdata, na.rm = TRUE)
  }))
  rownames(person_means) <- person_ids
  colnames(person_means) <- vars

  n_persons <- nrow(person_means)
  n_vars    <- length(vars)

  if (n_persons < n_vars + 1L) {
    warning("Number of persons (", n_persons, ") is less than number of variables + 1 (",
            n_vars + 1L, "). Cannot stably estimate between-subject network. Returning NULL.")
    return(.null_result())
  }

  # 2. Identify and remove zero-variance nodes
  #    "Zero between-person variance" means all subjects have the same mean
  #    activation rate for that node -- no individual differences to correlate.
  #    These nodes are dropped from EBICglasso but recorded so downstream
  #    code can distinguish "not estimable" from "estimated as zero".
  col_var      <- apply(person_means, 2, var, na.rm = TRUE)
  dropped_vars <- vars[col_var < 1e-10]
  valid_vars   <- vars[col_var >= 1e-10]

  if (length(dropped_vars) > 0L) {
    if (verbose)
      message("  Between-subject network: ",
              length(dropped_vars), " node(s) dropped -- zero between-person variance",
              " (cannot estimate partial correlations): ",
              paste(dropped_vars, collapse = ", "),
              "\n  Edges involving these nodes will be NA in the returned matrix.")
    if (length(valid_vars) < 3L) {
      warning("Fewer than 3 estimable nodes remain. ",
              "Cannot estimate between-subject network. Returning NULL.")
      return(.null_result(dropped_vars))
    }
    person_means <- person_means[, valid_vars, drop = FALSE]
  }

  n_valid <- length(valid_vars)

  # 3. Compute correlation matrix
  cor_mat <- cor(person_means, use = "pairwise.complete.obs")

  if (any(is.na(cor_mat))) {
    warning("Correlation matrix contains NAs (possibly due to pairwise missing persons). ",
            "Affected node-pair partial correlations set to 0.")
    cor_mat[is.na(cor_mat)] <- 0
  }

  # Ensure positive definiteness
  eig_vals <- eigen(cor_mat, only.values = TRUE)$values
  if (any(eig_vals <= 0)) {
    if (verbose) message("  Correlation matrix is not positive definite; projecting via nearPD...")
    if (!requireNamespace("Matrix", quietly = TRUE)) {
      cor_mat <- cor_mat + diag(n_valid) * (abs(min(eig_vals)) + 0.001)
      D_inv   <- diag(1 / sqrt(diag(cor_mat)))
      cor_mat <- D_inv %*% cor_mat %*% D_inv
    } else {
      pd_result <- Matrix::nearPD(cor_mat, corr = TRUE, keepDiag = TRUE)
      cor_mat   <- as.matrix(pd_result$mat)
    }
    dimnames(cor_mat) <- list(valid_vars, valid_vars)
  }

  # 4. EBICglasso estimation
  between_net_small <- tryCatch(
    qgraph::EBICglasso(
      S      = cor_mat,
      n      = n_persons,
      gamma  = gamma,
      returnAllResults = FALSE
    ),
    error = function(e) {
      warning("EBICglasso estimation failed: ", e$message, ". Returning NULL.")
      NULL
    }
  )

  if (is.null(between_net_small)) return(.null_result(dropped_vars))

  # 5. Embed the estimated sub-matrix back into a full p x p matrix.
  #    Cells that involve dropped nodes are set to NA (not 0) because:
  #      - 0 would be misread as "no association detected by the model"
  #      - NA correctly signals "association could not be estimated"
  if (length(dropped_vars) > 0L) {
    between_net <- matrix(NA_real_, n_vars, n_vars,
                          dimnames = list(vars, vars))
    between_net[valid_vars, valid_vars] <- between_net_small
    # Diagonal of dropped nodes: set to NA (undefined, not a partial corr)
    # Diagonal of valid nodes: EBICglasso already returns 0 on diag
  } else {
    between_net <- between_net_small
    dimnames(between_net) <- list(vars, vars)
  }

  if (verbose) {
    # Count edges only among estimable pairs (ignore NA cells)
    n_edges <- sum(abs(between_net[valid_vars, valid_vars]) > 1e-6,
                   na.rm = TRUE) / 2L
    message("  Between-subject network: ", n_edges, " non-zero edge(s) ",
            "among ", length(valid_vars), " estimable node(s).")
  }

  list(network = between_net, dropped_nodes = dropped_vars)
}


# -----------------------------------------------------------------------------
# Main function
# -----------------------------------------------------------------------------

#' Fit a KIM Network Model
#'
#' @description
#' Fits a psychological network model based on a conditional Markov random field:
#' the standard dynamic Ising network (temporal J + between-subject).
#'
#' Nodes at the current time point are mutually independent (no contemporaneous
#' network). Activation of each node at t+1 is determined solely by all nodes
#' at t.
#'
#' @param data Long-format data frame, one row per EMA measurement.
#' @param vars Character vector of symptom variable names; must be binary \code{\{0, 1\}}.
#' @param id Character string, subject ID column name. Default \code{"id"}.
#' @param time Character string, time index column name. Default \code{"time"}.
#' @param lags Integer, number of lags. Default \code{1L}.
#' @param lambda_rule Character, \code{"lambda.1se"} or \code{"lambda.min"}.
#'   Default \code{"lambda.min"}.
#' @param n_folds Integer, number of cross-validation folds. Default \code{5L}.
#' @param penalize_ar Logical, whether to penalise the autoregressive term.
#'   Default \code{TRUE}.
#' @param use_offset Logical, whether to use logit(person mean) as offset.
#'   Default \code{FALSE}.
#' @param min_obs Integer, minimum number of valid lag pairs. Default \code{10L}.
#' @param max_time_gap Integer or NULL. Maximum allowed time gap. If the
#'   difference between consecutive time points exceeds this value the lag pair
#'   is treated as discontinuous and excluded. Default \code{NULL} (no check).
#' @param gamma Numeric, EBICglasso hyperparameter for between-subject network.
#'   Default \code{0.5}.
#' @param verbose Logical, whether to print progress messages. Default \code{TRUE}.
#'
#' @return A list of class \code{"dynaNIRA_kim"} containing (among others):
#' \describe{
#'   \item{\code{spectral_radius}}{Alias for \code{spectral_radius_star}
#'     (kept for backward compatibility).}
#'   \item{\code{spectral_radius_J}}{Raw spectral radius \eqn{\rho(J)} of the
#'     temporal edge matrix.  Provided for reference only; do \strong{not} use
#'     this to assess KIM stability -- it ignores the sigmoid compression.}
#'   \item{\code{spectral_radius_star}}{Spectral radius of the
#'     \emph{weighted Jacobian}
#'     \eqn{\rho(\mathrm{diag}(\sigma'(h)) J)}, the correct local-
#'     linearisation stability criterion for KIM.  Values below 1 indicate
#'     the chain contracts near the estimated stationary point.  Because
#'     \eqn{\sigma'(h) \le 0.25}, \code{rho_star} \eqn{\le 0.25 \times}
#'     \code{rho_J}.}
#'   \item{\code{mixing_time}}{Approximate mixing time (steps) derived from
#'     \eqn{1 / (1 - \rho^*)}. \code{Inf} when \code{rho_star >= 1}.
#'     Use this to inform the \code{n_steps} argument of
#'     \code{\link{runKIMNIRA}}.}
#'   \item{\code{between_network}}{Full p x p between-subject partial
#'     correlation matrix, or \code{NULL} if estimation failed entirely.
#'     Cells involving nodes listed in \code{between_dropped_nodes} are
#'     \code{NA} -- not \code{0} -- to indicate "not estimable" rather than
#'     "no association".}
#'   \item{\code{between_dropped_nodes}}{Character vector of node names
#'     excluded from between-subject estimation due to zero between-person
#'     variance. Length-0 vector when all nodes were estimable.}
#'   \item{\code{node_zero_var}}{Character vector of node names whose
#'     response vector had zero variance in the lag-pair dataset (always 0
#'     or always 1). These nodes have \code{h} set to marginal logit and
#'     outgoing temporal edges set to 0.}
#'   \item{\code{node_fit_failed}}{Character vector of nodes for which
#'     \code{cv.glmnet} raised an error. Same fallback as zero-variance nodes.}
#' }
#'
#' @examples
#' \dontrun{
#' fit <- fitKIM(data, vars = c("sad", "worry", "sleep"))
#' }
#' @importFrom stats cor median na.omit runif setNames plogis var
#' @export
fitKIM <- function(data,
                   vars,
                   id               = "id",
                   time             = "time",
                   lags             = 1L,
                   lambda_rule      = c("lambda.min", "lambda.1se"),
                   n_folds          = 5L,
                   penalize_ar      = TRUE,
                   use_offset       = FALSE,
                   min_obs          = 10L,
                   max_time_gap     = NULL,
                   gamma            = 0.5,
                   verbose          = TRUE) {

  lambda_rule <- match.arg(lambda_rule)

  # -- Dependency check --------------------------------------------------------
  if (!requireNamespace("glmnet", quietly = TRUE))
    stop("Package 'glmnet' is required: install.packages('glmnet')")

  # -- Input validation --------------------------------------------------------
  if (!is.data.frame(data))
    stop("`data` must be a data.frame.")

  miss <- setdiff(c(vars, id, time), names(data))
  if (length(miss) > 0)
    stop("Columns missing from data: ", paste(miss, collapse = ", "))

  non_binary <- vars[!vapply(vars, function(v)
    all(na.omit(unique(data[[v]])) %in% c(0, 1)), logical(1))]
  if (length(non_binary) > 0)
    stop("Non-binary variables: ", paste(non_binary, collapse = ", "))

  n_vars <- length(vars)

  if (verbose) {
    message("Model: KIM (temporal + between-subject)")
  }

  # -- Exclude subjects with insufficient observations -------------------------
  obs_counts <- tapply(seq_len(nrow(data)), data[[id]], length)
  excluded   <- names(obs_counts)[obs_counts < (min_obs + lags)]
  if (length(excluded) > 0 && verbose)
    message("Excluding ", length(excluded), " subject(s) with fewer than ",
            min_obs + lags, " observations.")
  data_clean <- data[!data[[id]] %in% excluded, ]
  n_persons  <- length(unique(data_clean[[id]]))
  if (n_persons == 0L) stop("No subjects with sufficient observations.")


  # ==========================================================================
  # Section 1: Build lag-pair matrices
  # ==========================================================================

  if (verbose) message("Building lag-pair matrices...")

  X_lag          <- NULL   # s(t)
  Y_t1           <- NULL   # s(t+1)
  person_id_pool <- NULL
  total_na_removed  <- 0L
  total_gap_removed <- 0L
  excluded_post_clean <- character(0)

  for (pid in unique(data_clean[[id]])) {
    pdata <- data_clean[data_clean[[id]] == pid, ]
    pdata <- pdata[order(pdata[[time]]), ]
    n_t   <- nrow(pdata)
    if (n_t <= lags) next

    X_p <- as.matrix(pdata[1:(n_t - lags),  vars])
    Y_p <- as.matrix(pdata[(1 + lags):n_t,  vars])

    if (!is.null(max_time_gap)) {
      time_vals  <- pdata[[time]]
      time_diffs <- diff(time_vals, lag = lags)
      stopifnot(length(time_diffs) == nrow(X_p))
      gap_rows <- time_diffs > max_time_gap
      total_gap_removed <- total_gap_removed + sum(gap_rows)
      if (any(gap_rows)) {
        X_p <- X_p[!gap_rows, , drop = FALSE]
        Y_p <- Y_p[!gap_rows, , drop = FALSE]
      }
      if (nrow(X_p) == 0L) next
    }

    na_rows <- apply(cbind(X_p, Y_p), 1, anyNA)
    total_na_removed <- total_na_removed + sum(na_rows)
    if (all(na_rows)) next
    X_p <- X_p[!na_rows, , drop = FALSE]
    Y_p <- Y_p[!na_rows, , drop = FALSE]
    if (nrow(X_p) < min_obs) {
      excluded_post_clean <- c(excluded_post_clean, as.character(pid))
      next
    }

    storage.mode(X_p) <- "double"
    storage.mode(Y_p) <- "double"

    X_lag          <- rbind(X_lag, X_p)
    Y_t1           <- rbind(Y_t1,  Y_p)
    person_id_pool <- c(person_id_pool,
                        rep(as.character(pid), nrow(X_p)))
  }

  if (verbose) {
    message("  Valid lag pairs: ", nrow(X_lag))
    if (total_na_removed > 0L)
      message("  Lag pairs removed due to missing values: ", total_na_removed)
    if (!is.null(max_time_gap) && total_gap_removed > 0L)
      message("  Lag pairs removed due to time-gap discontinuities: ", total_gap_removed)
  }
  if (is.null(X_lag) || nrow(X_lag) < n_vars * 2L)
    stop("Insufficient lag pairs.")
  if (verbose && length(excluded_post_clean) > 0L)
    message("  Subjects additionally excluded after cleaning (insufficient lag pairs): ",
            length(excluded_post_clean))


  # ==========================================================================
  # Section 2: Compute offsets
  # ==========================================================================

  pids_clean <- unique(person_id_pool)
  person_means_all <- do.call(rbind, lapply(pids_clean, function(pid) {
    pdata <- data_clean[data_clean[[id]] == pid, vars, drop = FALSE]
    colMeans(pdata, na.rm = TRUE)
  }))
  rownames(person_means_all) <- pids_clean
  colnames(person_means_all) <- vars
  eps              <- 0.01
  person_means_all <- pmax(pmin(person_means_all, 1 - eps), eps)

  compute_offset <- function(node_name) {
    p <- person_means_all[person_id_pool, node_name]
    log(p / (1 - p))
  }


  # ==========================================================================
  # Section 3: Pseudo-likelihood estimation of temporal network
  #            (node-wise logistic regression, lagged predictors only)
  # ==========================================================================

  if (verbose) {
    message("Pseudo-likelihood estimation of temporal network (predictors = lag only)...")
  }

  J_mat <- matrix(0, n_vars, n_vars, dimnames = list(vars, vars))
  h_vec           <- setNames(numeric(n_vars), vars)
  lambda_selected <- setNames(numeric(n_vars), vars)
  node_fit_failed <- character(0)
  node_zero_var   <- character(0)   # nodes skipped due to zero response variance

  # Subject-grouped cross-validation foldid to prevent data leakage
  set.seed(42)
  unique_pids  <- unique(person_id_pool)
  n_unique     <- length(unique_pids)
  actual_folds <- min(n_folds, n_unique)
  if (actual_folds < n_folds && verbose)
    message("  Number of subjects (", n_unique, ") < n_folds (", n_folds,
            "); using ", actual_folds, " folds.")
  pid_folds   <- sample(rep(1:actual_folds, length.out = n_unique))
  names(pid_folds) <- unique_pids
  foldid <- pid_folds[person_id_pool]

  for (i in seq_len(n_vars)) {

    if (verbose) message("  Node ", i, "/", n_vars, ": ", vars[i])

    y_i <- Y_t1[, i]

    if (is.na(var(y_i)) || var(y_i, na.rm = TRUE) < 1e-10) {
      # Node has zero variance in the response: all lag pairs produced the same
      # outcome (always 0 or always 1). Model fitting is undefined.
      # Set intercept to marginal logit -- identical to the cv_fit failure path
      # below -- so all fallback intercepts are computed by the same formula.
      node_zero_var <- c(node_zero_var, vars[i])
      p_marginal    <- mean(y_i, na.rm = TRUE)
      p_marginal    <- max(min(p_marginal, 1 - eps), eps)
      h_vec[i]      <- log(p_marginal / (1 - p_marginal))
      warning("Node '", vars[i], "' has zero response variance; skipping. ",
              "Intercept set to marginal logit = ", round(h_vec[i], 4),
              "; all outgoing temporal edge weights set to 0.")
      next
    }

    X_full <- X_lag
    colnames(X_full) <- paste0("lag_", vars)

    pf <- rep(1, ncol(X_full))
    if (!penalize_ar) pf[i] <- 0

    storage.mode(X_full) <- "double"

    offset_i <- if (use_offset) compute_offset(vars[i]) else NULL

    cv_fit <- tryCatch(
      glmnet::cv.glmnet(
        x              = X_full,
        y              = y_i,
        family         = "binomial",
        alpha          = 1,
        nfolds         = actual_folds,
        penalty.factor = pf,
        offset         = offset_i,
        type.measure   = "deviance",
        foldid         = foldid
      ),
      error = function(e) {
        warning("Node '", vars[i], "' estimation failed: ", e$message)
        NULL
      }
    )

    if (is.null(cv_fit)) {
      node_fit_failed <- c(node_fit_failed, vars[i])
      p_marginal <- mean(y_i, na.rm = TRUE)
      p_marginal <- max(min(p_marginal, 1 - eps), eps)
      h_vec[i]   <- log(p_marginal / (1 - p_marginal))
      warning("Node '", vars[i], "' estimation failed; intercept set to marginal logit = ",
              round(h_vec[i], 4), "; temporal edge weights set to 0.")
      next
    }

    lam   <- cv_fit[[lambda_rule]]
    coefs <- as.numeric(coef(cv_fit, s = lam))

    h_vec[i]           <- coefs[1L]
    J_mat[i, ]         <- coefs[2L:(n_vars + 1L)]
    lambda_selected[i] <- lam
  }

  if (length(node_zero_var) > 0L && verbose)
    message("  Zero-variance nodes (response always 0 or always 1): ",
            paste(node_zero_var, collapse = ", "),
            " -- intercept = marginal logit, outgoing edges = 0.")
  if (length(node_fit_failed) > 0 && verbose) {
    message("  Warning: the following nodes failed estimation; marginal logit used as intercept: ",
            paste(node_fit_failed, collapse = ", "))
  }


  # ==========================================================================
  # Section 4: Derived statistics
  # ==========================================================================

  # rho_J: raw spectral radius of J.
  # This is the linear-system quantity (VAR(1) stability criterion).
  # It is retained for reference but should NOT be used as the stability
  # criterion for KIM, because the sigmoid activation function compresses
  # all outputs into (0, 1) regardless of how large rho_J is.  Using rho_J
  # directly leads to systematic false-positive "unstable network" warnings
  # whenever rho_J > 1, even though the chain may mix very quickly.
  rho_J <- max(abs(eigen(J_mat, only.values = TRUE)$values))

  # rho_star: spectral radius of the weighted Jacobian diag(d) %*% J,
  # where d_i = sigma'(h_i) = sigma(h_i) * (1 - sigma(h_i)) in (0, 0.25].
  #
  # This is the LOCAL linearisation of the KIM update around the stationary
  # point implied by the estimated thresholds h.  It correctly accounts for
  # the sigmoid compression:
  #
  #   d/ds_j  E[s_i(t+1) | s(t)] = sigma'(eta_i) * J_ij = d_i * J_ij
  #
  # Stability criterion: rho_star < 1  (not rho_J < 1).
  # Because d_i <= 0.25 for all i, we always have rho_star <= 0.25 * rho_J,
  # so a network with rho_J = 2 can still have rho_star = 0.4 and mix fast.
  #
  # Mixing time approximation 1 / (1 - rho_star) is derived from the
  # geometric decay rate of autocorrelations in the linearised system;
  # it remains an approximation (the true mixing time of the nonlinear
  # chain requires the spectral gap of the full transition kernel), but
  # it is far more accurate than 1 / (1 - rho_J).
  s_bar    <- stats::plogis(h_vec)          # sigma(h), marginal activation prob
  d_bar    <- s_bar * (1 - s_bar)           # sigma'(h), in (0, 0.25]
  rho_star <- max(abs(eigen(diag(d_bar) %*% J_mat,
                            only.values = TRUE)$values))

  mixing_time_approx <- if (rho_star < 1) round(1 / (1 - rho_star)) else {
    warning(
      "Weighted Jacobian spectral radius >= 1 (rho_star = ",
      round(rho_star, 4), "). ",
      "The KIM chain may not mix reliably near the estimated stationary point. ",
      "Consider increasing n_steps in runKIMNIRA() or inspecting the J matrix ",
      "for unusually large edge weights."
    )
    Inf
  }

  if (verbose) {
    message("  Spectral radius rho(J)    = ", round(rho_J,    4),
            "  (raw J matrix; reference only)")
    message("  Spectral radius rho_star  = ", round(rho_star, 4),
            "  (weighted Jacobian; stability criterion)")
    message("  Mixing time approx.       = ",
            if (is.finite(mixing_time_approx)) mixing_time_approx else "> 1000",
            " steps  (based on rho_star; local linear approximation)")
  }


  # ==========================================================================
  # Section 5: Between-subject network
  # ==========================================================================

  if (verbose) message("Estimating between-subject network...")
  between_est           <- .estimate_between_kim(data_clean, vars, id, gamma, verbose)
  between_net           <- between_est$network        # p x p matrix or NULL; NA = not estimable
  between_dropped_nodes <- between_est$dropped_nodes  # nodes excluded due to zero between-person variance


  # ==========================================================================
  # Section 6: Assemble return object
  # ==========================================================================

  result <- list(
    J                       = J_mat,
    h                       = h_vec,

    temporal_network        = J_mat,
    thresholds              = h_vec,

    between_network         = between_net,
    between_dropped_nodes   = between_dropped_nodes,
    lambda_selected         = lambda_selected,

    # Spectral diagnostics
    # rho_J    : raw spectral radius of J; reference only (linear VAR(1))
    # rho_star : spectral radius of weighted Jacobian diag(sigma'(h)) %*% J;
    #            the correct stability criterion for KIM
    spectral_radius         = rho_star,       # primary: backward-compatible name
    spectral_radius_J       = rho_J,          # raw rho(J) for reference
    spectral_radius_star    = rho_star,        # explicit name
    mixing_time             = mixing_time_approx,
    mixing_time_approx      = mixing_time_approx,

    n_lag_pairs             = nrow(X_lag),
    n_na_removed            = total_na_removed,
    n_gap_removed           = total_gap_removed,
    excluded_post_clean     = excluded_post_clean,
    vars                    = vars,
    n_persons               = n_persons,
    excluded                = excluded,
    node_fit_failed         = node_fit_failed,
    node_zero_var           = node_zero_var,
    model                   = "kim",
    estimation_method       = "pseudo-likelihood (temporal only)",
    use_offset              = use_offset,
    penalize_ar             = penalize_ar,
    max_time_gap            = max_time_gap,
    call                    = match.call()
  )

  class(result) <- "dynaNIRA_kim"
  if (verbose) message("fitKIM complete.")
  result
}
