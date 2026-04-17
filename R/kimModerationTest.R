# =============================================================================
# kimModerationTest.R  ----  Moderation Test for KIM Temporal Network
#
# Purpose:
#   Tests whether any node k moderates temporal edges J[i,j] in the KIM.
#   This is a PREREQUISITE check before running the NodeIdentifyR Algorithm
#   (NIRA): the NodeIdentifyR Algorithm (NIRA) assumes that perturbing node k's threshold leaves all J[i,j]
#   unchanged.  If node k moderates the j -> i temporal relationship, then
#   altering k's intercept implicitly reshapes the edge weight, making
#   NodeIdentifyR Algorithm (NIRA) results theoretically unsound.
#
# Model tested for each (response node i, moderator candidate k):
#
#   Full model:
#     logit P(s_i(t+1) = 1) =
#       sum_j  J_ij  * s_j(t)          [main temporal edges]
#       + sum_{j != k}  M_ijk * s_j(t) * s_k(t)   [moderation terms]
#       + h_i
#
#   Null model (no moderation by k):
#     logit P(s_i(t+1) = 1) = sum_j J_ij * s_j(t) + h_i
#
# Inference:
#   Likelihood-ratio test (LRT) comparing full vs. null model, followed by
#   permutation-based p-value correction for each interaction coefficient.
#   A coefficient M_ijk is flagged as significant when its permutation
#   p-value (corrected for multiple comparisons) falls below alpha.
#
# References:
#   Besag, J. (1975). Statistical analysis of non-lattice data.
#     The Statistician, 24(3), 179-195.
#   Haslbeck, J. M. B., & Waldorp, L. J. (2020). mgm: Estimating
#     time-varying mixed graphical models in high-dimensional data.
#     Journal of Statistical Software, 93(8), 1-46.
# =============================================================================


# -----------------------------------------------------------------------------
# Internal helpers
# -----------------------------------------------------------------------------

#' Build lag-pair matrices from a fitted KIM object
#'
#' Re-extracts X_lag / Y_t1 from the original data using the same exclusion
#' criteria as fitKIM().  Exported so users can inspect the lag pairs if needed.
#'
#' @param data Original long-format data frame passed to fitKIM().
#' @param kim_fit A \code{"dynaNIRA_kim"} object from \code{\link{fitKIM}}.
#' @param max_time_gap Passed through from kim_fit if not overridden.
#'
#' @return A list with \code{X_lag} (n_pairs x p), \code{Y_t1} (n_pairs x p),
#'   and \code{person_id_pool} (length n_pairs character vector).
#' @keywords internal
.kim_build_lag_pairs <- function(data, kim_fit) {

  vars         <- kim_fit$vars
  id_col       <- as.character(kim_fit$call$id  %||% "id")
  time_col     <- as.character(kim_fit$call$time %||% "time")
  lags         <- as.integer(kim_fit$call$lags   %||% 1L)
  min_obs      <- as.integer(kim_fit$call$min_obs %||% 10L)
  max_time_gap <- kim_fit$max_time_gap
  excluded     <- c(kim_fit$excluded, kim_fit$excluded_post_clean)

  # Try to recover id / time column names robustly
  # (call$id may be a symbol, not a string)
  safe_name <- function(x, default) {
    if (is.null(x)) return(default)
    x <- tryCatch(as.character(x), error = function(e) default)
    if (length(x) != 1L || !nzchar(x)) return(default)
    x
  }
  id_col   <- safe_name(kim_fit$call$id,   "id")
  time_col <- safe_name(kim_fit$call$time, "time")

  data_clean <- data[!data[[id_col]] %in% excluded, ]

  X_lag <- NULL; Y_t1 <- NULL; person_id_pool <- NULL

  for (pid in unique(data_clean[[id_col]])) {
    pdata <- data_clean[data_clean[[id_col]] == pid, ]
    pdata <- pdata[order(pdata[[time_col]]), ]
    n_t   <- nrow(pdata)
    if (n_t <= lags) next

    X_p <- as.matrix(pdata[1:(n_t - lags),       vars])
    Y_p <- as.matrix(pdata[(1L + lags):n_t,      vars])

    if (!is.null(max_time_gap)) {
      time_vals  <- pdata[[time_col]]
      time_diffs <- diff(time_vals, lag = lags)
      gap_rows   <- time_diffs > max_time_gap
      if (any(gap_rows)) {
        X_p <- X_p[!gap_rows, , drop = FALSE]
        Y_p <- Y_p[!gap_rows, , drop = FALSE]
      }
      if (nrow(X_p) == 0L) next
    }

    na_rows <- apply(cbind(X_p, Y_p), 1, anyNA)
    if (all(na_rows)) next
    X_p <- X_p[!na_rows, , drop = FALSE]
    Y_p <- Y_p[!na_rows, , drop = FALSE]
    if (nrow(X_p) < min_obs) next

    storage.mode(X_p) <- "double"
    storage.mode(Y_p) <- "double"

    X_lag          <- rbind(X_lag, X_p)
    Y_t1           <- rbind(Y_t1,  Y_p)
    person_id_pool <- c(person_id_pool,
                        rep(as.character(pid), nrow(X_p)))
  }

  list(X_lag = X_lag, Y_t1 = Y_t1, person_id_pool = person_id_pool)
}

# Null-coalescing operator (avoid rlang dependency)


#' Single-node moderation test (internal)
#'
#' For response node i and moderator candidate k, compares:
#'   Full model  : y_i ~ X_lag + X_lag[, j] * X_lag[, k]  for all j != k
#'   Null model  : y_i ~ X_lag
#' via likelihood-ratio test and per-coefficient permutation test.
#'
#' @param y_i        Binary response vector (length n_pairs).
#' @param X_lag      Predictor matrix (n_pairs x p), columns named by vars.
#' @param k          Integer index of the moderator column in X_lag.
#' @param vars       Character vector of node names.
#' @param n_perm     Number of permutations for coefficient-level p-values.
#' @param alpha      Significance threshold (before correction).
#' @param method     Multiple-comparison correction method (stats::p.adjust).
#' @param foldid     Cross-validation fold IDs (not used here; glm-based).
#'
#' @return A list: lrt_p (LRT p-value), coef_df (data frame of interaction
#'   coefficients with permutation p-values and adjusted p-values).
#' @keywords internal
.kim_moderation_node <- function(y_i, X_lag, k, vars, n_perm, alpha,
                                 method) {

  p      <- ncol(X_lag)
  n      <- length(y_i)
  k_name <- vars[k]

  # Skip if response has zero variance
  if (var(y_i, na.rm = TRUE) < 1e-10)
    return(list(lrt_p = NA_real_, coef_df = NULL, skipped = TRUE))

  # Skip if moderator column has zero variance (always 0 or always 1)
  if (var(X_lag[, k], na.rm = TRUE) < 1e-10)
    return(list(lrt_p = NA_real_, coef_df = NULL, skipped = TRUE))

  # -- Null model: main effects only -----------------------------------------
  null_fit <- tryCatch(
    suppressWarnings(
      stats::glm(y_i ~ X_lag, family = stats::binomial(), maxit = 100)
    ),
    error = function(e) NULL
  )
  if (is.null(null_fit)) return(list(lrt_p = NA_real_, coef_df = NULL,
                                     skipped = TRUE))

  # -- Interaction terms: s_j(t) * s_k(t) for all j != k -------------------
  # Each column is named "int_<j>_x_<k>"
  j_idx      <- setdiff(seq_len(p), k)
  int_names  <- paste0("int_", vars[j_idx], "_x_", k_name)
  X_int      <- X_lag[, j_idx, drop = FALSE] * X_lag[, k]
  colnames(X_int) <- int_names

  X_full <- cbind(X_lag, X_int)

  # -- Full model -------------------------------------------------------------
  full_fit <- tryCatch(
    suppressWarnings(
      stats::glm(y_i ~ X_full, family = stats::binomial(), maxit = 100)
    ),
    error = function(e) NULL
  )
  if (is.null(full_fit)) return(list(lrt_p = NA_real_, coef_df = NULL,
                                     skipped = TRUE))

  # -- Likelihood-ratio test (node-level omnibus) ----------------------------
  lrt    <- tryCatch(stats::anova(null_fit, full_fit, test = "Chisq"),
                     error = function(e) NULL)
  lrt_p  <- if (!is.null(lrt)) lrt[2, "Pr(>Chi)"] else NA_real_

  # -- Extract observed interaction coefficients -----------------------------
  full_coef  <- stats::coef(full_fit)
  int_coef   <- full_coef[int_names]     # named vector, length p-1

  # -- Permutation test for each interaction coefficient ---------------------
  # Null distribution: permute s_k(t) (break the moderator-response link)
  # and refit the full model to get the permuted coefficient.
  # This preserves the marginal distribution of s_k(t) while destroying
  # any moderation signal.

  perm_mat <- matrix(NA_real_, nrow = n_perm, ncol = length(j_idx))
  colnames(perm_mat) <- int_names

  for (b in seq_len(n_perm)) {
    # Permute the moderator column only
    X_lag_perm        <- X_lag
    X_lag_perm[, k]   <- sample(X_lag[, k])
    X_int_perm        <- X_lag_perm[, j_idx, drop = FALSE] * X_lag_perm[, k]
    colnames(X_int_perm) <- int_names
    X_full_perm       <- cbind(X_lag_perm, X_int_perm)

    fit_perm <- tryCatch(
      suppressWarnings(
        stats::glm(y_i ~ X_full_perm, family = stats::binomial(), maxit = 50)
      ),
      error = function(e) NULL
    )
    if (!is.null(fit_perm)) {
      pc <- stats::coef(fit_perm)
      perm_mat[b, ] <- pc[int_names]
    }
  }

  # Two-sided permutation p-value: proportion of |perm| >= |obs|
  perm_p <- vapply(seq_along(j_idx), function(q) {
    obs_abs  <- abs(int_coef[q])
    perm_abs <- abs(perm_mat[, q])
    valid    <- !is.na(perm_abs)
    if (sum(valid) < 10L) return(NA_real_)
    (sum(perm_abs[valid] >= obs_abs, na.rm = TRUE) + 1L) /
      (sum(valid) + 1L)
  }, numeric(1))

  perm_p_adj <- stats::p.adjust(perm_p, method = method)

  coef_df <- data.frame(
    response_node   = vars[attr(y_i, "node_idx") %||%
                           which(colnames(X_lag) == colnames(X_lag)[1])],
    moderator_node  = k_name,
    edge_from       = vars[j_idx],
    edge_to         = NA_character_,   # filled by caller
    int_term        = int_names,
    coef_obs        = as.numeric(int_coef),
    perm_p          = perm_p,
    perm_p_adj      = perm_p_adj,
    significant     = !is.na(perm_p_adj) & perm_p_adj < alpha,
    stringsAsFactors = FALSE
  )

  list(lrt_p = lrt_p, coef_df = coef_df, skipped = FALSE)
}


# -----------------------------------------------------------------------------
# Main exported function
# -----------------------------------------------------------------------------

#' Test Moderation Effects in the KIM Temporal Network
#'
#' @description
#' Tests whether any node \eqn{k} acts as a \emph{moderator} of temporal
#' edges \eqn{J_{ij}} in the fitted KIM.  This is a \strong{prerequisite}
#' check before running \code{\link{runKIMNIRA}}: NIRA assumes that
#' perturbing node \eqn{k}'s threshold leaves all edge weights unchanged.
#' If \eqn{k} moderates the \eqn{j \to i} relationship (i.e. the
#' interaction term \eqn{M_{ijk} \neq 0}), then the NodeIdentifyR Algorithm (NIRA) simulation
#' implicitly reshapes that edge, making results theoretically unsound.
#'
#' \strong{Model compared for each (response node \eqn{i}, moderator \eqn{k}):}
#'
#' \emph{Null:}
#' \deqn{\mathrm{logit}\,P(s_i(t+1)=1) = \sum_j J_{ij}\,s_j(t) + h_i}
#'
#' \emph{Full:}
#' \deqn{\mathrm{logit}\,P(s_i(t+1)=1) = \sum_j J_{ij}\,s_j(t)
#'   + \sum_{j \neq k} M_{ijk}\,s_j(t)\,s_k(t) + h_i}
#'
#' \strong{Inference:}
#' \enumerate{
#'   \item An omnibus likelihood-ratio test (LRT) assesses whether the
#'     full model significantly improves fit over the null model for each
#'     (i, k) pair.
#'   \item For each individual interaction coefficient \eqn{M_{ijk}}, a
#'     permutation p-value is computed by permuting the moderator column
#'     \eqn{s_k(t)} (destroying the moderator-outcome link while
#'     preserving marginal distributions).
#'   \item Permutation p-values are corrected for multiple comparisons
#'     across all (i, j, k) triples using \code{method}.
#' }
#'
#' \strong{Interpretation:}
#' \itemize{
#'   \item An \strong{empty} \code{significant_moderators} list means no
#'     significant moderation was detected -- the NodeIdentifyR Algorithm (NIRA) can proceed.
#'   \item A \strong{non-empty} list identifies which nodes moderate which
#'     temporal edges.  NodeIdentifyR Algorithm (NIRA) results should be interpreted with caution, or
#'     the analysis should be re-specified to account for moderation.
#' }
#'
#' @param data The \strong{original} long-format data frame passed to
#'   \code{\link{fitKIM}} (same object, not pre-processed).
#' @param kim_fit A \code{"dynaNIRA_kim"} object from \code{\link{fitKIM}}.
#' @param moderator_nodes Integer vector or \code{NULL}.  Indices of nodes
#'   to test as moderators.  \code{NULL} (default) tests all \eqn{p} nodes.
#' @param n_perm Positive integer.  Permutations per (i, k) pair for
#'   coefficient-level p-values.  Default \code{500L}.
#'   \itemize{
#'     \item \code{500}:  acceptable for exploratory screening.
#'     \item \code{2000}: recommended for publication.
#'   }
#' @param alpha Numeric significance threshold applied \emph{after}
#'   multiple-comparison correction.  Default \code{0.05}.
#' @param method Character.  Multiple-comparison correction method passed to
#'   \code{\link[stats]{p.adjust}}.  Default \code{"BH"} (Benjamini-Hochberg).
#'   Use \code{"bonferroni"} for conservative family-wise control.
#' @param correct_across Character, one of \code{"moderator"} or \code{"global"}.
#'   \describe{
#'     \item{\code{"moderator"}}{Correction applied separately within each
#'       moderator candidate \eqn{k}: p-values for all (i, j) pairs sharing
#'       the same \eqn{k} are corrected together.  Less conservative.}
#'     \item{\code{"global"}}{All interaction p-values across every
#'       (i, j, k) triple are corrected in a single pool.  Most conservative;
#'       recommended when a single omnibus statement is needed.}
#'   }
#'   Default \code{"moderator"}.
#' @param seed Integer random seed for reproducibility.  Default \code{2025L}.
#' @param verbose Logical, print progress.  Default \code{TRUE}.
#'
#' @return A named list of class \code{"kim_moderation_result"}:
#' \describe{
#'   \item{\code{summary_table}}{Data frame with one row per (response i,
#'     moderator k, edge j) triple.  Columns: \code{response_node},
#'     \code{moderator_node}, \code{edge_from}, \code{coef_obs},
#'     \code{perm_p}, \code{perm_p_adj}, \code{significant}.}
#'   \item{\code{lrt_table}}{Data frame with one row per (i, k) pair.
#'     Columns: \code{response_node}, \code{moderator_node}, \code{lrt_p},
#'     \code{lrt_p_adj}, \code{lrt_significant}.}
#'   \item{\code{significant_moderators}}{Character vector of node names
#'     that show at least one significant moderation effect.  Empty character
#'     vector means the NodeIdentifyR Algorithm (NIRA) is appropriate.}
#'   \item{\code{significant_triples}}{Data frame of only the significant
#'     (i, j, k) rows from \code{summary_table}.  Zero rows = no moderation.}
#'   \item{\code{params}}{List of call parameters for reproducibility.}
#' }
#'
#' @examples
#' \dontrun{
#' fit <- fitKIM(data, vars = c("dep","anx","sleep","fatigue"),
#'               id = "id", time = "time")
#'
#' # Test all nodes as moderators (exploratory)
#' mod_test <- kimModerationTest(data, fit, n_perm = 500)
#' print(mod_test)
#'
#' # Check outcome: empty = NodeIdentifyR Algorithm (NIRA) is appropriate
#' mod_test$significant_moderators
#'
#' # Test only specific nodes
#' mod_test2 <- kimModerationTest(data, fit,
#'                                moderator_nodes = c(1, 3),
#'                                n_perm = 2000,
#'                                method = "bonferroni")
#' }
#'
#' @seealso \code{\link{fitKIM}}, \code{\link{runKIMNIRA}},
#'   \code{\link{print.kim_moderation_result}}
#' @importFrom stats glm binomial anova p.adjust coef var
#' @export
kimModerationTest <- function(data,
                              kim_fit,
                              moderator_nodes = NULL,
                              n_perm          = 500L,
                              alpha           = 0.05,
                              method          = "BH",
                              correct_across  = c("moderator", "global"),
                              seed            = 2025L,
                              verbose         = TRUE) {

  # -- Argument validation --------------------------------------------------
  if (!inherits(kim_fit, "dynaNIRA_kim"))
    stop("`kim_fit` must be produced by fitKIM() (class 'dynaNIRA_kim').")
  if (!is.data.frame(data) && !is.matrix(data))
    stop("`data` must be a data.frame or matrix.")

  correct_across <- match.arg(correct_across)

  valid_methods <- c("holm","hochberg","hommel","bonferroni",
                     "BH","BY","fdr","none")
  if (!method %in% valid_methods)
    stop("`method` must be one of: ", paste(valid_methods, collapse = ", "))

  n_perm <- as.integer(n_perm)
  if (n_perm < 50L)
    warning("`n_perm` < 50 will yield very unstable permutation p-values.")

  vars   <- kim_fit$vars
  p      <- length(vars)

  if (is.null(moderator_nodes)) {
    moderator_nodes <- seq_len(p)
  } else {
    moderator_nodes <- as.integer(moderator_nodes)
    bad <- moderator_nodes[moderator_nodes < 1L | moderator_nodes > p]
    if (length(bad) > 0L)
      stop("moderator_nodes indices out of range [1, ", p, "]: ",
           paste(bad, collapse = ", "))
  }

  if (verbose) {
    message("-- KIM Temporal Network Moderation Test --")
    message("  Nodes (p)         : ", p)
    message("  Moderators tested : ", length(moderator_nodes),
            " (", paste(vars[moderator_nodes], collapse = ", "), ")")
    message("  n_perm            : ", n_perm)
    message("  Correction        : method = ", method,
            ", scope = ", correct_across)
    message("  alpha             : ", alpha)
  }

  # -- Rebuild lag-pair matrices --------------------------------------------
  if (verbose) message("  Rebuilding lag-pair matrices from original data...")
  lp <- .kim_build_lag_pairs(data, kim_fit)

  if (is.null(lp$X_lag) || nrow(lp$X_lag) < p * 2L)
    stop("Insufficient lag pairs reconstructed from data. ",
         "Ensure `data` is the same object passed to fitKIM().")

  X_lag <- lp$X_lag
  Y_t1  <- lp$Y_t1
  colnames(X_lag) <- vars
  colnames(Y_t1)  <- vars
  n_pairs <- nrow(X_lag)

  if (verbose) message("  Lag pairs available: ", n_pairs)

  # -- Main loop: for each moderator k, test all response nodes i -----------
  set.seed(seed)

  all_coef_rows <- list()   # accumulates coefficient-level results
  lrt_rows      <- list()   # accumulates LRT results

  n_k <- length(moderator_nodes)
  for (ki in seq_along(moderator_nodes)) {

    k      <- moderator_nodes[ki]
    k_name <- vars[k]

    if (verbose)
      message("  [", ki, "/", n_k, "] Moderator: ", k_name)

    # Check moderator variance across all lag pairs
    if (var(X_lag[, k], na.rm = TRUE) < 1e-10) {
      if (verbose)
        message("    Skipped -- moderator node '", k_name,
                "' has zero variance in lag-pair matrix.")
      next
    }

    for (i in seq_len(p)) {

      y_i <- Y_t1[, i]
      # Attach node index as attribute for use inside helper
      attr(y_i, "node_idx") <- i

      res <- .kim_moderation_node(
        y_i    = y_i,
        X_lag  = X_lag,
        k      = k,
        vars   = vars,
        n_perm = n_perm,
        alpha  = alpha,
        method = method
      )

      # LRT row
      lrt_rows[[length(lrt_rows) + 1L]] <- data.frame(
        response_node  = vars[i],
        moderator_node = k_name,
        lrt_p          = res$lrt_p,
        skipped        = isTRUE(res$skipped),
        stringsAsFactors = FALSE
      )

      if (!isTRUE(res$skipped) && !is.null(res$coef_df)) {
        df <- res$coef_df
        df$response_node <- vars[i]
        df$edge_to       <- vars[i]
        all_coef_rows[[length(all_coef_rows) + 1L]] <- df
      }
    }
  }

  # -- Assemble tables ------------------------------------------------------
  lrt_table <- do.call(rbind, lrt_rows)
  rownames(lrt_table) <- NULL

  # LRT multiple-comparison correction (across all non-skipped (i,k) pairs)
  valid_lrt   <- !is.na(lrt_table$lrt_p) & !lrt_table$skipped
  lrt_p_adj   <- rep(NA_real_, nrow(lrt_table))
  if (any(valid_lrt))
    lrt_p_adj[valid_lrt] <- stats::p.adjust(lrt_table$lrt_p[valid_lrt],
                                             method = method)
  lrt_table$lrt_p_adj     <- lrt_p_adj
  lrt_table$lrt_significant <- !is.na(lrt_p_adj) & lrt_p_adj < alpha

  if (length(all_coef_rows) == 0L) {
    # All nodes skipped
    summary_table <- data.frame(
      response_node  = character(0), moderator_node = character(0),
      edge_from      = character(0), edge_to        = character(0),
      int_term       = character(0), coef_obs       = numeric(0),
      perm_p         = numeric(0),   perm_p_adj     = numeric(0),
      significant    = logical(0),   stringsAsFactors = FALSE
    )
  } else {
    summary_table <- do.call(rbind, all_coef_rows)
    rownames(summary_table) <- NULL

    # -- Multiple-comparison correction on permutation p-values -------------
    if (correct_across == "global") {

      # Single pool across all (i, j, k) triples
      valid <- !is.na(summary_table$perm_p)
      adj   <- rep(NA_real_, nrow(summary_table))
      if (any(valid))
        adj[valid] <- stats::p.adjust(summary_table$perm_p[valid],
                                      method = method)
      summary_table$perm_p_adj <- adj

    } else {
      # Correct separately within each moderator k
      adj <- rep(NA_real_, nrow(summary_table))
      for (k_name in vars[moderator_nodes]) {
        rows  <- summary_table$moderator_node == k_name
        valid <- rows & !is.na(summary_table$perm_p)
        if (any(valid))
          adj[valid] <- stats::p.adjust(summary_table$perm_p[valid],
                                        method = method)
      }
      summary_table$perm_p_adj <- adj
    }

    summary_table$significant <- !is.na(summary_table$perm_p_adj) &
                                   summary_table$perm_p_adj < alpha
  }

  # -- Identify significant moderators -------------------------------------
  sig_triples <- summary_table[!is.na(summary_table$significant) &
                                  summary_table$significant, ]
  sig_mods    <- unique(sig_triples$moderator_node)

  if (verbose) {
    message("  Done.")
    message("  Significant interaction triples (perm p_adj < ", alpha, "): ",
            nrow(sig_triples))
    if (length(sig_mods) > 0L) {
      message("  Significant moderator nodes: ",
              paste(sig_mods, collapse = ", "))
      message("  => NodeIdentifyR Algorithm (NIRA) should be interpreted with caution for these nodes.")
    } else {
      message("  No significant moderation detected.",
              " NodeIdentifyR Algorithm (NIRA) assumption is supported.")
    }
  }

  # -- Return ----------------------------------------------------------------
  result <- list(
    summary_table           = summary_table,
    lrt_table               = lrt_table,
    significant_moderators  = sig_mods,
    significant_triples     = sig_triples,
    params = list(
      moderator_nodes = moderator_nodes,
      n_perm          = n_perm,
      alpha           = alpha,
      method          = method,
      correct_across  = correct_across,
      seed            = seed,
      n_pairs         = n_pairs,
      vars            = vars
    )
  )
  class(result) <- "kim_moderation_result"
  result
}


# -----------------------------------------------------------------------------
# S3 print method
# -----------------------------------------------------------------------------

#' @export
print.kim_moderation_result <- function(x, ...) {

  p <- x$params
  cat("\n-- KIM Temporal Network Moderation Test --\n\n")
  cat("  Nodes tested as moderators :",
      length(p$moderator_nodes), "\n")
  cat("  Lag pairs used             :", p$n_pairs, "\n")
  cat("  Permutations               :", p$n_perm, "\n")
  cat("  Correction                 :", p$method,
      "(scope =", p$correct_across, ")\n")
  cat("  Alpha                      :", p$alpha, "\n\n")

  # LRT summary
  cat("  Omnibus LRT (per response x moderator pair):\n")
  lrt <- x$lrt_table[!x$lrt_table$skipped, ]
  n_lrt_sig <- sum(lrt$lrt_significant, na.rm = TRUE)
  cat("    Significant pairs (LRT p_adj <", p$alpha, "):",
      n_lrt_sig, "/", nrow(lrt), "\n\n")

  # Coefficient-level summary
  n_sig <- nrow(x$significant_triples)
  cat("  Significant interaction triples (perm p_adj <",
      p$alpha, "):", n_sig, "\n\n")

  if (n_sig == 0L) {
    cat("  Result: NO significant moderation detected.\n")
    cat("          NodeIdentifyR Algorithm (NIRA) assumption (constant edge weights) is SUPPORTED.\n")
    cat("          => runKIMNIRA() can proceed.\n\n")
  } else {
    cat("  Result: Significant moderation detected.\n")
    cat("          NodeIdentifyR Algorithm (NIRA) results should be interpreted with caution\n")
    cat("          for the following moderator node(s):\n")
    cat("          ", paste(x$significant_moderators, collapse = ", "), "\n\n")
    cat("  Significant triples (sorted by perm_p_adj):\n")

    out <- x$significant_triples[
      order(x$significant_triples$perm_p_adj), ,
      drop = FALSE
    ]
    out <- out[, c("moderator_node", "edge_from", "edge_to",
                   "coef_obs", "perm_p", "perm_p_adj"), drop = FALSE]
    out$coef_obs    <- round(out$coef_obs,    4)
    out$perm_p      <- round(out$perm_p,      4)
    out$perm_p_adj  <- round(out$perm_p_adj,  4)
    print(out, row.names = FALSE)
    cat("\n  Interpretation of coef_obs:\n")
    cat("    M_ijk > 0: moderator node STRENGTHENS the j->i edge\n")
    cat("    M_ijk < 0: moderator node WEAKENS  the j->i edge\n\n")
  }

  invisible(x)
}


# -----------------------------------------------------------------------------
# Convenience summary plot
# -----------------------------------------------------------------------------

#' Plot Moderation Test Results as a Heatmap
#'
#' @description
#' Draws a heatmap of adjusted permutation p-values for all tested
#' (moderator, edge) combinations.  Cells with significant moderation
#' (p_adj < alpha) are highlighted.
#'
#' @param x A \code{"kim_moderation_result"} object.
#' @param type Character, \code{"coef"} (default) colours cells by the
#'   observed interaction coefficient; \code{"pvalue"} colours by
#'   \code{-log10(perm_p_adj)}.
#' @param ... Currently ignored.
#'
#' @return A \code{ggplot2} object (invisibly).
#' @export
plot.kim_moderation_result <- function(x,
                                       type = c("coef", "pvalue"),
                                       ...) {

  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required for plotting.")

  type <- match.arg(type)
  st   <- x$summary_table
  if (nrow(st) == 0L) {
    message("No coefficient results to plot (all nodes may have been skipped).")
    return(invisible(NULL))
  }

  st$edge_label <- paste0(st$edge_from, "\u2192", st$edge_to)
  st$sig_label  <- ifelse(!is.na(st$significant) & st$significant, "*", "")

  fill_var <- if (type == "coef") "coef_obs" else {
    st$neg_log_p <- -log10(pmax(st$perm_p_adj, 1e-4))
    "neg_log_p"
  }

  fill_label <- if (type == "coef") "Interaction\ncoefficient"  else
    "-log10(p_adj)"

  gg <- ggplot2::ggplot(
    st,
    ggplot2::aes(x = .data$moderator_node,
                 y = .data$edge_label,
                 fill = .data[[fill_var]])
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.4) +
    ggplot2::geom_text(ggplot2::aes(label = .data$sig_label),
                       color = "black", size = 5, vjust = 0.75) +
    ggplot2::scale_fill_gradient2(
      low      = "#2166AC",
      mid      = "white",
      high     = "#D6604D",
      midpoint = 0,
      name     = fill_label
    ) +
    ggplot2::labs(
      x        = "Moderator Node",
      y        = "Temporal Edge (j \u2192 i)",
      title    = "KIM Temporal Network -- Moderation Test",
      subtitle = paste0("* = p_adj < ", x$params$alpha,
                        " (", x$params$method, " correction",
                        ", scope = ", x$params$correct_across, ")")
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      axis.text.x      = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid       = ggplot2::element_blank(),
      legend.position  = "right"
    )

  print(gg)
  invisible(gg)
}
