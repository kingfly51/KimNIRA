# -----------------------------------------------------------------------------
# kimPermutationNIRAtest() -- Batch permutation + multiple-comparison correction
# -----------------------------------------------------------------------------

#' Batch Paired Permutation Tests Across All KIM NodeIdentifyR Algorithm (NIRA) Nodes
#'
#' @description
#' Applies \code{\link{kimPermutationTest}} to every intervention node in a
#' \code{runKIMNIRA()} result, corrects for multiple comparisons, and
#' produces both a statistics table and a plotting-ready data frame
#' consumable by \code{\link{kimPlotNIRA}}.
#'
#' This is the paired-design analogue of
#' \code{NIRApost::permutationNIRAtest}.
#'
#' @param nira_result A \code{"kim_nira_result"} object from
#'   \code{\link{runKIMNIRA}}.
#' @param nPerm Positive integer, permutations per node. Default \code{5000}.
#' @param method Correction method for \code{\link[stats]{p.adjust}}:
#'   \code{"bonferroni"} (default), \code{"BH"}, \code{"holm"}, etc.
#' @param seed Integer random seed. Default \code{2025L}.
#' @param verbose Logical. Default \code{TRUE}.
#'
#' @return A named \code{list} of class \code{"kim_perm_nira_result"}:
#' \describe{
#'   \item{\code{stat}}{Data frame, one row per node. Columns:
#'     \code{mean_other}, \code{sd_other}, \code{se_other}, \code{ci_other},
#'     \code{ciLower_other}, \code{ciUpper_other}, \code{cohen_d},
#'     \code{p}, \code{p.adjust}.}
#'   \item{\code{plot_data}}{Data frame with columns: \code{mean},
#'     \code{ciLower}, \code{ciUpper}, \code{id}, \code{node}.
#'     First row = \code{"baseline"}. Ready for \code{kimPlotNIRA()}.}
#'   \item{\code{nira_result}}{Reference to input.}
#' }
#'
#' @examples
#' \dontrun{
#' res <- runKIMNIRA(kim_fit, intervention_type = "alleviating")
#' stat_result <- kimPermutationNIRAtest(res)
#' print(stat_result)
#' kimPlotNIRA(stat_result, "alleviating")
#' }
#'
#' @seealso \code{\link{kimPermutationTest}}, \code{\link{kimPlotNIRA}}
#' @importFrom stats p.adjust
#' @export
kimPermutationNIRAtest <- function(nira_result,
                                   nPerm   = 5000L,
                                   method  = "bonferroni",
                                   seed    = 2025L,
                                   verbose = TRUE) {

  if (!inherits(nira_result, "kim_nira_result"))
    stop("`nira_result` must be produced by runKIMNIRA() ",
         "(class = 'kim_nira_result').")

  valid_m <- c("holm", "hochberg", "hommel", "bonferroni",
               "BH", "BY", "fdr", "none")
  if (!method %in% valid_m)
    stop("Invalid method. Choose from: ", paste(valid_m, collapse = ", "))

  base_sc <- nira_result$baseline_scores
  int_sc  <- nira_result$intervention_scores
  vars    <- nira_result$kim_fit$vars
  n_vars  <- length(vars)

  if (verbose) {
    message("-- KIM NodeIdentifyR Algorithm (NIRA) Permutation Test --")
    message("  nodes = ", n_vars, "  |  nPerm = ", nPerm,
            "  |  correction = ", method)
    message("  paired = ", isTRUE(nira_result$paired))
  }

  set.seed(seed)

  sl <- data.frame(
    mean_other    = numeric(n_vars),
    sd_other      = numeric(n_vars),
    se_other      = numeric(n_vars),
    ci_other      = numeric(n_vars),
    ciLower_other = numeric(n_vars),
    ciUpper_other = numeric(n_vars),
    cohen_d       = numeric(n_vars),
    p             = numeric(n_vars),
    stringsAsFactors = FALSE
  )
  rownames(sl) <- vars

  last_pl <- NULL
  for (j in seq_len(n_vars)) {
    vname <- vars[j]
    if (verbose) message("  [", j, "/", n_vars, "] ", vname)

    pl <- kimPermutationTest(base_sc, int_sc[[vname]], nPerm = nPerm)

    sl$mean_other[j]    <- pl$mean_other
    sl$sd_other[j]      <- pl$sd_other
    sl$se_other[j]      <- pl$se_other
    sl$ci_other[j]      <- pl$ci_other
    sl$ciLower_other[j] <- pl$ciLower_other
    sl$ciUpper_other[j] <- pl$ciUpper_other
    sl$cohen_d[j]       <- pl$cohens_d
    sl$p[j]             <- pl$p_value

    last_pl <- pl
  }

  sl$p.adjust <- stats::p.adjust(sl$p, method = method)

  if (verbose)
    message("  Significant (p.adjust < 0.05): ",
            sum(sl$p.adjust < 0.05), " / ", n_vars)

  # -- plot_data --------------------------------------------------------------
  d1 <- data.frame(
    mean    = sl$mean_other,
    ciLower = sl$ciLower_other,
    ciUpper = sl$ciUpper_other,
    stringsAsFactors = FALSE
  )
  rownames(d1) <- vars

  d2 <- data.frame(
    mean    = last_pl$mean_original,
    ciLower = last_pl$ciLower_original,
    ciUpper = last_pl$ciUpper_original,
    stringsAsFactors = FALSE
  )
  rownames(d2) <- "baseline"

  d3      <- rbind(d2, d1)
  d3$id   <- 0:(nrow(d3) - 1L)
  d3$node <- rownames(d3)
  rownames(d3) <- NULL

  result <- list(
    stat        = sl,
    plot_data   = d3,
    nira_result = nira_result
  )
  class(result) <- "kim_perm_nira_result"
  result
}


#' @export
print.kim_perm_nira_result <- function(x, ...) {
  cat("\n-- KIM NodeIdentifyR Algorithm (NIRA) Permutation Test Result --\n\n")
  df <- x$stat
  df$sig <- ifelse(df$p.adjust < 0.001, "***",
                   ifelse(df$p.adjust < 0.01,  "**",
                          ifelse(df$p.adjust < 0.05,  "*", "")))
  out <- data.frame(
    node     = rownames(x$stat),
    mean     = round(df$mean_other, 4),
    cohen_d  = round(df$cohen_d, 4),
    p        = round(df$p, 4),
    p.adjust = round(df$p.adjust, 4),
    sig      = df$sig,
    stringsAsFactors = FALSE
  )
  out <- out[order(out$p.adjust), ]
  print(out, row.names = FALSE)
  cat("\n  Sig: '***' < 0.001, '**' < 0.01, '*' < 0.05\n\n")
  invisible(x)
}
