#' KimNIRA: Kinetic Ising Model with NodeIdentifyR Algorithm for Dynamic
#' Network Intervention Analysis
#'
#' @description
#' \pkg{KimNIRA} integrates the \strong{Kinetic Ising Model (KIM)} with the
#' \strong{NodeIdentifyR Algorithm (NIRA)} to provide a complete six-stage
#' workflow for dynamic binary psychological network analysis on ecological
#' momentary assessment (EMA) data.
#'
#' ---
#'
#' \strong{Stage 1 -- Network Estimation}
#' \itemize{
#'   \item \code{\link{fitKIM}} -- Fits a KIM to long-format binary EMA data.
#'     Estimates the \emph{temporal} (within-person) network via
#'     LASSO-penalised pseudo-likelihood logistic regression, and the
#'     \emph{between-subject} network via EBICglasso on person-level node
#'     means. Returns edge weight matrix \eqn{J}, threshold vector \eqn{h},
#'     weighted Jacobian spectral radius \eqn{\rho^*} (the correct KIM
#'     stability criterion), and diagnostic metadata.
#' }
#'
#' \strong{Stage 2 -- Network Centrality}
#' \itemize{
#'   \item \code{\link{kimCentrality}} -- Computes and plots centrality indices
#'     for both the temporal and between-subject networks. For the directed
#'     temporal network, reports InStrength, OutStrength, and signed Expected
#'     Influence (in/out), correctly accounting for the fitKIM direction
#'     convention (\code{J[i,j]} = j drives i). For the undirected
#'     between-subject network, reports Strength, Closeness, Betweenness, and
#'     Expected Influence. Optionally returns custom \pkg{ggplot2} panels.
#' }
#'
#' \strong{Stage 3 -- Moderation Testing (NodeIdentifyR Algorithm Prerequisite)}
#' \itemize{
#'   \item \code{\link{kimModerationTest}} -- Tests whether any node moderates
#'     temporal edge weights before running the NodeIdentifyR Algorithm
#'     (NIRA). The NodeIdentifyR Algorithm (NIRA) assumes that perturbing
#'     one node's threshold leaves all \eqn{J_{ij}} unchanged; this function
#'     verifies that assumption via likelihood-ratio tests and
#'     permutation-based inference on three-way interaction terms
#'     \eqn{M_{ijk} \cdot s_j(t) \cdot s_k(t)}. An empty
#'     \code{significant_moderators} result confirms the NodeIdentifyR
#'     Algorithm (NIRA) is appropriate.
#' }
#'
#' \strong{Stage 4 -- NodeIdentifyR Algorithm (NIRA) Intervention Simulation}
#' \itemize{
#'   \item \code{\link{runKIMNIRA}} -- Runs the NodeIdentifyR Algorithm (NIRA)
#'     via paired Monte Carlo simulation. For each node, an alleviating or
#'     aggravating threshold shift \eqn{\delta} is applied and the resulting
#'     change in projected network activation rate is recorded. Nodes are
#'     ranked by absolute effect size \eqn{|\Delta\bar{s}|}.
#'   \item \code{\link{kim_compute_delta}} -- Computes the intervention
#'     threshold shift \eqn{\delta = \pm k \cdot \mathrm{SD}(h)}.
#'   \item \code{\link{kim_delta_prob_table}} -- Tabulates per-node direct
#'     marginal probability changes from the threshold shift.
#'   \item \code{\link{kim_sigmoid}}, \code{\link{kim_init_states}},
#'     \code{\link{kim_step}}, \code{\link{kim_trajectory}} -- Low-level
#'     simulation primitives (sigmoid activation, initial state sampling,
#'     single-step update, full trajectory aggregation).
#' }
#'
#' \strong{Stage 5 -- Significance Testing}
#' \itemize{
#'   \item \code{\link{kimPermutationTest}} -- Paired sign-flip permutation
#'     test for a single intervention node; returns two-sided permutation
#'     p-value and Cohen's \eqn{d}.
#'   \item \code{\link{kimPermutationNIRAtest}} -- Applies
#'     \code{kimPermutationTest} to every node in a
#'     \code{\link{runKIMNIRA}} result, corrects for multiple comparisons,
#'     and produces a plotting-ready data frame for
#'     \code{\link{kimPlotNIRA}}.
#' }
#'
#' \strong{Stage 6 -- Stability Assessment}
#' \itemize{
#'   \item \code{\link{kimStabilityNIRAtest}} -- Repeats the NodeIdentifyR
#'     Algorithm (NIRA) across \eqn{n_{\text{runs}}} independent random
#'     seeds and tabulates, for each node, the proportion of runs in which
#'     it occupied each rank position. Derives a position-weighted consensus
#'     ranking (\code{stable_ranking}).
#' }
#'
#' \strong{Visualisation}
#' \itemize{
#'   \item \code{\link{kimPlotNetwork}} -- Plots a single KIM network
#'     (temporal or between-subject) via \pkg{qgraph}.
#'   \item \code{\link{kimPlotAllNetworks}} -- Side-by-side temporal and
#'     between-subject network panels with shared node layout.
#'   \item \code{\link{kimCentrality}} -- Centrality plots via
#'     \code{qgraph::centralityPlot}, with an optional \pkg{ggplot2} mode
#'     (\code{use_ggplot = TRUE}).
#'   \item \code{\link{kimPlotNIRA}} -- Mean projected activation rate plot
#'     with 95\% CI error bars and corrected significance stars; mirrors the
#'     style of Lunansky et al. (2022) Figure 3.
#'   \item \code{\link{kimPlotStability}} -- Horizontal stacked bar chart of
#'     rank-proportion stability across NodeIdentifyR Algorithm (NIRA) runs.
#'   \item \code{\link{kimPlotTrajectory}} -- Longitudinal intervention
#'     trajectory plots: \code{"summary"} (group mean), \code{"node"}
#'     (per-node facets), \code{"diff"} (difference curves),
#'     \code{"panel"} (top-\eqn{n} target comparison), and
#'     \code{"spaghetti"} (individual trajectories).
#'   \item \code{\link{plot.kim_moderation_result}} -- Heatmap of moderation
#'     test results, colour-coded by interaction coefficient or
#'     \eqn{-\log_{10}(p_{\text{adj}})}.
#' }
#'
#' ---
#'
#' @section Recommended workflow:
#' \preformatted{
#' # 1. Fit the KIM
#' fit <- fitKIM(data, vars = c("dep", "anx", "sleep"),
#'               id = "id", time = "time")
#'
#' # 2. Inspect centrality of temporal and between-subject networks
#' cent <- kimCentrality(fit, use_ggplot = TRUE)
#' cent$temporal   # OutStrength identifies the strongest drivers
#' cent$between    # Strength in the between-subject network
#'
#' # 3. Verify NodeIdentifyR Algorithm (NIRA) prerequisite
#' mod <- kimModerationTest(data, fit, n_perm = 500)
#' print(mod)           # empty significant_moderators => proceed
#' plot(mod)            # moderation heatmap
#'
#' # 4. Run the NodeIdentifyR Algorithm (NIRA)
#' nira <- runKIMNIRA(fit, intervention_type = "alleviating")
#' print(nira)
#'
#' # 5. Significance testing
#' perm <- kimPermutationNIRAtest(nira)
#' kimPlotNIRA(perm, "alleviating")
#'
#' # 6. Stability assessment
#' stab <- kimStabilityNIRAtest(fit, n_runs = 20)
#' kimPlotStability(stab)
#' }
#'
#' @references
#' Besag, J. (1975). Statistical analysis of non-lattice data.
#'   \emph{The Statistician}, 24(3), 179-195.
#'
#' Campajola, C., Lillo, F., Tantari, D., & Viola, M. (2022).
#'   Unveiling the relation between herding and liquidity with trader
#'   lead-lag networks. \emph{Scientific Reports}, 12, 19339.
#'
#' Epskamp, S., & Fried, E. I. (2018). A tutorial on regularized partial
#'   correlation networks. \emph{Psychological Methods}, 23(4), 617-634.
#'
#' Lunansky, G., van Borkulo, C. D., Blanken, T. F., Cramer, A. O. J., &
#'   Borsboom, D. (2022). Intervening on psychopathology networks:
#'   Evaluating intervention targets through simulations.
#'   \emph{Methods}, 204, 29-37.
#'
#' @author Fei Wang \email{bjnwangfei0501@@outlook.com}
#'
#' @seealso
#' \itemize{
#'   \item Package repository: \url{https://github.com/kingfly51/KimNIRA}
#'   \item Bug reports: \url{https://github.com/kingfly51/KimNIRA/issues}
#' }
#'
"_PACKAGE"
