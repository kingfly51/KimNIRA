# KimNIRA 0.3.1 (revision)

## New
* `kimBistability()`, `kimPlotLandscape()`: exact stationary distribution of the
  number of active nodes, potential U(K) = -log P(K), modes and barriers, basin
  masses, exact mean first-passage times between low/high basins, within-person
  vs between-person origin of bimodality, mean-field fixed points, and how each
  intervention reshapes the landscape.
* `kimMixing()`: exact relaxation time -1/log|lambda_2| per person.
* `kimDoseResponse()`, `kimPlotDoseResponse()`: ranking over a logit dose grid
  and under equal-probability doses; `runKIMNIRA(dose_scale = "probability",
  target_dp = )` and per-target `delta` vectors.
* `kimPathwayDecomposition()`: knock-out contrasts (outgoing pathways, feedback
  loops, self-persistence), direct/spillover split, half-effect time, and rank
  agreement with out-strength / out-EI.
* `kimContemporaneousCheck()`: test of the conditional-independence assumption.
* `kimChanceLevel()`: chance levels for target-recovery metrics.

## Changed / fixed
* fitKIM: `between_person = "fe"` (Mundlak within/between decomposition) removed,
  together with `J_between` and its cross-validation helper: its autoregressive
  terms carry an O(1/T) Nickell-type bias that is large at EMA lengths. The
  options are now `"re"` (default), `"pooled"` and `"res"`; cross-validation and
  glmnet are only used by `"pooled"`. `cv_gap` raised to at least `lags`;
  slope-structure rule `< 6` nodes (mlVAR default); constant-outcome
  message for re/res corrected; heuristic relaxation now at the mean-field fixed
  point and maximised over persons; `transition_types` reports what one step is.
* Exact engine: soft size limit `options(KimNIRA.max_exact_p = 12)`; corrected
  statements on uniqueness (P > 0) and on numerical accuracy; optional rcond.
* runKIMNIRA: R-hat pilot chains share the same persons across starts; exact
  result reports `paired = NA`; "best target" wording replaced by "top-ranked".
* kimPlotTrajectory: no silent switch from MC to exact; step axis labelled as
  observation intervals.
* kimPlotStability: sequential palette; x axis states person-bootstrap with
  re-estimation. kimMCError: +/- 2 SD labelled as single-run range, not a CI.
* kimStabilityBootstrap: manifest includes `person_J`; optional per-resample
  recalibration of equal-probability doses.
* kimCentrality: reports the autoregressive weight next to the indices.

# KimNIRA 0.4.0

## Estimation is restructured around two person-level models plus edge screening

* `fitKIM()` now offers only `between_person = "re"` (random intercept) and
  `"res"` (random intercept and random slopes). The former `"pooled"` option,
  the glmnet L1 path, `.kim_block_cv()` and the associated `lambda_rule`,
  `n_folds`, `cv_scheme`, `cv_gap`, `nlambda` and `penalize_ar` arguments have
  been removed. `pooled` remains available as `screen = "none"` only in the
  sense that no cell is screened out; it is not a separate estimator.

* `fitKIM()` gained `screen = c("fdr", "significance", "none")`. The default,
  `"fdr"`, keeps cells whose unpenalised multilevel coefficient is significant
  under Benjamini-Hochberg control at `fdr_q` over all p^2 candidate cells, and
  then refits the surviving cells WITHOUT a penalty. This two-stage structure
  removes the L1 shrinkage from the reported weights.

  Why: on a simulated 6-node network, an information-criterion penalty path
  (rpql/glmnet style) put its minimum at the least-penalised end of the path and
  gave specificity 0.36-0.50, because the likelihood gap between the sparse and
  dense model is far larger than log(n) per parameter. The path itself contained
  a lambda with specificity and sensitivity both 1.00. BH-FDR screening reached
  specificity 0.96 with sensitivity 0.94, and matched Bonferroni wherever edges
  were detectable while keeping usable power when they were not (F1 0.79 versus
  0.54 at |J| in 0.25-0.45).

* The fit now carries `support`, `stage1_beta`, `stage1_se`, `screen` and
  `screen_threshold` so the screening step is inspectable.

* New exported helper `kimSelectionStats(fit, truth)` reports TP/FP/FN/TN,
  sensitivity, specificity, precision, F1 and sign accuracy against a known
  transition matrix.

* `kimMarkovCheck()` no longer has a penalised branch; its `inner_folds` and
  `lambda_grid` arguments are gone. It works for `re` and `res` fits only.

* `glmnet` is no longer an imported dependency.

## Upgrading

`fitKIM(..., between_person = "pooled")` errors. Use `between_person = "re"`
(or `"res"`), and `screen = "none"` if you want a dense J for comparison.
