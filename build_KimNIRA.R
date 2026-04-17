#!/usr/bin/env Rscript
# =============================================================================
# build_KimNIRA.R
# Script to build, check, and install the KimNIRA package
#
# Usage (from the directory containing the KimNIRA/ folder):
#   Rscript build_KimNIRA.R
#
# Or from within R:
#   source("build_KimNIRA.R")
# =============================================================================

# ── 0. Verify required build tools ──────────────────────────────────────────

required_build_pkgs <- c("devtools", "roxygen2", "usethis")
missing_build <- required_build_pkgs[
  !sapply(required_build_pkgs, requireNamespace, quietly = TRUE)
]

if (length(missing_build) > 0) {
  message("Installing missing build tools: ", paste(missing_build, collapse = ", "))
  install.packages(missing_build, repos = "https://cloud.r-project.org")
}

library(devtools)
library(roxygen2)

# ── 1. Set working directory to the package root ────────────────────────────
# Adjust this path to where your KimNIRA/ folder lives.
pkg_path <- "D:/Rdaima/dynaNIRA/dytemp_nira_1.0/KimNIRA_package/KimNIRA"

if (!dir.exists(pkg_path))
  stop("Package directory not found: ", pkg_path,
       "\nMake sure KimNIRA/ is in your working directory.")

# ── 2. Install runtime dependencies ─────────────────────────────────────────

runtime_deps <- c("glmnet", "ggplot2", "tidyr", "ggsci")
optional_deps <- c("qgraph", "Matrix")

missing_runtime <- runtime_deps[
  !sapply(runtime_deps, requireNamespace, quietly = TRUE)
]
if (length(missing_runtime) > 0) {
  message("Installing runtime dependencies: ", paste(missing_runtime, collapse = ", "))
  install.packages(missing_runtime, repos = "https://cloud.r-project.org")
}

missing_optional <- optional_deps[
  !sapply(optional_deps, requireNamespace, quietly = TRUE)
]
if (length(missing_optional) > 0) {
  message("Installing optional dependencies: ", paste(missing_optional, collapse = ", "))
  install.packages(missing_optional, repos = "https://cloud.r-project.org")
}

# ── 3. Generate documentation from roxygen2 tags ────────────────────────────

message("\n[Step 1/4] Generating documentation with roxygen2...")

# Clean up any functions sourced directly into the global environment
# that would conflict with the package namespace during roxygenise().
kim_fns <- c("fitKIM", "runKIMNIRA", "kim_sigmoid", "kim_compute_delta",
             "kim_delta_prob_table", "kim_init_states", "kim_step",
             "kim_trajectory", "kimCentrality", "kimPermutationTest",
             "kimPermutationNIRAtest", "kimStabilityNIRAtest",
             "kimModerationTest", "kimPlotNetwork", "kimPlotAllNetworks",
             "kimPlotNIRA", "kimPlotStability", "kimPlotTrajectory")
to_remove <- kim_fns[sapply(kim_fns, exists, envir = .GlobalEnv,
                            inherits = FALSE)]
if (length(to_remove) > 0L) {
  message("  Removing ", length(to_remove),
          " global object(s) that would conflict with package namespace: ",
          paste(to_remove, collapse = ", "))
  rm(list = to_remove, envir = .GlobalEnv)
}

roxygen2::roxygenise(pkg_path)
message("  Documentation generated.")

# ── 4. Check the package ─────────────────────────────────────────────────────

message("\n[Step 2/4] Running R CMD check...")
check_result <- devtools::check(
  pkg        = pkg_path,
  cran       = FALSE,
  error_on   = "warning",  # Change to "error" for stricter checks
  quiet      = FALSE
)

# ── 5. Build the source tarball ──────────────────────────────────────────────

message("\n[Step 3/4] Building source package tarball...")
tarball <- devtools::build(
  pkg       = pkg_path,
  path      = "D:/Rdaima/dynaNIRA/dytemp_nira_1.0",
  binary    = FALSE,
  vignettes = FALSE,
  quiet     = FALSE
)
message("  Built: ", tarball)

# ── 6. Install the package ───────────────────────────────────────────────────

message("\n[Step 4/4] Installing KimNIRA...")
devtools::install(
  pkg             = pkg_path,
  reload          = TRUE,
  quick           = FALSE,
  build_vignettes = FALSE,
  quiet           = FALSE
)
message("  KimNIRA installed successfully.\n")

# ── 7. Quick smoke test ──────────────────────────────────────────────────────

message("Running quick smoke test...")
library(KimNIRA)

# Verify all exported functions are accessible
exported_fns <- c(
  "fitKIM", "runKIMNIRA",
  "kim_sigmoid", "kim_compute_delta", "kim_delta_prob_table",
  "kim_init_states", "kim_step", "kim_trajectory",
  "kimCentrality",
  "kimModerationTest",
  "kimPermutationTest", "kimPermutationNIRAtest",
  "kimStabilityNIRAtest",
  "kimPlotNetwork", "kimPlotAllNetworks",
  "kimPlotNIRA", "kimPlotStability", "kimPlotTrajectory"
)

ok <- sapply(exported_fns, existsMethod, signature = character(0),
             USE.NAMES = TRUE)
# simpler check:
ok <- sapply(exported_fns,
             function(f) exists(f, envir = asNamespace("KimNIRA")),
             USE.NAMES = TRUE)

if (all(ok)) {
  message("  All ", length(exported_fns), " exported functions found. PASS\n")
} else {
  message("  WARNING: Missing functions: ",
          paste(names(ok)[!ok], collapse = ", "), "\n")
}

message("=================================================")
message("  KimNIRA package build complete!")
message("  Tarball: ", tarball)
message("  Load with: library(KimNIRA)")
message("=================================================\n")


install.packages("D:/Rdaima/dynaNIRA/dytemp_nira_1.0/KimNIRA_1.1.0.tar.gz", repos = NULL, type = "source")
