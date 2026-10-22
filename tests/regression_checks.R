# Regression checks that need no estimation packages (hand-built fits).
library(KimNIRA)
ns <- asNamespace("KimNIRA")
mkfit <- function(J, h, shifts = NULL) {
  vars <- colnames(J)
  f <- list(J = J, h = stats::setNames(h, vars), vars = vars, temporal_network = J,
            person_shifts = shifts, person_J = NULL)
  class(f) <- "dynaNIRA_kim"; f
}
set.seed(3); p <- 6; v <- paste0("S", 1:p)
J <- matrix(rnorm(p * p, 0, .4), p, p, dimnames = list(v, v)); diag(J) <- .8
sh <- matrix(rnorm(20 * p, 0, .5), 20, p, dimnames = list(NULL, v)); sh <- sweep(sh, 2, colMeans(sh))
f <- mkfit(J, rnorm(p, -1, .5), sh)

# 1. Exact and Monte Carlo engines estimate the same stationary contrast.
re <- runKIMNIRA(f, method = "exact", delta = -1, verbose = FALSE)
rm <- runKIMNIRA(f, method = "mc", delta = -1, n_sim = 4000, n_steps = 60, verbose = FALSE)
z <- (rm$effects - re$effects) / rm$mc_diagnostics$mc_se
stopifnot(all(abs(z) < 4))

# 2. Direct + spillover contributions add up to the total effect.
s <- re$summary
stopifnot(all(abs(s$direct_contribution + s$spillover_contribution - s$mean_diff) < 1e-12))

# 3. Stored MC trajectory reproduces the reported effects.
tr <- rm$trajectory; w <- (tr$burn_in + 2):(tr$n_steps + 1)
stopifnot(abs(mean(tr$interventions[[1]]$summary[w]) - mean(tr$baseline$summary[w]) -
              rm$effects[1]) < 1e-12)

# 4. Equal-probability calibration hits its target.
rp <- runKIMNIRA(f, method = "exact", dose_scale = "probability", target_dp = .05, verbose = FALSE)
m0 <- ns$.kim_pop_marginals(f)
chk <- vapply(seq_len(p), function(j) ns$.kim_pop_marginals(f, j, rp$delta[j])[j] - m0[j], 0)
stopifnot(all(abs(chk + .05) < 1e-6))

# 5. Finite-time response converges to the stationary effect.
resp <- ns$.kim_response(f, rp$delta, horizon = 200)$mean_changes
stopifnot(all(abs(resp[201, ] - rp$effects) < 1e-8))

# 6. Independent nodes give a unimodal count distribution; a strongly coupled
#    network is metastable, and exact first-passage times match simulation.
J0 <- matrix(0, p, p, dimnames = list(v, v))
stopifnot(!kimBistability(mkfit(J0, seq(-2, 1, length.out = p)), verbose = FALSE)$population_bimodal)
Jb <- matrix(0.9, p, p, dimnames = list(v, v)); diag(Jb) <- 1.5; Jb[, "S1"] <- 1.6; Jb["S1", "S1"] <- 2
fb <- mkfit(Jb, c(-3.2, rep(-3.4, 5)))
bb <- kimBistability(fb, delta = -1, verbose = FALSE)
stopifnot(bb$population_bimodal, nrow(bb$mean_field) == 2L)
tk <- ns$.kim_transition(fb$J, fb$h); pi0 <- ns$.kim_stationary(tk$P)$pi; K <- rowSums(tk$states)
cut <- bb$settings$cut; lowp <- pi0 * (K <= cut[1]); cs <- t(apply(tk$P, 1, cumsum))
set.seed(1); hits <- replicate(2000, { st <- sample.int(length(pi0), 1, prob = lowp); t <- 0
  repeat { st <- findInterval(runif(1), cs[st, ]) + 1; t <- t + 1; if (K[st] >= cut[2]) break }; t })
stopifnot(abs(mean(hits) - bb$basins$mfpt_low_to_high) < 4 * sd(hits) / sqrt(2000))
cat("All regression checks passed.\n")
