source("~/Desktop/Rdaima/dynaNIRA/dytemp_nira_1.0/fitKIM.R")
source("~/Desktop/Rdaima/dynaNIRA/dytemp_nira_1.0/runKIMNIRA.R")
source("~/Desktop/Rdaima/dynaNIRA/dytemp_nira_1.0/simulationStudy_new.R")

s2 <- run_study2_target_recovery(p=4L, n_reps=50L, n_sim=5000L, n_sim_true=5000L,
                                 n_cores=20, verbose=TRUE)
s2
s2_10 <- run_study2_target_recovery(p=10L, n_reps=50L, n_sim=5000L, n_sim_true=5000L,
                                    n_cores=20, verbose=TRUE)
s2_10
s2_20 <- run_study2_target_recovery(p=20L, n_reps=50L, n_sim=5000L, n_sim_true=5000L,
                                    n_cores=20, verbose=TRUE)
s2_20


source("~/Desktop/Rdaima/dynaNIRA/dytemp_nira_1.0/plot_study2.R")
plot_study2_results(s2,s2_10,s2_20)



