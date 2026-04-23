# =============================================================================
# simulationStudy.R  ——  KIM-NIRA 模拟研究（parallel）
#
# Study 2: Intervention Target Recovery
#
# 依赖：source("fitKIM.R")、source("runKIMNIRA.R")、parallel（内置）、glmnet、ggplot2
#
# 并行说明（Linux）：
#   全部使用 parallel::mclapply（fork 多进程，子进程自动继承父环境）。
#   n_cores 参数统一控制，默认 = detectCores() - 1。
#   Study 1/2：外层 36 个条件并行。
# =============================================================================


# ─────────────────────────────────────────────────────────────────────────────
# 内部工具：解析核数
# ─────────────────────────────────────────────────────────────────────────────

#' @keywords internal
.resolve_cores <- function(n_cores) {
  if (is.null(n_cores))
    n_cores <- max(1L, parallel::detectCores() - 1L)
  as.integer(n_cores)
}


# ─────────────────────────────────────────────────────────────────────────────
# 真实参数生成
# ─────────────────────────────────────────────────────────────────────────────

#' @keywords internal
generate_true_params <- function(p         = 6L,
                                 density   = c("sparse", "medium"),
                                 strength  = c("weak", "strong"),
                                 hub_boost = 2.0,
                                 seed      = 42L) {
  density  <- match.arg(density)
  strength <- match.arg(strength)
  set.seed(seed)

  vars           <- paste0("V", seq_len(p))
  target_density <- switch(density, sparse = 0.30, medium = 0.55)
  weight_range   <- switch(strength, weak = c(0.05, 0.20), strong = c(0.30, 0.80))

  J <- matrix(0, p, p, dimnames = list(vars, vars))
  diag(J) <- runif(p, weight_range[1], weight_range[2])

  off_idx    <- which(row(J) != col(J))
  n_off      <- max(round(length(off_idx) * target_density), p)
  sel        <- sample(off_idx, min(n_off, length(off_idx)))
  sgn        <- sample(c(1,-1), length(sel), replace=TRUE, prob=c(0.7,0.3))
  J[sel]     <- runif(length(sel), weight_range[1], weight_range[2]) * sgn

  hub      <- sample(p, 1)
  J_base   <- J
  for (i in seq_len(p)) {
    if (i == hub) next
    if (abs(J_base[i, hub]) < 1e-6)
      J_base[i, hub] <- runif(1, weight_range[1], weight_range[2]) *
        sample(c(1,-1), 1, prob=c(0.7,0.3))
  }

  boost <- hub_boost
  for (k in 1:20) {
    J_try <- J_base
    for (i in seq_len(p)) { if (i != hub) J_try[i, hub] <- J_base[i, hub] * boost }
    if (max(abs(eigen(J_try, only.values=TRUE)$values)) < 0.95) break
    boost <- boost * 0.85
  }
  J        <- J_try
  h        <- runif(p, -1.5, -0.3);  h[hub] <- runif(1, -0.8, -0.1)
  names(h) <- vars

  list(J=J, h=h, vars=vars, p=p, hub_node=hub, hub_var=vars[hub],
       density=density, strength=strength,
       true_rho = max(abs(eigen(J, only.values=TRUE)$values)),
       hub_boost_actual=boost, weight_range=weight_range)
}


#' @keywords internal
generate_kim_data <- function(true_params, n_persons=50L, n_time=60L,
                              h_person_sd=0.3, seed=123L) {
  J <- true_params$J;  h <- true_params$h
  vars <- true_params$vars;  p <- true_params$p
  set.seed(seed)

  all_data <- vector("list", n_persons)
  for (pid in seq_len(n_persons)) {
    h_pid <- h + rnorm(p, 0, h_person_sd)
    s     <- as.integer(runif(p) < kim_sigmoid(h_pid))
    states <- matrix(NA_integer_, nrow=n_time, ncol=p)
    states[1,] <- s
    for (t in 2:n_time) {
      eta <- as.numeric(s %*% t(J)) + h_pid
      s   <- as.integer(runif(p) < kim_sigmoid(eta))
      states[t,] <- s
    }
    df <- as.data.frame(states);  colnames(df) <- vars
    df$id <- pid;  df$time <- seq_len(n_time)
    all_data[[pid]] <- df
  }
  data <- do.call(rbind, all_data)
  data[, c("id","time",vars)]
}


# ─────────────────────────────────────────────────────────────────────────────
# 评估指标
# ─────────────────────────────────────────────────────────────────────────────

#' @keywords internal
compute_recovery_metrics <- function(J_true, J_est, h_true, h_est) {
  te <- abs(J_true)>1e-6;  ee <- abs(J_est)>1e-6
  tp <- sum(te&ee); fp <- sum(!te&ee); fn <- sum(te&!ee); tn <- sum(!te&!ee)
  sens <- if ((tp+fn)>0) tp/(tp+fn) else NA_real_
  spec <- if ((tn+fp)>0) tn/(tn+fp) else NA_real_
  prec <- if ((tp+fp)>0) tp/(tp+fp) else NA_real_
  f1   <- if (!is.na(prec)&&!is.na(sens)&&(prec+sens)>0)
            2*prec*sens/(prec+sens) else NA_real_
  sgn  <- if (sum(te)>0) mean(sign(J_est[te])==sign(J_true[te])) else NA_real_
  data.frame(
    J_RMSE=round(sqrt(mean((J_est-J_true)^2)),4),
    J_MAE =round(mean(abs(J_est-J_true)),4),
    J_cor =round(cor(as.vector(J_true),as.vector(J_est)),4),
    sensitivity=round(sens,4), specificity=round(spec,4),
    precision=round(prec,4), F1=round(f1,4),
    sign_accuracy=round(sgn,4),
    h_RMSE=round(sqrt(mean((h_est-h_true)^2)),4),
    h_MAE =round(mean(abs(h_est-h_true)),4),
    h_cor =round(cor(h_true,h_est),4),
    stringsAsFactors=FALSE)
}


#' 大样本单次仿真确定真实 NIRA 排名（ground truth）
#' @keywords internal
compute_true_nira_ranking <- function(true_params, n_sim=20000L, n_steps=30L,
                                      intervention_sd=2, type="alleviating",
                                      seed=999L) {
  J <- true_params$J;  h <- true_params$h
  vars <- true_params$vars;  p <- true_params$p

  delta       <- kim_compute_delta(h, intervention_sd, type)
  init_states <- kim_init_states(h, n_sim, vars, seed=seed)
  U_list      <- .kim_pregenerate_U(n_sim, p, n_steps, seed=seed+1000L)

  baseline <- kim_trajectory(init_states, J, h, n_steps,
                             target_node=NULL, delta=0, U_list=U_list)
  effects <- setNames(numeric(p), vars)
  for (j in seq_len(p)) {
    iv <- kim_trajectory(init_states, J, h, n_steps,
                         target_node=j, delta=delta, U_list=U_list)
    effects[j] <- mean(iv - baseline)
  }

  abs_eff <- abs(effects)
  ord     <- order(abs_eff, decreasing=TRUE)
  srt     <- sort(abs_eff, decreasing=TRUE)
  list(effects=effects, abs_effects=abs_eff,
       node_to_rank = rank(-abs_eff),
       ranked_vars  = vars[ord],
       gap_ratio    = if (length(srt)>=2 && srt[2]>1e-10)
                        (srt[1]-srt[2])/srt[2] else Inf)
}


# =============================================================================
# Study 2: Intervention Target Recovery
# =============================================================================

#' Study 2: 干预靶点恢复  |  外层 36 条件并行
#'
#' 每次 rep = 大样本真实排名（单次 n_sim_true）+ 单次 runKIMNIRA。
#' 无嵌套稳定性检验，无 rank_stability 过滤。
#'
#' @param p 节点数，默认 6L
#' @param n_reps 每条件重复，默认 50L
#' @param n_sim 估计排名仿真次数，默认 5000L
#' @param n_sim_true 真实排名仿真次数，默认 20000L
#' @param n_steps 轨迹步数，默认 30L
#' @param n_cores 并行核数，NULL = detectCores()-1
#' @param verbose 打印进度
#' @export
run_study2_target_recovery <- function(p=6L, n_reps=50L,
                                       n_sim=5000L, n_sim_true=20000L,
                                       n_steps=30L, n_cores=NULL, verbose=TRUE) {
  n_cores <- .resolve_cores(n_cores)

  conditions <- expand.grid(
    density=c("sparse","medium"), n_persons=c(50L,100L,200L),
    n_time=c(14L,30L,60L), strength=c("weak","strong"),
    stringsAsFactors=FALSE)
  n_cond <- nrow(conditions)

  if (verbose) {
    message("═══ Study 2 (v4): Intervention Target Recovery ═══")
    message("  p=",p,"  cores=",n_cores,"  conds=",n_cond,
            "  reps=",n_reps,"  total=",n_cond*n_reps)
    message("  n_sim=",n_sim,"  n_sim_true=",n_sim_true,"  n_steps=",n_steps)
  }

  cond_results <- parallel::mclapply(seq_len(n_cond), function(ic) {
    cond <- conditions[ic,]
    rep_list <- vector("list", n_reps)

    for (rep in seq_len(n_reps)) {
      s <- ic*10000L + rep

      tp <- tryCatch(generate_true_params(p=p, density=cond$density,
                       strength=cond$strength, hub_boost=1.0, seed=s),
                     error=function(e) NULL)
      if (is.null(tp)) next

      tn <- tryCatch(compute_true_nira_ranking(tp, n_sim=n_sim_true,
                       n_steps=n_steps, seed=s+9000L), error=function(e) NULL)
      if (is.null(tn)) next

      dat <- tryCatch(generate_kim_data(tp, n_persons=cond$n_persons,
                        n_time=cond$n_time, seed=s+5000L), error=function(e) NULL)
      if (is.null(dat)) next

      fit <- tryCatch(fitKIM(dat, vars=tp$vars, id="id", time="time",
                       verbose=FALSE), error=function(e) NULL)
      if (is.null(fit)) next

      en <- tryCatch(runKIMNIRA(fit, intervention_type="alleviating",
                       n_sim=n_sim, n_steps=n_steps,
                       seed=s+7000L, verbose=FALSE), error=function(e) NULL)
      if (is.null(en)) next

      vars    <- tp$vars
      tr      <- tn$ranked_vars
      er      <- en$ranked_targets
      k3      <- min(3L, p)
      enr     <- match(vars, er);  names(enr) <- vars
      tau     <- tryCatch(cor(tn$node_to_rank, enr, method="kendall"),
                          error=function(e) NA_real_)

      rep_list[[rep]] <- data.frame(
        density=cond$density, n_persons=cond$n_persons,
        n_time=cond$n_time, strength=cond$strength, rep=rep,
        top1_hit    =as.integer(er[1]==tr[1]),
        top3_overlap=round(sum(tr[1:k3] %in% er[1:k3])/k3, 4),
        kendall_tau =round(tau, 4),
        gap_ratio   =round(tn$gap_ratio, 4),
        stringsAsFactors=FALSE)
    }

    do.call(rbind, Filter(Negate(is.null), rep_list))
  }, mc.cores=n_cores)

  results <- do.call(rbind, Filter(Negate(is.null), cond_results))

  if (is.null(results) || nrow(results)==0) {
    warning("Study 2 无有效结果！")
    return(list(results=data.frame(), summary=NULL, conditions=conditions))
  }

  av  <- c("top1_hit","top3_overlap","kendall_tau","gap_ratio")
  grp <- list(density=results$density, n_persons=results$n_persons,
              n_time=results$n_time, strength=results$strength)
  summary_df <- aggregate(results[,av], by=grp, FUN=function(x) round(mean(x,na.rm=TRUE),3))
  summary_sd <- aggregate(results[,av], by=grp, FUN=function(x) round(sd(x,na.rm=TRUE),3))

  if (verbose) {
    message("═══ Study 2 完成  n=",nrow(results)," ═══")
    message("  Top-1 命中率：",round(mean(results$top1_hit,    na.rm=TRUE),3))
    message("  Top-3 重叠比：",round(mean(results$top3_overlap,na.rm=TRUE),3))
    message("  Kendall tau： ",round(mean(results$kendall_tau, na.rm=TRUE),3))
    print(summary_df)
  }
  list(results=results, summary=summary_df, summary_sd=summary_sd,
       conditions=conditions,
       params=list(p=p, n_reps=n_reps, n_sim=n_sim,
                   n_sim_true=n_sim_true, n_steps=n_steps))
}

