# =============================================================================
# simulationStudy.R  ——  KIM-NIRA 模拟研究（修订版 v2）
#
# 本文件实现论文 Simulation Study 部分的全部三项研究：
#   Study 1: Parameter Recovery（参数恢复）
#
# 依赖：
#   source("fitKIM.R")
#   source("runKIMNIRA.R")
#   glmnet, ggplot2
#
# =============================================================================


# ─────────────────────────────────────────────────────────────────────────────
# 辅助函数：真实网络生成
# ─────────────────────────────────────────────────────────────────────────────

#' 生成真实的 KIM 参数（带枢纽节点设计）
#'
#' @description
#' 刻意构造一个"枢纽节点"（hub），使其对网络中其他节点的传出影响
#' 最强，确保 NIRA 排名有明确分离。具体策略：
#'   - 随机选一个节点作为 hub，其传出边权乘以 hub_boost
#'   - 其余节点按常规密度和强度生成
#'   - 截距 h 有层次分化（hub 节点更居中，更容易受干预影响）
#'
#' 边权范围修正为逻辑回归系数尺度：
#'   weak:   J_ij ∈ [0.05, 0.20]  →  ΔP ≈ 0.01–0.05（概率尺度）
#'   strong: J_ij ∈ [0.30, 0.80]  →  ΔP ≈ 0.07–0.20（概率尺度）
#'
#' 谱半径控制：自适应缩小 hub_boost，直到加权 Jacobian 谱半径
#'   rho_star = rho(diag(sigma'(h)) J) < 0.95。
#'   注意：正确的 KIM 稳定性准则是 rho_star 而非 rho(J)；
#'   由于 sigma'(h_i) <= 0.25，rho_star <= 0.25 * rho(J)，
#'   用 rho(J) 控制会过度压缩 hub_boost。
#'
#' @param p 节点数
#' @param density 网络密度
#' @param strength 边权强度
#' @param hub_boost hub 节点传出边的放大倍数，默认 2.0
#' @param seed 随机种子
#' @keywords internal
generate_true_params <- function(p           = 6L,
                                 density     = c("sparse", "medium"),
                                 strength    = c("weak", "strong"),
                                 hub_boost   = 2.0,
                                 seed        = 42L) {
  density  <- match.arg(density)
  strength <- match.arg(strength)
  set.seed(seed)
  
  vars <- paste0("V", seq_len(p))
  
  target_density <- switch(density,
                           sparse = 0.30,
                           medium = 0.55)
  
  # [v2-1] 修正边权范围：逻辑回归系数尺度
  weight_range <- switch(strength,
                         weak   = c(0.05, 0.20),
                         strong = c(0.30, 0.80))
  
  J <- matrix(0, p, p, dimnames = list(vars, vars))
  
  # 自回归边（对角线）：正值，带分化
  diag(J) <- runif(p, weight_range[1], weight_range[2])
  
  # 非对角线边：按密度随机选择
  off_diag_idx <- which(row(J) != col(J))
  n_off_diag <- round(length(off_diag_idx) * target_density)
  n_off_diag <- max(n_off_diag, p)  # 至少 p 条边
  selected <- sample(off_diag_idx, min(n_off_diag, length(off_diag_idx)))
  
  signs <- sample(c(1, -1), length(selected), replace = TRUE, prob = c(0.7, 0.3))
  weights <- runif(length(selected), weight_range[1], weight_range[2]) * signs
  J[selected] <- weights
  
  # 枢纽节点：随机选一个，放大其传出边
  hub_node <- sample(p, 1)
  
  # 保存 hub 放大前的 J（用于谱半径控制时回退）
  J_base <- J
  
  # 确保 hub 对所有其他节点都有非零传出边
  for (i in seq_len(p)) {
    if (i == hub_node) next
    if (abs(J_base[i, hub_node]) < 1e-6) {
      J_base[i, hub_node] <- runif(1, weight_range[1], weight_range[2]) *
        sample(c(1, -1), 1, prob = c(0.7, 0.3))
    }
  }
  
  # [v2-2] 自适应 hub_boost：逐步缩小直到 rho_star < rho_target
  #
  # 修正说明：
  #   原代码用 rho(J)（J 矩阵自身的谱半径）作为稳定性准则，这是错误的。
  #   KIM 的 one-step 更新为 E[s(t+1)|s(t)] = sigma(J s(t) + h)，
  #   其局部线性化 Jacobian 为 diag(sigma'(h)) %*% J，其中
  #     sigma'(h_i) = sigma(h_i) * (1 - sigma(h_i))  in (0, 0.25]
  #   正确的稳定性准则是：
  #     rho_star = max|eigenvalue(diag(sigma'(h)) %*% J)| < 1
  #   由于 sigma'(h_i) <= 0.25，rho_star <= 0.25 * rho(J)，
  #   用 rho(J) 做控制会严重高估不稳定性，导致 hub_boost 被过度压缩。
  #
  #   注意：此时 h 尚未生成，用分布中值 h_tmp = -0.9 近似计算 rho_star，
  #   待 h 生成后再重新计算精确的 final_rho_star。
  current_boost <- hub_boost
  rho_target    <- 0.95    # rho_star 的目标上界

  h_tmp <- rep(-0.9, p)    # 近似截距（用于稳定性控制）
  s_tmp <- 1 / (1 + exp(-h_tmp))
  d_tmp <- s_tmp * (1 - s_tmp)   # sigma'(h_tmp), in (0, 0.25]

  for (attempt in 1:20) {
    J_try <- J_base
    for (i in seq_len(p)) {
      if (i == hub_node) next
      J_try[i, hub_node] <- J_base[i, hub_node] * current_boost
    }

    # rho_star：加权 Jacobian 谱半径（正确的 KIM 稳定性准则）
    rho_star_try <- max(abs(eigen(diag(d_tmp) %*% J_try,
                                  only.values = TRUE)$values))
    if (rho_star_try < rho_target) break

    # 缩小 boost（每次乘 0.85）
    current_boost <- current_boost * 0.85
  }

  J <- J_try
  actual_boost <- current_boost

  # 截距向量：hub 节点更居中（更容易被干预"推动"）
  h <- runif(p, -1.5, -0.3)
  h[hub_node] <- runif(1, -0.8, -0.1)  # hub 更接近 0
  names(h) <- vars

  # 用真实 h 精确计算两个谱半径并记录
  s_bar          <- 1 / (1 + exp(-h))
  d_bar          <- s_bar * (1 - s_bar)          # sigma'(h), in (0, 0.25]
  final_rho_J    <- max(abs(eigen(J,                    only.values = TRUE)$values))
  final_rho_star <- max(abs(eigen(diag(d_bar) %*% J,   only.values = TRUE)$values))

  list(J = J, h = h, vars = vars, p = p,
       hub_node = hub_node, hub_var = vars[hub_node],
       density = density, strength = strength,
       true_rho       = final_rho_star,   # 主要稳定性指标（rho_star）
       true_rho_J     = final_rho_J,      # 保留原始 rho(J) 供参考
       true_rho_star  = final_rho_star,   # 显式命名，与 fitKIM 输出一致
       hub_boost_actual = actual_boost,
       weight_range = weight_range)
}


#' 基于 KIM 参数生成模拟 EMA 数据
#'
#' @param true_params generate_true_params() 的返回值
#' @param n_persons 被试数量
#' @param n_time 每人时间点数
#' @param h_person_sd 个体间截距偏移的标准差，默认 0.3
#' @param seed 随机种子
#' @keywords internal
generate_kim_data <- function(true_params,
                              n_persons   = 50L,
                              n_time      = 60L,
                              h_person_sd = 0.3,
                              seed        = 123L) {
  J    <- true_params$J
  h    <- true_params$h
  vars <- true_params$vars
  p    <- true_params$p
  
  set.seed(seed)
  
  all_data <- vector("list", n_persons)
  
  for (pid in seq_len(n_persons)) {
    # 每个被试有自己的截距偏移
    h_pid <- h + rnorm(p, 0, h_person_sd)
    
    # 初始状态
    s <- as.integer(runif(p) < kim_sigmoid(h_pid))
    
    states <- matrix(NA_integer_, nrow = n_time, ncol = p)
    states[1, ] <- s
    
    for (t in 2:n_time) {
      eta <- as.numeric(s %*% t(J)) + h_pid
      probs <- kim_sigmoid(eta)
      s <- as.integer(runif(p) < probs)
      states[t, ] <- s
    }
    
    df <- as.data.frame(states)
    colnames(df) <- vars
    df$id   <- pid
    df$time <- seq_len(n_time)
    all_data[[pid]] <- df
  }
  
  data <- do.call(rbind, all_data)
  data <- data[, c("id", "time", vars)]
  data
}


# ─────────────────────────────────────────────────────────────────────────────
# 辅助函数：评估指标
# ─────────────────────────────────────────────────────────────────────────────

#' 计算参数恢复指标
#' @keywords internal
compute_recovery_metrics <- function(J_true, J_est, h_true, h_est) {
  
  p <- nrow(J_true)
  
  j_rmse <- sqrt(mean((J_est - J_true)^2))
  j_mae  <- mean(abs(J_est - J_true))
  j_cor  <- cor(as.vector(J_true), as.vector(J_est))
  
  true_edges <- abs(J_true) > 1e-6
  est_edges  <- abs(J_est) > 1e-6
  
  tp <- sum(true_edges & est_edges)
  fp <- sum(!true_edges & est_edges)
  fn <- sum(true_edges & !est_edges)
  tn <- sum(!true_edges & !est_edges)
  
  sensitivity <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
  specificity <- if ((tn + fp) > 0) tn / (tn + fp) else NA_real_
  precision   <- if ((tp + fp) > 0) tp / (tp + fp) else NA_real_
  f1          <- if (!is.na(precision) && !is.na(sensitivity) &&
                     (precision + sensitivity) > 0)
    2 * precision * sensitivity / (precision + sensitivity) else NA_real_
  
  if (sum(true_edges) > 0) {
    sign_agree <- mean(sign(J_est[true_edges]) == sign(J_true[true_edges]))
  } else {
    sign_agree <- NA_real_
  }
  
  h_rmse <- sqrt(mean((h_est - h_true)^2))
  h_mae  <- mean(abs(h_est - h_true))
  h_cor  <- cor(h_true, h_est)
  
  data.frame(
    J_RMSE        = round(j_rmse, 4),
    J_MAE         = round(j_mae, 4),
    J_cor         = round(j_cor, 4),
    sensitivity   = round(sensitivity, 4),
    specificity   = round(specificity, 4),
    precision     = round(precision, 4),
    F1            = round(f1, 4),
    sign_accuracy = round(sign_agree, 4),
    h_RMSE        = round(h_rmse, 4),
    h_MAE         = round(h_mae, 4),
    h_cor         = round(h_cor, 4),
    stringsAsFactors = FALSE
  )
}


#' 计算 KIM-NIRA 的"真实"干预效果排名
#'
#' @description
#' 使用真实参数和大样本多次仿真获得 ground truth 排名。
#' 同时输出效应分离度（gap ratio），便于判断排名是否可辨。
#'
#' @keywords internal
compute_true_nira_ranking <- function(true_params,
                                      n_sim            = 5000L,
                                      n_steps          = 30L,
                                      intervention_sd  = 2,
                                      type             = "alleviating",
                                      n_replicates     = 3L,
                                      seed             = 999L) {
  J    <- true_params$J
  h    <- true_params$h
  vars <- true_params$vars
  p    <- true_params$p
  
  delta <- kim_compute_delta(h, intervention_sd, type)
  
  # 多次独立仿真取平均，确保排名稳定
  all_effects <- matrix(NA_real_, nrow = n_replicates, ncol = p)
  colnames(all_effects) <- vars
  
  for (r in seq_len(n_replicates)) {
    r_seed <- seed + (r - 1L) * 77777L
    set.seed(r_seed)
    init_states <- kim_init_states(h, n_sim, vars, seed = r_seed)
    U_list <- .kim_pregenerate_U(n_sim, p, n_steps, seed = r_seed + 1000L)
    
    baseline <- kim_trajectory(init_states, J, h, n_steps,
                               target_node = NULL, delta = 0,
                               U_list = U_list)
    
    for (j in seq_len(p)) {
      intervention <- kim_trajectory(init_states, J, h, n_steps,
                                     target_node = j, delta = delta,
                                     U_list = U_list)
      all_effects[r, j] <- mean(intervention - baseline)
    }
  }
  
  # 取各次仿真的平均效应
  avg_effects <- colMeans(all_effects)
  abs_eff <- abs(avg_effects)
  
  # 排名（节点→排名映射）
  ranking_order <- order(abs_eff, decreasing = TRUE)  # 名次→节点索引
  node_to_rank  <- rank(-abs_eff)                      # 节点→排名（1=最强）
  
  # 效应分离度：前两名效应差 / 第二名效应
  sorted_abs <- sort(abs_eff, decreasing = TRUE)
  if (length(sorted_abs) >= 2 && sorted_abs[2] > 1e-10) {
    gap_ratio <- (sorted_abs[1] - sorted_abs[2]) / sorted_abs[2]
  } else {
    gap_ratio <- Inf
  }
  
  # 排名稳定性：各次仿真的排名是否一致
  rank_stability <- mean(apply(all_effects, 1, function(e) {
    cor(rank(-abs(e)), node_to_rank, method = "kendall")
  }))
  
  list(
    effects        = avg_effects,
    abs_effects    = abs_eff,
    node_to_rank   = node_to_rank,         # V1排第几, V2排第几, ...
    ranked_vars    = vars[ranking_order],   # 排名第1的是谁, 第2的是谁, ...
    gap_ratio      = gap_ratio,
    rank_stability = rank_stability,
    all_effects    = all_effects
  )
}


# =============================================================================
# Study 1: Parameter Recovery
# =============================================================================

#' Study 1: 参数恢复模拟研究
#' @export
run_study1_parameter_recovery <- function(p       = 6L,
                                          n_reps  = 50L,
                                          verbose = TRUE) {
  
  conditions <- expand.grid(
    density   = c("sparse", "medium"),
    n_persons = c(50L, 100L, 200L),
    n_time    = c(14L, 30L, 60L),
    strength  = c("weak", "strong"),
    stringsAsFactors = FALSE
  )
  n_cond <- nrow(conditions)
  
  if (verbose) {
    message("═══ Study 1: Parameter Recovery ═══")
    message("  节点数 p = ", p)
    message("  条件数 = ", n_cond)
    message("  每条件重复 = ", n_reps)
    message("  总运行次数 = ", n_cond * n_reps)
  }
  
  all_results <- vector("list", n_cond * n_reps)
  idx <- 0L
  
  for (ic in seq_len(n_cond)) {
    cond <- conditions[ic, ]
    
    if (verbose)
      message("\n  条件 ", ic, "/", n_cond,
              ": density=", cond$density,
              " n_persons=", cond$n_persons,
              " n_time=", cond$n_time,
              " strength=", cond$strength)
    
    for (rep in seq_len(n_reps)) {
      idx <- idx + 1L
      rep_seed <- ic * 10000L + rep
      
      true_params <- generate_true_params(
        p = p, density = cond$density,
        strength = cond$strength, seed = rep_seed
      )
      
      data <- generate_kim_data(
        true_params, n_persons = cond$n_persons,
        n_time = cond$n_time, seed = rep_seed + 5000L
      )
      
      fit <- tryCatch(
        fitKIM(data, vars = true_params$vars,
               id = "id", time = "time",
               verbose = FALSE),
        error = function(e) NULL
      )
      
      if (is.null(fit)) {
        if (verbose) message("    rep ", rep, " 拟合失败，跳过")
        next
      }
      
      metrics <- compute_recovery_metrics(
        true_params$J, fit$J,
        true_params$h, fit$h
      )
      
      metrics$density          <- cond$density
      metrics$n_persons        <- cond$n_persons
      metrics$n_time           <- cond$n_time
      metrics$strength         <- cond$strength
      metrics$rep              <- rep
      # true_rho / est_rho 均使用加权 Jacobian 谱半径 rho_star
      # （正确的 KIM 稳定性指标，已考虑 sigmoid 压缩）
      metrics$true_rho_star    <- true_params$true_rho_star
      metrics$true_rho_J       <- true_params$true_rho_J
      metrics$est_rho_star     <- fit$spectral_radius_star
      metrics$est_rho_J        <- fit$spectral_radius_J
      # 向后兼容：保留 true_rho / est_rho 字段，均指 rho_star
      metrics$true_rho         <- true_params$true_rho_star
      metrics$est_rho          <- fit$spectral_radius_star
      metrics$hub_boost_actual <- true_params$hub_boost_actual
      
      all_results[[idx]] <- metrics
      
      if (verbose && rep %% 10 == 0)
        message("    rep ", rep, "/", n_reps, " done")
    }
  }
  
  results <- do.call(rbind, Filter(Negate(is.null), all_results))
  
  agg_vars <- c("J_RMSE", "J_cor", "sensitivity", "specificity",
                "F1", "sign_accuracy", "h_RMSE", "h_cor")
  
  summary_df <- aggregate(
    results[, agg_vars],
    by = list(
      density   = results$density,
      n_persons = results$n_persons,
      n_time    = results$n_time,
      strength  = results$strength
    ),
    FUN = function(x) round(mean(x, na.rm = TRUE), 3)
  )
  
  summary_sd <- aggregate(
    results[, agg_vars],
    by = list(
      density   = results$density,
      n_persons = results$n_persons,
      n_time    = results$n_time,
      strength  = results$strength
    ),
    FUN = function(x) round(sd(x, na.rm = TRUE), 3)
  )
  
  if (verbose) {
    message("\n═══ Study 1 完成 ═══")
    message("  成功运行数：", nrow(results))
    message("\n汇总（均值）：")
    print(summary_df)
  }
  
  list(
    results    = results,
    summary    = summary_df,
    summary_sd = summary_sd,
    conditions = conditions,
    params     = list(p = p, n_reps = n_reps)
  )
}


#' 绘制 Study 1 参数恢复结果热力图
#' @export
plot_study1 <- function(study1_result, metric = "J_cor") {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("需要 ggplot2")

  df            <- study1_result$summary
  df$condition  <- paste0(df$density, " / ", df$strength)
  df$sample     <- paste0("n=", df$n_persons, ", T=", df$n_time)

  ggplot2::ggplot(df, ggplot2::aes(
    x = .data$sample, y = .data$condition,
    fill = .data[[metric]]
  )) +
    ggplot2::geom_tile(color = "white", linewidth = 0.5) +
    ggplot2::geom_text(ggplot2::aes(label = round(.data[[metric]], 2)),
                       size = 3.5, color = "black") +
    ggplot2::scale_fill_gradient2(
      low = "#D73027", mid = "#FFFFBF", high = "#1A9850",
      midpoint = median(df[[metric]], na.rm = TRUE)
    ) +
    ggplot2::labs(
      title = '',
      x = "Sample Size", y = "Network Condition", fill = metric
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      plot.title  = ggplot2::element_text(face = "bold")
    )
}

