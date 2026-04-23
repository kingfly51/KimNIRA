source("D:/Rdaima/dynaNIRA/dytemp_nira_1.0/KimNIRA_package/KimNIRA/R/fitKIM.R")
source("D:/Rdaima/dynaNIRA/dytemp_nira_1.0/KimNIRA_package/KimNIRA/R/runKIMNIRA.R")
source("D:/Rdaima/dynaNIRA/dytemp_nira_1.0/study1/simulationStudy_fixed.R")


s1_4 <- run_study1_parameter_recovery(p = 4L, n_reps = 20L, verbose = TRUE)
s1_4
plot_study1(s1_4, metric = "J_cor")
saveRDS(s1_4, file = "D:/Rdaima/dynaNIRA/dytemp_nira_1.0/study1/s1_4_parameter_recovery_p4.rds")


s1_10 <- run_study1_parameter_recovery(p = 10L, n_reps = 20L, verbose = TRUE)
plot_study1(s1_10, metric = "J_cor")
s1_10
saveRDS(s1_10, file = "D:/Rdaima/dynaNIRA/dytemp_nira_1.0/study1/s1_4_parameter_recovery_p10.rds")


s1_20 <- run_study1_parameter_recovery(p = 20L, n_reps = 20L, verbose = TRUE)
plot_study1(s1_20, metric = "J_cor")
s1_20
saveRDS(s1_20, file = "s1_4_parameter_recovery_p20.rds")



plot_study1_combined <- function(study1_results_list, metric = "J_cor") {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("需要 ggplot2")
  
  # 合并多个结果，添加 p 标识
  df_list <- list()
  for (p_val in names(study1_results_list)) {
    df <- study1_results_list[[p_val]]$summary
    df$p_vars <- paste0("p = ", p_val)
    df_list[[p_val]] <- df
  }
  df <- do.call(rbind, df_list)
  
  # 转换为有序因子，按数字顺序排列
  p_order <- as.character(sort(as.numeric(names(study1_results_list))))
  df$p_vars <- factor(df$p_vars, levels = paste0("p = ", p_order))
  
  df$condition <- paste0(df$density, " / ", df$strength)
  df$sample    <- paste0("n=", df$n_persons, ", T=", df$n_time)
  
  ggplot2::ggplot(df, ggplot2::aes(
    x = .data$sample, y = .data$condition,
    fill = .data[[metric]]
  )) +
    ggplot2::geom_tile(color = "white", linewidth = 0.5) +
    ggplot2::geom_text(ggplot2::aes(label = round(.data[[metric]], 2)),
                       size = 3, color = "black") +
    ggplot2::scale_fill_gradient2(
      low = "#D73027", mid = "#FFFFBF", high = "#1A9850",
      midpoint = median(df[[metric]], na.rm = TRUE)
    ) +
    ggplot2::facet_wrap(~ p_vars, ncol = 3) +
    ggplot2::labs(
      title = 'Parameter Recovery Heatmap',
      x = "Sample Size (n persons, T time points)", 
      y = "Network Condition", 
      fill = metric
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      plot.title  = ggplot2::element_text(face = "bold", hjust = 0.5),
      strip.background = ggplot2::element_rect(fill = "lightgray", color = NA),
      strip.text = ggplot2::element_text(face = "bold")
    )
}

# 使用时注意列表名称是字符格式
results_combined <- list("4" = s1_4, "10" = s1_10, "20" = s1_20)
plot_study1_combined(results_combined, metric = "J_cor")