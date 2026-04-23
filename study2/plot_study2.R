# =============================================================================
# plot_study2.R  ——  Study 2 三指标结果可视化（出版级）
#
# 用法：
#   source("plot_study2.R")
#   fig <- plot_study2_results(s2, s2_10, s2_20)
#   ggsave("study2_results.pdf", fig, width = 16, height = 11, dpi = 300)
# =============================================================================

#' 绘制 Study 2 三指标结果图（p=4/10/20 合并）
#'
#' @param s2    p=4  的 run_study2_target_recovery() 返回值
#' @param s2_10 p=10 的返回值
#' @param s2_20 p=20 的返回值
#' @return ggplot 对象
plot_study2_results <- function(s2, s2_10, s2_20) {

  for (pkg in c("ggplot2","dplyr","tidyr","forcats","scales")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      stop("请先安装：install.packages('", pkg, "')")
  }
  library(ggplot2); library(dplyr); library(tidyr)
  library(forcats); library(scales)

  # ── 1. 整理数据 ──────────────────────────────────────────────────────────────
  prep <- function(s, p_val) {
    s$summary %>%
      mutate(
        p_lab     = factor(paste0("p = ", p_val),
                           levels = c("p = 4","p = 10","p = 20")),
        n_lab     = factor(paste0("n = ", n_persons),
                           levels = c("n = 50","n = 100","n = 200")),
        T_lab     = factor(paste0("T = ", n_time),
                           levels = c("T = 14","T = 30","T = 60")),
        signal    = factor(ifelse(strength == "strong",
                                  "Strong signal","Weak signal"),
                           levels = c("Strong signal","Weak signal")),
        density_f = factor(ifelse(density == "medium",
                                   "Medium density","Sparse density"),
                           levels = c("Medium density","Sparse density"))
      )
  }

  df_all <- bind_rows(prep(s2, 4), prep(s2_10, 10), prep(s2_20, 20))

  # 长格式：三指标
  df_long <- df_all %>%
    pivot_longer(cols      = c(top1_hit, top3_overlap, kendall_tau),
                 names_to  = "metric",
                 values_to = "value") %>%
    mutate(
      metric = factor(metric,
                      levels = c("top1_hit","top3_overlap","kendall_tau"),
                      labels = c("Top-1 Hit Rate",
                                 "Top-3 Overlap",
                                 "Kendall's \u03c4"))
    )

  # ── 2. 配色与主题 ─────────────────────────────────────────────────────────────
  # p 用颜色，density 用线型，T 用 x 轴，n 用子图列
  p_colors <- c("p = 4"  = "#1D6FA4",
                "p = 10" = "#E07B28",
                "p = 20" = "#2CA25F")

  density_lt <- c("Medium density" = "solid",
                  "Sparse density" = "dashed")

  signal_strip_fill <- c("Strong signal" = "#1D4E89",
                          "Weak signal"   = "#8B1A1A")

  # ── 3. 绘图 ─────────────────────────────────────────────────────────────────
  # 分面结构：行 = metric × signal，列 = n
  # x 轴 = T，颜色 = p，线型 = density
  g <- ggplot(df_long,
              aes(x        = T_lab,
                  y        = value,
                  colour   = p_lab,
                  linetype = density_f,
                  group    = interaction(p_lab, density_f))) +

    # ── 参考线 ──
    geom_hline(yintercept = 1/3, colour = "grey80",
               linewidth = 0.4, linetype = "dotted") +
    geom_hline(yintercept = 1,   colour = "grey70",
               linewidth = 0.4, linetype = "dotted") +

    # ── 误差带（±1 SD，需要 sd 数据） ──
    # （此处先略去，若需要可以 join summary_sd）

    # ── 折线 ──
    geom_line(linewidth = 0.9, alpha = 0.9) +

    # ── 点 ──
    geom_point(size = 2.4, alpha = 0.95, stroke = 0.3) +

    # ── 分面：行 = metric + signal，列 = n ──
    facet_grid(
      metric + signal ~ n_lab,
      scales = "free_y",
      switch = "y"
    ) +

    # ── 坐标轴 ──
    scale_y_continuous(
      limits = c(-0.02, 1.05),
      breaks = c(0, 0.25, 0.5, 0.75, 1.0),
      labels = c("0","0.25","0.50","0.75","1.00"),
      expand = c(0, 0)
    ) +
    scale_x_discrete(expand = c(0.15, 0)) +

    # ── 颜色、线型 ──
    scale_colour_manual(
      name   = "Network size (p)",
      values = p_colors,
      guide  = guide_legend(
        title.position = "top",
        override.aes   = list(linewidth = 1.5, size = 3)
      )
    ) +
    scale_linetype_manual(
      name   = "Network density",
      values = density_lt,
      guide  = guide_legend(
        title.position = "top",
        override.aes   = list(linewidth = 1.2, colour = "grey30")
      )
    ) +

    # ── 标签 ──
    labs(
      title   = "Study 2: Intervention Target Recovery",
      subtitle= paste0(
        "",
        ""
      ),
      x       = "Time points per participant (T)",
      y       = NULL,
      caption = paste0(
        "",
        ""
      )
    ) +

    # ── 主题 ──
    theme_minimal(base_size = 11.5) +
    theme(
      # 标题
      plot.title    = element_text(face   = "bold", size  = 15,
                                   colour = "#1A1A2E",
                                   margin = margin(b = 3)),
      plot.subtitle = element_text(size   = 9,  colour = "grey40",
                                   margin = margin(b = 10)),
      plot.caption  = element_text(size   = 8,  colour = "grey55",
                                   hjust  = 0,  margin = margin(t = 8)),

      # 分面标签
      strip.text.y.left = element_text(
        face   = "bold", size = 9, angle = 90,
        colour = "grey20", margin = margin(r = 4)
      ),
      strip.text.x = element_text(
        face   = "bold", size = 10.5,
        colour = "grey15", margin = margin(t = 4, b = 4)
      ),
      strip.background = element_rect(fill = "#F0F4F8", colour = NA),
      strip.placement  = "outside",

      # 面板
      panel.background = element_rect(fill = "white", colour = NA),
      panel.grid.major.y = element_line(colour = "#E8ECF0", linewidth = 0.45),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.spacing.x    = unit(1.0, "lines"),
      panel.spacing.y    = unit(0.6, "lines"),
      panel.border       = element_rect(colour = "#D0D8E4",
                                        fill = NA, linewidth = 0.3),

      # 坐标轴
      axis.text.x  = element_text(size = 9,   colour = "grey30"),
      axis.text.y  = element_text(size = 8.5, colour = "grey30"),
      axis.title.x = element_text(size = 10,  colour = "grey20",
                                  face = "bold", margin = margin(t = 8)),

      # 图例
      legend.position  = "right",
      legend.title     = element_text(face = "bold", size = 9.5),
      legend.text      = element_text(size = 9),
      legend.key.width = unit(2, "lines"),
      legend.spacing.y = unit(0.8, "lines"),
      legend.box       = "vertical",
      legend.margin    = margin(l = 6),

      # 外边距
      plot.margin      = margin(12, 16, 10, 12),
      plot.background  = element_rect(fill = "#FAFBFC", colour = NA)
    )

  g
}


# =============================================================================
# 附加：带误差带的增强版（同时展示 SD）
# =============================================================================

#' 带 ±1SD 误差带的 Study 2 结果图
#'
#' 需要同时传入 summary 和 summary_sd（均来自 run_study2_target_recovery 返回值）
plot_study2_with_sd <- function(s2, s2_10, s2_20) {

  for (pkg in c("ggplot2","dplyr","tidyr","forcats","scales")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      stop("请先安装：install.packages('", pkg, "')")
  }
  library(ggplot2); library(dplyr); library(tidyr)
  library(forcats); library(scales)

  prep_both <- function(s, p_val) {
    grp_vars <- c("density","n_persons","n_time","strength")

    mean_df <- s$summary %>%
      pivot_longer(c(top1_hit, top3_overlap, kendall_tau),
                   names_to = "metric", values_to = "mean_val")
    sd_df <- s$summary_sd %>%
      pivot_longer(c(top1_hit, top3_overlap, kendall_tau),
                   names_to = "metric", values_to = "sd_val")

    left_join(mean_df, sd_df,
              by = c(grp_vars, "metric")) %>%
      mutate(
        p_lab     = factor(paste0("p = ", p_val),
                           levels = c("p = 4","p = 10","p = 20")),
        n_lab     = factor(paste0("n = ", n_persons),
                           levels = c("n = 50","n = 100","n = 200")),
        T_lab     = factor(paste0("T = ", n_time),
                           levels = c("T = 14","T = 30","T = 60")),
        signal    = factor(ifelse(strength == "strong",
                                  "Strong signal","Weak signal"),
                           levels = c("Strong signal","Weak signal")),
        density_f = factor(ifelse(density == "medium",
                                   "Medium density","Sparse density"),
                           levels = c("Medium density","Sparse density")),
        metric    = factor(metric,
                           levels = c("top1_hit","top3_overlap","kendall_tau"),
                           labels = c("Top-1 Hit Rate",
                                      "Top-3 Overlap",
                                      "Kendall's \u03c4")),
        ylo = pmax(mean_val - sd_val, 0),
        yhi = pmin(mean_val + sd_val, 1)
      )
  }

  df_long <- bind_rows(
    prep_both(s2,    4),
    prep_both(s2_10, 10),
    prep_both(s2_20, 20)
  )

  p_colors <- c("p = 4"  = "#1D6FA4",
                "p = 10" = "#E07B28",
                "p = 20" = "#2CA25F")
  density_lt <- c("Medium density" = "solid",
                  "Sparse density" = "dashed")

  ggplot(df_long,
         aes(x        = T_lab,
             y        = mean_val,
             colour   = p_lab,
             fill     = p_lab,
             linetype = density_f,
             group    = interaction(p_lab, density_f))) +

    geom_hline(yintercept = 1/3, colour = "grey80",
               linewidth = 0.4, linetype = "dotted") +
    geom_hline(yintercept = 1,   colour = "grey70",
               linewidth = 0.4, linetype = "dotted") +

    # 误差带
    geom_ribbon(aes(ymin = ylo, ymax = yhi),
                alpha = 0.12, colour = NA) +

    # 折线 + 点
    geom_line(linewidth = 0.85, alpha = 0.9) +
    geom_point(size = 2.3, stroke = 0.25, alpha = 0.95) +

    facet_grid(metric + signal ~ n_lab,
               scales = "free_y", switch = "y") +

    scale_y_continuous(limits = c(-0.02, 1.05),
                       breaks = c(0, 0.25, 0.5, 0.75, 1),
                       labels = c("0","0.25","0.50","0.75","1.00"),
                       expand = c(0, 0)) +
    scale_x_discrete(expand = c(0.15, 0)) +

    scale_colour_manual(name   = "Network size (p)",
                        values = p_colors,
                        guide  = guide_legend(
                          title.position = "top",
                          override.aes   = list(linewidth=1.5, size=3, fill=NA))) +
    scale_fill_manual(values = p_colors, guide = "none") +
    scale_linetype_manual(name   = "Network density",
                          values = density_lt,
                          guide  = guide_legend(
                            title.position = "top",
                            override.aes   = list(linewidth=1.2,
                                                  colour="grey30"))) +

    labs(
      title    = "Study 2: Intervention Target Recovery",
      subtitle = paste0(
        "Mean \u00b1 1 SD across 50 replications per cell  \u2014  ",
        "rows: metric \u00d7 signal strength; columns: sample size (n)"
      ),
      x       = "Time points per participant (T)",
      y       = NULL,
      caption = "Shaded bands show \u00b11 SD. Dotted lines at chance (1/3) and ceiling (1.0)."
    ) +

    theme_minimal(base_size = 11.5) +
    theme(
      plot.title    = element_text(face="bold", size=15, colour="#1A1A2E",
                                   margin=margin(b=3)),
      plot.subtitle = element_text(size=9, colour="grey40",
                                   margin=margin(b=10)),
      plot.caption  = element_text(size=8, colour="grey55",
                                   hjust=0, margin=margin(t=8)),
      strip.text.y.left = element_text(face="bold", size=9, angle=90,
                                        colour="grey20"),
      strip.text.x      = element_text(face="bold", size=10.5,
                                        colour="grey15",
                                        margin=margin(t=4,b=4)),
      strip.background  = element_rect(fill="#F0F4F8", colour=NA),
      strip.placement   = "outside",
      panel.background  = element_rect(fill="white", colour=NA),
      panel.grid.major.y= element_line(colour="#E8ECF0", linewidth=0.45),
      panel.grid.major.x= element_blank(),
      panel.grid.minor  = element_blank(),
      panel.spacing.x   = unit(1.0,"lines"),
      panel.spacing.y   = unit(0.6,"lines"),
      panel.border      = element_rect(colour="#D0D8E4", fill=NA,
                                        linewidth=0.3),
      axis.text.x  = element_text(size=9,   colour="grey30"),
      axis.text.y  = element_text(size=8.5, colour="grey30"),
      axis.title.x = element_text(size=10,  colour="grey20", face="bold",
                                   margin=margin(t=8)),
      legend.position  = "right",
      legend.title     = element_text(face="bold", size=9.5),
      legend.text      = element_text(size=9),
      legend.key.width = unit(2,"lines"),
      legend.spacing.y = unit(0.8,"lines"),
      legend.box       = "vertical",
      legend.margin    = margin(l=6),
      plot.margin      = margin(12,16,10,12),
      plot.background  = element_rect(fill="#FAFBFC", colour=NA)
    )
}

