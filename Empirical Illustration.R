devtools::install_github("kingfly51/KimNIRA")

library(KimNIRA)
library(readxl)

help(package="KimNIRA")

#载入数据
ema_long <- read_excel("ema_long.xlsx")

#估计网络
fit_net1<- fitKIM(data = ema_long,
                    vars = c("Dep","Loneliness","RIS","Stress","SA","Anx"),
                    id            = "BH",
                    time          = "day_num",
                    lags             = 1L,
                    lambda_rule      = "lambda.min",
                    n_folds          = 5L,
                    penalize_ar      = TRUE,
                    use_offset       = TRUE,
                    min_obs          = 18L,
                    max_time_gap     = 1,
                    gamma            = 0.5,
                    verbose          = TRUE)

kimPlotNetwork(fit_net1)
kimPlotAllNetworks(fit_net1)

fit_net1$temporal_network
fit_net1$between_network

#中心性指数
cent <- kimCentrality(fit_net1)
cent$temporal   # InStrength, OutStrength, InEI, OutEI
cent$between    # Strength, Closeness, Betweenness, EI

cent <- kimCentrality(fit_net1, use_ggplot = TRUE)
cent$plots$temporal


# 调节效应检验

kmt <- kimModerationTest(data = ema_long,
                  kim_fit = fit_net1,
                  n_perm = 2000
                  )
print(kmt)
kmt$significant_moderators


#NIRA
nira <- runKIMNIRA(fit_net1,
                   "alleviating")
nira

p1 <- kimPlotTrajectory(nira,type = "summary")
print(p1)

p2 <- kimPlotTrajectory(nira,type = "node")
print(p2)

p2_1 <- kimPlotTrajectory(nira,type = "diff")
print(p2_1)

p2_2 <- kimPlotTrajectory(nira,type = "panel")
print(p2_2)

p2_3 <- kimPlotTrajectory(nira,type = "spaghetti")
print(p2_3)


# 4. Significance testing
perm <- kimPermutationNIRAtest(nira)
p3 <- kimPlotNIRA(perm, "alleviating")


# 5. Stability assessment
stab <- kimStabilityNIRAtest(fit_net1, n_runs = 100)
p4 <- kimPlotStability(stab,high =6)
p4


library(patchwork)
library(ggplot2)
combined_plot <- (p1 + p2) / (p3 + p4) +
  plot_annotation(
    tag_levels = 'A',                    # 添加子图标签 (A, B, C, D)
    title = "",
    theme = theme(plot.title = element_text(hjust = 0.9, face = "bold"))
  ) &
  theme_minimal(base_size = 10)          # 统一主题样式（可选）

# 显示组合图
print(combined_plot)
