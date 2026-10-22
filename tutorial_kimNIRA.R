# =============================================================================
# KimNIRA 使用教程 —— 以 ema_long.xlsx 为例

pkgload::load_all("D:/Rdaima/dynaNIRA/multiple behabiour research/KimNIRA1.2.8/KimNIRA_0.3.1/KimNIRA")



library(readxl)  
ema_long <- read_excel("D:/Rdaima/dynaNIRA/multiple behabiour research/KimNIRA1.2.8/CC/ema_long.xlsx")
#length(unique(ema_long$BH))

symptoms <- c("Dep", "Loneliness", "RIS", "Stress", "SA", "Anx")

##拟合 KIM

fit <- fitKIM(dat              = ema_long,
              vars             = symptoms,
              id               = "BH",
              dayvar           = "dayvar",
              beepvar          = "beepvar",
              allow_overnight  = FALSE,      # 只保留日内转移
              between_person   = "re",
              screen           = "fdr",
              fdr_q            = 0.10,
              estimate_between = TRUE,      
              verbose          = TRUE)

## 拟合结果里最需要先看的字段
fit$temporal_network     # 时间转移网络 J；J[to, from] = from 对 to 的影响
fit$between_network      # 个体均值层面的 between 网络（qgraph::EBICglasso）
fit$h                    # 节点截距（总体平均 logit 水平）
fit$support              # p x p 逻辑矩阵：哪些边通过筛选
fit$n_lag_pairs          # 实际用于估计的滞后对数量
fit$transition_types     # 保留下来的转移类型；确认“一步”是什么间隔
fit$person_sd            # 每个节点个体截距的 SD
fit$person_constant      # 哪些人-节点结果恒定（截距主要靠收缩估计）
fit$node_fit_failed      # 哪些节点拟合失败
fit$spectral_radius_star # 线性化谱半径；>= 1 时提示可能接近临界状态

# 重要：如果 transition_types 显示多种转移（例如 beep 1->2 与 2->1 混合），
# 说明“一步”不是单一间隔；此时多步轨迹的解释要谨慎。



## 网络图与中心性
#单张网络图：type = "temporal" 或 "between"
p_temporal <- kimPlotNetwork(fit, type = "temporal")
p_between  <- kimPlotNetwork(fit, type = "between")

# 同时画时间网络与 between 网络
kimPlotAllNetworks(fit)

# 中心性指标
#   temporal 默认包含 InStrength、OutStrength、InEI、OutEI
#   between  默认包含 Strength、Closeness、Betweenness、ExpectedInfluence
cent <- kimCentrality(fit, use_ggplot = TRUE)

cent$temporal    # 时间网络中心性数据框
cent$between     # between 网络中心性数据框
cent$plots$temporal
cent$plots$between


##精确 NIRA：固定 logit 位移
#基于 2 个标准差（Lunansky et al., 2022 的默认规则），用默认 intervention_sd = 2
res <- runKIMNIRA(fit, method = "exact",
                  intervention_type = "alleviating",
                  intervention_sd = 2, verbose = FALSE)

cat("排名第一的靶点：", res$ranked_targets[1], "\n\n")
print(res$summary[, c("node", "benefit", "direct_benefit", "spillover_benefit", "rank")],
      digits = 4, row.names = FALSE)


##蒙罗卡罗模拟NIRA：固定logit位移
res_mc <- runKIMNIRA(fit, method = "mc",
                  intervention_type = "alleviating",
                  intervention_sd = 2, verbose = FALSE)

cat("排名第一的靶点：", res_mc$ranked_targets[1], "\n\n")
print(res_mc$summary[, c("node", "benefit", "direct_benefit", "spillover_benefit", "rank")],
      digits = 4, row.names = FALSE)



##精确 NIRA：等概率干预
res1 <- runKIMNIRA(fit, method = "exact",
                  dose_scale = "probability",
                  target_dp = 0.10,
                  calibrate_on = "stationary",   # 或 "intercept"
                  intervention_type = "alleviating",
                  verbose = FALSE)

cat("排名第一的靶点：", res1$ranked_targets[1], "\n\n")
print(res1$summary[, c("node", "benefit", "direct_benefit", "spillover_benefit", "rank")],
      digits = 4, row.names = FALSE)




### 路径分解：为什么这个节点排第一
# kimPathwayDecomposition() 给出：
#   effect_full     ：完整网络的效应（负号表示缓解）
#   direct          ：目标节点自身变化 / p
#   spillover       ：其余节点变化 / p
#   propagation     ：删除目标出边后效应减少的量
#   feedback        ：删除指向目标的入边后效应减少的量
#   self_persistence：删除目标自回归边后效应减少的量
#   out_strength / out_ei / autoregressive / baseline_probability
#   rank_nira / rank_out_strength / half_time

dc <- kimPathwayDecomposition(res)
print(dc[, c("node", "effect_full", "direct", "spillover", "propagation",
             "feedback", "self_persistence", "out_strength", "out_ei")],
      digits = 3, row.names = FALSE)


##NIRA 图
p_nira <- kimPlotNIRA(
  res,
  perturbation_type = "alleviating"
)
print(p_nira)


p_n <- kimPlotNIRA(
  res,
  perturbation_type = "alleviating",
  effect =  "spillover"
)
print(p_n)

p_nira_mc <- kimPlotNIRA(
  res_mc,
  perturbation_type = "alleviating"
)
print(p_nira_mc)


##轨迹图

p_traj_summary <- kimPlotTrajectory(res, type = "summary")
p_traj_node    <- kimPlotTrajectory(res, type = "node", target = res$ranked_targets[1])
p_traj_diff    <- kimPlotTrajectory(res, type = "diff")
p_traj_panel   <- kimPlotTrajectory(res, type = "panel")
p_spaghetti <- kimPlotTrajectory(res_mc, type = "spaghetti", n_traces = 50)#exact 结果没有随机轨迹，不支持 spaghetti；MC 结果可以画单条轨迹的“毛线图”


print(p_traj_summary)
print(p_traj_node)
print(p_traj_diff)
print(p_traj_panel)
print(p_spaghetti)


#稳定性检验

boot <- kimStabilityBootstrap(
  obj       = res,
  data      = ema_long,
  n_boot    = 100,          # 先小规模试跑；正式分析改大
  seed      = 2025,
  parallel  = TRUE,        # Windows 下可设 TRUE 并给 ncores
  ncores = 6,
  progress  = TRUE
)

p_stab <- kimPlotStability(boot, high = length(symptoms))
print(p_stab)


## Δ 效应柱状图：bootstrap 的效应与 CI --------------------------------
# kimDeltaPlot 用 kimStabilityBootstrap 的结果画效应柱状图 + 95% 区间 + 星号。
# 与 kimPlotNIRA 不同：这是每个节点效应的柱状图（非“基线→干预”的连线图）。
p_delta <- kimDeltaPlot(boot)   # 默认 stars="auto" -> vs_zero，只标 vs 0 的星号
print(p_delta)



## 蒙特卡洛误差：单次 MC 结果的数值精度
# kimMCError 重复跑 n_runs 次 MC，衡量“一次 runKIMNIRA(method='mc')”本身的
# 随机波动（run-to-run SD）。注意：这是数值精度，不是抽样不确定性；
# 抽样不确定性用 kimStabilityBootstrap 才有意义。
mc_err <- kimMCError(res_mc, n_runs = 20)
print(mc_err$stat)      # effect、mc_sd、mc_se_mean、±2 SD 区间
print(mc_err$plot)      # 柱状图 + 误差棒



##模型诊断
# 混合/弛豫：精确链的第二大特征值模与弛豫时间
mix <- kimMixing(fit)
print(mix)   # 关注 relaxation；它决定多步轨迹/稳态需要多少步



# 双稳/景观：如果关心是否存在多个稳定盆地
land <- kimBistability(fit)
print(land)

p_land2 <- kimPlotLandscape(
  land,
  per_state = TRUE,          # 除以 C(p, K)，去掉组合数导致的中间态虚高
  type      = "probability"  # 直接画 P(K)，而不是势能 U(K)，potential
)
print(p_land2)

# 景观字段说明：
land$count_distribution   # 各条件(基线/干预)下按激活数 K 的稳态概率
land$basins               # 低/高激活盆地及其质量、平均首达时间
land$spectral             # 特征值/谱隙（双稳性的谱证据）
land$interventions        # 各干预目标对景观的影响
#想看干预如何改变景观,用下面的代码
#land2 <- kimBistability(fit, delta = -1, targets = symptoms)
#land2$interventions    # 每个目标干预后的景观变化

# 同期依赖：控制滞后后，其他节点的当期状态是否仍有额外解释力
# 显著不代表模型错误，但说明“仅滞后”的转移模型可能遗漏同期结构。
cont <- kimContemporaneousCheck(fit, ema_long)
print(cont)

# 一阶马尔可夫检验：比较一阶与更高阶滞后
#     注意：每天两次、排除跨夜的数据通常凑不出连续三点，
#     函数会明确报告“无法检验”，这不是错误。
mk <- kimMarkovCheck(fit, ema_long)
print(mk)




###拼图
library(patchwork)
library(ggplot2)

combined_plot <- (p_traj_summary | p_traj_node) / 
  (p_nira         | p_stab) +
  plot_annotation(
    tag_levels = 'A', # 自动打上 A, B, C, D 的标签
    tag_prefix = '(', 
    tag_suffix = ')',
    title = "",
    theme = theme(plot.title = element_text(size = 6, face = "bold", hjust = 0.5))
  )

# 查看拼图
print(combined_plot)
