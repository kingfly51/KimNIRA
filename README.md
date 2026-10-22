# KimNIRA 0.4.0

**Kinetic Ising Model + NodeIdentifyR Algorithm**  
用于**二元密集纵向数据**（如经验取样法 EMA / ESM）的动力学网络建模与稳态干预分析。

KimNIRA 把"网络中心性"替换为"**干预效应**"：先拟合 kinetic Ising 模型（一阶、同步、带随机效应的滞后 logistic 转移模型），再用**精确稳态分布**或**配对蒙特卡洛**模拟对每个节点施加持续干预，按干预后系统平均激活率的下降幅度对靶点排序。

---

## 目录

- [安装](#安装)
- [数据要求](#数据要求)
- [快速上手](#快速上手)
- [拟合结果里最该先看的字段](#拟合结果里最该先看的字段)
- [完整工作流](#完整工作流)
- [导出函数速查表](#导出函数速查表)
- [重要说明与常见坑](#重要说明与常见坑)
- [参考教程](#参考教程)
- [版本说明](#版本说明)

---

## 安装

### 1. 安装依赖

```r
install.packages(c("lme4", "ggplot2", "tidyr", "graphics", "grDevices", "stats", "utils"))

# 可选：画网络图、并行 bootstrap、测试
install.packages(c("qgraph", "Matrix", "dplyr", "parallel", "testthat"))
```

### 2. 安装 KimNIRA

**从源码 tarball 安装：**

```r
install.packages("KimNIRA_0.4.0.tar.gz", repos = NULL, type = "source")
```

**或在命令行：**

```bash
R CMD INSTALL KimNIRA_0.4.0.tar.gz
```

**Windows 二进制（.zip）：**

```r
install.packages("KimNIRA_0.4.0.zip", repos = NULL, type = "win.binary")
```

### 3. 开发模式（不安装，直接加载源码）

```r
pkgload::load_all("path/to/KimNIRA")
```

> 注意：开发模式下 `kimStabilityBootstrap(parallel = TRUE)` 需要 worker 能
> `library(KimNIRA)`，即必须安装同版本；否则请用 `parallel = FALSE`。

### 4. 验证安装

```r
library(KimNIRA)
packageVersion("KimNIRA")   # 0.4.0
```

---

## 数据要求

输入是**长格式**数据框，每行一次测量：

| 要求 | 说明 |
|---|---|
| 节点变量 | 多列 **0/1 二值**数值变量（缺失用 `NA`） |
| 个体 ID | 一列被试标识（`id`） |
| 时间索引 | `dayvar` + `beepvar`（推荐），或全局 `time` 列 |
| 观测密度 | 每人至少 `min_obs`（默认 10）个有效滞后对 |
| 转移类型 | 一步 = 相邻两次观测；`allow_overnight = FALSE` 只保留日内转移 |

`fitKIM()` 会报告 `transition_types`，明确"一步"到底是什么间隔。**如果混有多种转移，多步轨迹的解释要谨慎。**

---

## 快速上手

```r
library(KimNIRA)
library(readxl)

ema_long <- read_excel("ema_long.xlsx")
symptoms <- c("Dep", "Loneliness", "RIS", "Stress", "SA", "Anx")

fit <- fitKIM(
  data            = ema_long,
  vars            = symptoms,
  id              = "BH",
  dayvar          = "dayvar",
  beepvar         = "beepvar",
  allow_overnight = FALSE,      # 只保留日内转移
  between_person  = "re",       # 随机截距（默认）
  screen          = "fdr",      # BH-FDR 边筛选（默认）
  fdr_q           = 0.10,
  estimate_between = TRUE,      # 顺便估计 between 网络
  verbose         = TRUE
)

# 对每个节点做精确稳态干预，按缓解效应排序
res <- runKIMNIRA(
  fit,
  method             = "exact",
  intervention_type  = "alleviating",
  intervention_sd    = 2,       # δ = 2 SD(h)，Lunansky et al. (2022) 默认
  verbose            = FALSE
)

res$ranked_targets
res$summary[, c("node", "benefit", "direct_benefit", "spillover_benefit", "rank")]
```

---

## 拟合结果里最该先看的字段

| 字段 | 含义 |
|---|---|
| `fit$temporal_network` | 时间转移网络 **J**；`J[to, from]` = from 在 t 对 to 在 t+1 的影响 |
| `fit$between_network` | 个体均值层面的 between 网络（`qgraph::EBICglasso`） |
| `fit$h` | 节点截距（总体平均 logit 水平） |
| `fit$support` | p×p 逻辑矩阵：哪些边通过筛选 |
| `fit$stage1_beta` / `stage1_se` | 筛选阶段（未惩罚）的系数与标准误 |
| `fit$screen` / `screen_threshold` | 实际使用的筛选方式与阈值 |
| `fit$n_lag_pairs` | 实际用于估计的滞后对数量 |
| `fit$transition_types` | 保留下来的转移类型；确认"一步"是什么间隔 |
| `fit$person_sd` | 每个节点个体截距的 SD |
| `fit$person_constant` | 哪些人-节点结果恒定（截距主要靠收缩估计） |
| `fit$node_fit_failed` | 哪些节点拟合失败 |
| `fit$spectral_radius_star` | 线性化谱半径；≥ 1 提示可能接近临界状态 |

---

## 完整工作流

### 1. 拟合模型

```r
fit <- fitKIM(
  data = ema_long, vars = symptoms, id = "BH",
  dayvar = "dayvar", beepvar = "beepvar",
  allow_overnight = FALSE,
  between_person = "re",        # 或 "res"（随机截距 + 随机斜率）
  screen = "fdr", fdr_q = 0.10,
  estimate_between = TRUE, verbose = TRUE
)
```

**`between_person`：**

- `"re"`：每个节点一个随机截距（默认，部分 pooling）。
- `"res"`：在随机截距基础上加随机斜率（每人有自己的转移矩阵）。6 个节点以上默认 `slope_structure = "diag"`。

**`screen`：**

- `"fdr"`（默认）：BH-FDR 控制，阈值 `fdr_q`；
- `"significance"`：Bonferroni 校正，阈值 `alpha`；
- `"none"`：不筛选，J 为稠密矩阵（兼容旧版行为）。

### 2. 网络图与中心性

```r
p_temporal <- kimPlotNetwork(fit, type = "temporal")
p_between  <- kimPlotNetwork(fit, type = "between")
kimPlotAllNetworks(fit)                 # 两张图并排

cent <- kimCentrality(fit, use_ggplot = TRUE)
cent$temporal   # InStrength / OutStrength / InEI / OutEI
cent$between    # Strength / Closeness / Betweenness / ExpectedInfluence
```

> `kimPlotNetwork()` 依赖 `qgraph`；未安装时会提示。

### 3. 干预排序（NIRA）

**固定 logit 剂量（默认）：**

```r
res    <- runKIMNIRA(fit, method = "exact", intervention_type = "alleviating",
                     intervention_sd = 2, verbose = FALSE)
res_mc <- runKIMNIRA(fit, method = "mc",    intervention_type = "alleviating",
                     intervention_sd = 2, verbose = FALSE)

res$ranked_targets
res$summary[, c("node", "benefit", "direct_benefit", "spillover_benefit", "rank")]
```

**等概率剂量（每个靶点自身稳态激活概率下降 `target_dp`）：**

```r
res_prob <- runKIMNIRA(
  fit,
  method            = "exact",
  dose_scale        = "probability",
  target_dp         = 0.10,
  calibrate_on      = "stationary",   # 或 "intercept"
  intervention_type = "alleviating",
  verbose           = FALSE
)
```

`exact` 精确求解有限状态链的稳态分布；`mc` 用配对蒙特卡洛近似。节点数较多（默认 `p > 12`）时精确引擎会拒绝运行，可改用 `method = "mc"` 或提高 `options(KimNIRA.max_exact_p = ...)`。

### 4. 路径分解：为什么这个节点排第一

```r
dc <- kimPathwayDecomposition(res)
print(dc[, c("node", "effect_full", "direct", "spillover", "propagation",
             "feedback", "self_persistence", "out_strength", "out_ei")],
      digits = 3, row.names = FALSE)
```

| 列 | 含义 |
|---|---|
| `effect_full` | 完整网络的效应（负号 = 缓解） |
| `direct` | 目标节点自身变化 / p |
| `spillover` | 其余节点变化 / p |
| `propagation` | 删除目标出边后效应减少的量 |
| `feedback` | 删除指向目标的入边后效应减少的量 |
| `self_persistence` | 删除目标自回归边后效应减少的量 |
| `rank_nira` / `rank_out_strength` / `half_time` | 排序对照与半效应时间 |

### 5. NIRA 图与轨迹图

```r
p_nira <- kimPlotNIRA(res, perturbation_type = "alleviating")
p_nira_spill <- kimPlotNIRA(res, perturbation_type = "alleviating",
                            effect = "spillover")

# MC 结果可用 p_adjust 控制多重比较（默认 "BH"）
p_nira_mc <- kimPlotNIRA(res_mc, perturbation_type = "alleviating",
                         p_adjust = "bonferroni")

p_summary <- kimPlotTrajectory(res, type = "summary")
p_node    <- kimPlotTrajectory(res, type = "node", target = res$ranked_targets[1])
p_diff    <- kimPlotTrajectory(res, type = "diff")
p_panel   <- kimPlotTrajectory(res, type = "panel")
p_spag    <- kimPlotTrajectory(res_mc, type = "spaghetti", n_traces = 50)
```

> `type = "spaghetti"` 只有 MC 结果支持；exact 结果没有随机轨迹。
> `type = "node"` 图中**蓝色实线 = 干预靶点**（图例显示 `Baseline / Intervention / Target`）。

### 6. 稳定性（被试重抽样 bootstrap）

```r
boot <- kimStabilityBootstrap(
  obj      = res,
  data     = ema_long,
  n_boot   = 100,          # 正式分析建议 ≥ 1000
  seed     = 2025,
  parallel = TRUE,
  ncores   = 6,
  progress = TRUE
)

kimPlotStability(boot, high = length(symptoms))
kimDeltaPlot(boot)         # 效应柱状图 + 95% 区间 + 星号
```

这是**抽样不确定性**的正确来源（对被试重抽样并重新拟合）。`parallel = TRUE` 要求 worker 已安装同版本 KimNIRA。

### 7. 蒙特卡洛误差（数值精度）

```r
mc_err <- kimMCError(res_mc, n_runs = 20)
mc_err$stat    # effect / mc_sd / mc_se_mean / ±2 SD 区间
mc_err$plot    # 柱状图 + 误差棒
```

> `kimMCError` 衡量的是**一次 MC 运行的随机波动**，不是抽样不确定性。

### 8. 模型诊断

```r
# 混合 / 弛豫时间
mix <- kimMixing(fit)
print(mix)                 # 关注 relaxation；决定多步轨迹需要多少步

# 稳态景观 / 双稳性
land <- kimBistability(fit)
print(land)

kimPlotLandscape(land, per_state = TRUE, type = "probability")

# 同期依赖（条件独立假设）
cont <- kimContemporaneousCheck(fit, ema_long)
print(cont)

# 一阶马尔可夫检验（预测性敏感性分析）
mk <- kimMarkovCheck(fit, ema_long)
print(mk)
```

`land$count_distribution`、`land$basins`、`land$spectral`、`land$interventions` 分别给出激活数分布、盆地质量与首达时间、谱证据、各干预对景观的影响。

### 9. 剂量-反应与仿真验证

```r
dose <- kimDoseResponse(fit, target_dp = 0.10, calibrate_on = "stationary")
kimPlotDoseResponse(dose)

# 仿真研究：已知真实转移矩阵时评估边筛选
kimSelectionStats(fit, truth = J_true)

# 随机排序下各恢复指标的期望（机会水平）
kimChanceLevel(p = length(symptoms), k = 3L)
```

### 10. 底层工具（一般不直接调用）

```r
kim_compute_delta(fit, intervention_sd = 2, intervention_type = "alleviating")
kim_delta_prob_table(fit, delta)
kim_init_states(h, n_sim = 5000)
kim_step(state, J, h)
kim_trajectory(J, h, n_steps = 30, n_sim = 5000)
```

---

## 导出函数速查表

| 函数 | 作用 |
|---|---|
| `fitKIM()` | 拟合 KIM（随机截距/随机斜率 + 边筛选） |
| `kimSelectionStats()` | 对照真实网络评估边筛选（TP/FP/FN/TN、灵敏度、特异度、F1） |
| `runKIMNIRA()` | 精确或 MC 的 NIRA 干预排序 |
| `kim_compute_delta()` | 计算固定 logit 剂量的截距位移 δ |
| `kim_delta_prob_table()` | 各节点激活概率变化表 |
| `kim_init_states()` | 按独立分布抽样初始状态 |
| `kim_step()` | KIM 单步状态更新 |
| `kim_trajectory()` | 模拟完整 KIM 轨迹 |
| `kimStabilityBootstrap()` | 被试重抽样稳定性（含 p 值） |
| `kimMCError()` | MC 运行间误差 |
| `kimMarkovCheck()` | 一阶 vs 高阶滞后的预测性检验 |
| `kimCentrality()` | 时间/个体间网络中心性 |
| `kimPlotNetwork()` | 单张网络图（temporal / between） |
| `kimPlotAllNetworks()` | 两张网络图并排 |
| `kimPlotNIRA()` | NIRA 效应连线图（总效应 / 溢出效应） |
| `kimPlotStability()` | 排名稳定性堆叠柱状图 |
| `kimPlotTrajectory()` | 轨迹图（summary / node / diff / panel / spaghetti） |
| `kimDeltaPlot()` | 效应柱状图 + CI + 显著性星号 |
| `kimBistability()` | 稳态景观、盆地、首达时间、谱证据 |
| `kimPlotLandscape()` | 景观图（势能 / 概率） |
| `kimMixing()` | 第二特征值模与弛豫时间 |
| `kimDoseResponse()` | 剂量-反应与等概率剂量排序 |
| `kimPlotDoseResponse()` | 剂量-反应图 |
| `kimPathwayDecomposition()` | 敲除分解（direct / spillover / propagation / feedback / self-persistence） |
| `kimChanceLevel()` | 仿真恢复指标的机会水平 |
| `kimContemporaneousCheck()` | 同期依赖检验 |

另有 5 个 S3 `print` 方法：`kim_nira_result`、`kim_stability_bootstrap`、`kim_markov_check`、`kim_bistability`、`kim_dose_response`。

---

## 重要说明与常见坑

1. **MC 的 p 值 ≠ 抽样推断。** `kimPlotNIRA()` 对 MC 结果的星号来自配对 z 检验，反映的是**蒙特卡洛精度**（轨迹数越多，p 越小）。抽样不确定性请用 `kimStabilityBootstrap()`。
2. **exact 结果没有 p 值、不画星号。** 这是确定性计算，`p_adjust` 会被忽略。
3. **"一步"是观测间隔，不是日历天。** 若排除跨夜转移，模型时间不能直接换算成天数。
4. **景观单峰 ≠ 才能做稳态干预。** 稳态干预对双稳系统同样可算；单峰真正支持的是"稳态均值有代表性、有限时间≈稳态、无阈值响应"。
5. **`type = "node"` 的蓝色线是靶点**；不再用 `[target]` 文字后缀。
6. **bootstrap 会重拟合模型**，耗时较长；先用小 `n_boot` 试跑。
7. **`options(KimNIRA.max_exact_p = 12)`** 是精确引擎的软上限（内存 8·4^p 字节）；超过请用 MC。

---

## 参考教程

完整可运行脚本随包安装，位于：

```r
system.file("tutorials/tutorial_kimNIRA.R", package = "KimNIRA")
```

脚本以 `CC/ema_long.xlsx`（117 名被试、6 个症状、日内 beep 数据）为例，覆盖从拟合、NIRA 排序、路径分解、轨迹图、bootstrap、诊断到拼图的全部流程。

---

## 版本说明

**0.4.0 的主要变化：**

- 估计结构改为 **两阶段**：先做边筛选（BH-FDR / Bonferroni / 不筛选），再对保留的边**无惩罚重估**；
- 移除 L1 罚路径、`glmnet` 依赖、`between_person = "pooled"`；
- 新增 `screen`、`stage1_beta`、`stage1_se`、`support`、`screen_threshold`；
- `kimMarkovCheck()` 不再有惩罚分支；
- `kimPlotNIRA()` 新增 `p_adjust`；
- `kimPlotTrajectory(type = "node")` 用蓝色线标识靶点并加图例；
- 并行 bootstrap 直接 `library(KimNIRA)`，不再依赖 `pkgload` 开发态加载。

详见 `NEWS.md`。

---

## License

GPL-3
