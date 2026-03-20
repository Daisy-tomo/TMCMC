# TMCMC 航空发动机热力循环参数贝叶斯反演

## 程序说明

文件 `TMCMC_engine_inversion.m` 实现了基于 **TMCMC（Transitional Markov Chain Monte Carlo）** 算法的贝叶斯后验反演，用于从航空发动机整机性能观测量（比推力 $R_{ud}$、比油耗 $C_{ud}$）反推热力循环关键修正参数的后验分布。

---

## 问题描述

| 类型 | 内容 |
|------|------|
| **待反演参数** | $\eta_k$（压气机效率）、$\eta_t$（涡轮气动效率）、$\eta_T$（涡轮膨胀效率）、$\eta_m$（机械效率） |
| **观测量** | 比推力 $R_{ud}\ [\text{N·s/kg}]$、比油耗 $C_{ud}\ [\text{kg/(N·h)}]$ |
| **先验** | 四参数均匀先验（见下表） |
| **似然** | 独立高斯似然，噪声水平 1% |

先验范围：

| 参数 | 下界 | 上界 | 真值（验证用） |
|------|------|------|----------------|
| $\eta_k$ | 0.840 | 0.860 | 0.860 |
| $\eta_t$ | 0.860 | 0.920 | 0.880 |
| $\eta_T$ | 0.980 | 0.995 | 0.988 |
| $\eta_m$ | 0.970 | 0.990 | 0.980 |

---

## TMCMC 算法原理

TMCMC 由 Ching & Chen (2007) 提出，通过构造一系列中间分布
$$
p_j(\boldsymbol\theta) \propto p(\boldsymbol\theta)\cdot[p(\mathbf{D}|\boldsymbol\theta)]^{\beta_j}, \quad 0=\beta_0<\beta_1<\cdots<\beta_J=1
$$
将先验 $p_0$ 逐步过渡到后验 $p_J$，避免了直接 MCMC 在高维/多峰后验中的混合困难。

### 算法流程

```
Stage 0
  └─ 从均匀先验采样 N 个粒子 {θ_i}，计算对数似然 log L(θ_i)

Stage j = 1, 2, ...  (直到 β_J = 1)
  ├─ [求 β_j]  二分法搜索 β_j，使重要性权重的变异系数（COV）= 目标值（默认 1.0）
  ├─ [加权]    w_i = exp((β_j - β_{j-1}) · log L(θ_i))，归一化得 w̃_i
  ├─ [重采样]  按 w̃ 多项式重采样 N 个粒子（bootstrap）
  ├─ [协方差]  计算加权样本协方差 Σ_j = scale² · Cov({θ_i}, w̃)
  └─ [MCMC]   以 Σ_j 为提议，从每个粒子出发做 N_MCMC 步 MH 采样
               目标分布：π_j(θ) ∝ 均匀先验 × L(θ)^β_j
```

### 关键参数说明

| 参数 | 含义 | 推荐值 |
|------|------|--------|
| `N` | 粒子数 | 1000–3000 |
| `COV_target` | 权重变异系数目标（越小阶段越多） | 0.5–1.5 |
| `N_MCMC` | 每粒子每阶段 MCMC 步数 | 3–10 |
| `scale` | 提议协方差缩放因子（$c_v$） | 0.15–0.30 |

---

## 程序结构

```
TMCMC_engine_inversion.m
│
├─ 主脚本（第 1–120 行）
│   ├─ 参数定义与初始化
│   ├─ 前向模型验证
│   ├─ 生成虚拟观测数据
│   ├─ TMCMC 配置与运行
│   ├─ 后验统计输出
│   └─ 绘图调用
│
└─ 局部函数
    ├─ engine_forward        — 发动机热力循环前向模型（→ R_ud, C_ud）
    ├─ piecewise_kT          — 分段燃气绝热指数
    ├─ piecewise_RT          — 分段燃气气体常数
    ├─ delta_cooling         — 涡轮冷却引气系数
    ├─ generate_virtual_data — 生成含高斯噪声的虚拟观测
    ├─ log_likelihood        — 对数似然（独立高斯）
    ├─ find_next_beta        — 二分法求下一个温度参数 β
    ├─ compute_cov           — 计算重要性权重的变异系数
    ├─ run_tmcmc             — TMCMC 主算法
    └─ plot_results          — 绘制三幅诊断图
```

---

## 输出说明

### 控制台输出

- 前向模型验证结果（$R_{ud}$、$C_{ud}$）
- 虚拟观测值与噪声水平
- 每个 TMCMC 阶段的 $\beta$ 演化与 MCMC 接受率
- 后验统计表：真值 / 后验均值 / MAP / 95% CI
- 后验预测与真值的相对误差

### 图形输出

| 图片文件 | 内容 |
|----------|------|
| `TMCMC_marginal_posterior.png` | 四个参数的边缘后验直方图（含真值线、均值线、MAP线） |
| `TMCMC_joint_posterior.png` | Corner plot：对角=边缘分布，下三角=散点图，上三角=Pearson 相关系数 |
| `TMCMC_diagnostics.png` | β 演化曲线 + 各阶段 MCMC 接受率 |

---

## 运行方法

1. 将 `TMCMC_engine_inversion.m` 复制到 MATLAB 工作目录
2. 在 MATLAB 命令窗口运行：
   ```matlab
   run('TMCMC_engine_inversion.m')
   ```
   或直接打开文件按 **F5** 运行

3. 典型运行时间（N=2000，COV_target=1.0，N_MCMC=3）：**约 30–90 秒**（普通笔记本）

> **提示**：如需更快运行，可将 `N` 降至 `1000`，或将 `N_MCMC` 降至 `2`。
> 如需更精确结果，可将 `N` 提高至 `3000–5000`，`N_MCMC` 提高至 `5`。

---

## 与 `san.m`（标准 MH-MCMC）的对比

| 对比项 | `san.m`（Metropolis-Hastings） | `TMCMC_engine_inversion.m`（TMCMC） |
|--------|-------------------------------|--------------------------------------|
| 算法核心 | 单链自适应 MH | 粒子群过渡采样 |
| 烧入处理 | 需手动设置 burn-in | 无需 burn-in（过渡过程自然收敛） |
| 后验代表性 | 单链，依赖混合 | N 个粒子并行，覆盖更全面 |
| 多峰鲁棒性 | 较差 | 较好（多粒子探索） |
| 对数证据 | 不直接输出 | 输出 `log_evidence`（可用于模型选择） |
| 联合分布图 | 无 | Corner plot（散点 + 相关系数） |
| 推荐场景 | 快速调试 | 正式后验分析 |

---

## 参考文献

> Ching, J., & Chen, Y. (2007). Transitional Markov Chain Monte Carlo Method for Bayesian Model Updating, Model Class Selection, and Model Averaging. *Journal of Engineering Mechanics*, 133(7), 816–832. https://doi.org/10.1061/(ASCE)0733-9399(2007)133:7(816)
