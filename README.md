# PRF v2

细核组上的三个增量指标。它们描述同一次层级加细里、方向尺度变小的那批高斯，不替代 v1 的全场景 patch 统计。

$$
R_{\mathrm{phy}}^{\mathrm{detail}}=\max(R_{\mathrm{obs}}^{\mathrm{detail}},R_{\mathrm{kernel}}^{\mathrm{detail}},R_{\mathrm{spacing}}^{\mathrm{detail}})
$$

P50 与 P90 在软加权汇总之后分开取三项的较大值。三项使用同一批局部和同一套权重。完整计算式见下方。

## 安装

```bash
pip install numpy scipy
pip install -e .
```

## 调用

```python
from prf_detail import (
    covariance_world, directional_scale, directional_tangent,
    detail_kernel_quantiles, detail_spacing_quantiles, neighbor_spacing,
    summarize_detail,
)

tangent = directional_tangent(tangent_1, tangent_2, direction_deg=0.0)
child = directional_scale(covariance_world(rotations, scales), tangent)
parent = directional_scale(parent_covariances, tangent)
kernel = detail_kernel_quantiles(child, parent, expanded, visible_weights)

spacing, accepted = neighbor_spacing(pool_xy, query_xy, own_index, k=6, eps=eps, margin=margin)
local_spacing = detail_spacing_quantiles(spacing, accepted, visible_weights)

summary = summarize_detail(
    obs_p50, obs_p90, kernel_p50, kernel_p90, spacing_p50, spacing_p90,
    delta_e, eligible, lam=0.5, planned=45,
)
print(summary["obs_p50"], summary["kernel_p50"], summary["spacing_p50"], summary["phy_p50"])
```

## 计算公式

方向 $\theta$ 用 ROI 的两条切向合成，不再单位化：

$$
t=\cos\theta\, t_1+\sin\theta\, t_2.
$$

单个高斯的方向尺度是切向标准差乘 3.77。$\Sigma\_i$ 是该高斯的世界坐标协方差。

$$
\sigma_i=3.77\sqrt{t^{\mathsf T}\Sigma_i t},\qquad
\rho_i=\frac{\sigma_i^{\mathrm{child}}}{\sigma_i^{\mathrm{parent}}}.
$$

展开且 $\rho\_i<0.99$ 的子节点进入细核组。$\rho\_i\in[0.99,1.01]$ 为相近，$\rho\_i>1.01$ 为变粗，未展开为保留。后三类不进入下面的分位数。

组内分位数是可见权重的离散逆 CDF。把细核节点按数值从小到大排列，权重为 $w\_i$，

$$
F(k)=\frac{\sum_{i=1}^{k} w_{(i)}}{\sum_i w_i}.
$$

$Q\_q$ 取第一个满足 $F(k)\ge q$ 的数值。局部组合 $j$ 的核指标是

$$
R_{\mathrm{kernel},j,50}^{\mathrm{detail}}=Q_{0.5}^{\mathrm{visible}}(\sigma_i),\qquad
R_{\mathrm{kernel},j,90}^{\mathrm{detail}}=Q_{0.9}^{\mathrm{visible}}(\sigma_i).
$$

间距在同一张切平面上计算。对细核中心 $u\_i$，在当前层全部同表面中心中取 6 个其他中心；与查询点或已选邻居距离不超过 $\varepsilon$ 的中心不占名额。$\varepsilon$ 默认是 ROI 短边的 $10^{-4}$。

$$
s_i=\frac{1}{6}\sum_{r=1}^{6}\lVert u_i-u_{n_r}\rVert.
$$

不足 6 个时 $s\_i$ 留空，不减小 $k$。第六近邻距离达到切向缓冲边界时，该中心不进入分位数。局部间距是

$$
R_{\mathrm{spacing},j,50}^{\mathrm{detail}}=Q_{0.5}^{\mathrm{visible}}(s_i),\qquad
R_{\mathrm{spacing},j,90}^{\mathrm{detail}}=Q_{0.9}^{\mathrm{visible}}(s_i).
$$

场景汇总先做颜色软权重。$\Delta E\_j$ 未测量时该项不进入求和，也不记为 0。$N$ 是计划局部数，实验中为 45。

$$
g_j=\mathrm{clip}\!\left(\frac{\Delta E_j}{0.02},0,1\right),\qquad
w_j=\frac{1}{N}\big[(1-\lambda)+\lambda g_j\big].
$$

把有效局部值按从小到大排序。加权中位数是累积权重第一次达到总权重一半时的那个数值：

$$
R_{\cdot,q}^{\mathrm{detail}}
=\min\lbrace v:\sum_{R_{\cdot,j,q}^{\mathrm{detail}}\le v} w_j
\ge \frac{1}{2}\sum_j w_j\rbrace.
$$

物理项在汇总之后、按分位数分别取较大值：

$$
R_{\mathrm{phy},50}^{\mathrm{detail}}
=\max\big(R_{\mathrm{obs},50}^{\mathrm{detail}},R_{\mathrm{kernel},50}^{\mathrm{detail}},R_{\mathrm{spacing},50}^{\mathrm{detail}}\big),
$$

$$
R_{\mathrm{phy},90}^{\mathrm{detail}}
=\max\big(R_{\mathrm{obs},90}^{\mathrm{detail}},R_{\mathrm{kernel},90}^{\mathrm{detail}},R_{\mathrm{spacing},90}^{\mathrm{detail}}\big).
$$

可见权重、$\Delta E$ 和某一局部是否有效由调用方传入。本库不渲染。

## 观测项

$R\_{\mathrm{obs}}^{\mathrm{detail}}$ 是细核组所在表面位置上的训练图像采样支持。训练相机不随层级增加，所以它不是新增观测分辨率。

在固定切平面的点 $x$ 上，相机 $c$ 的投影雅可比为 $J\_c(x)$，最弱方向的世界尺度是 $r\_c(x)=1/\sigma\_{\min}(J\_c(x))$。在视锥内可见、且夹角至少 $12^\circ$ 的真实训练相机对中，

$$
r_{\mathrm{obs}}(x)=\min_{(a,b)}\max(r_a(x),r_b(x)).
$$

细核中心先投影到冻结平面再查询。局部值仍是可见权重逆 CDF，场景汇总仍用上面的 $g\_j,w\_j$。没有有效相机对就留空。只做视锥筛选时，结果必须标明遮挡尚未验证。

`detail_physical` 就是上面的三分量最大值。Mip-NeRF360 原六层、$\lambda=0.5$ 的 60 个有效局部里，观测 P50/P90 都小于核和间距，所以这一轮物理项的数值由核或间距决定。Kitchen 观测 P50 停在 0.00304 附近；Room 约 0.0068 到 0.0092；Garden 约 0.0059 到 0.0151。这些起伏来自细核组的位置和权重，不是相机变化。

## 测试结果

下面是 Mip-NeRF360 的 Kitchen、Room、Garden 原六层 L0–L5。主结果取 $\lambda=0.5$，各层对用自身有效集合做逆 CDF 加权中位数。分母固定为 45。单位是各场景自己的坐标长度，场景之间不能比较。Kitchen 五层的有效组合是同一批；Room 五层共同集合只有 1 项；Garden 五层交集为空。

### 增量核

$$R_{\mathrm{kernel}}^{\mathrm{detail}}$$

| 层级对 | Kitchen P50 / P90 | 有效 | Room P50 / P90 | 有效 | Garden P50 / P90 | 有效 |
|---|---:|---:|---:|---:|---:|---:|
| L0→L1 | 0.06566 / 0.17142 | 6/45 | 0.37713 / 0.45366 | 1/45 | 0.19631 / 0.37186 | 1/45 |
| L1→L2 | 0.04698 / 0.06275 | 6/45 | 0.13201 / 1.13241 | 3/45 | 0.28689 / 0.75924 | 2/45 |
| L2→L3 | 0.03693 / 0.09274 | 6/45 | 0.09821 / 0.13502 | 3/45 | 0.14904 / 0.70756 | 3/45 |
| L3→L4 | 0.02812 / 0.06259 | 6/45 | 0.07964 / 0.14066 | 3/45 | 0.07034 / 0.15527 | 6/45 |
| L4→L5 | 0.02338 / 0.04481 | 7/45 | 0.07638 / 0.09561 | 3/45 | 0.03986 / 0.17116 | 4/45 |

### 增量间距

$$R_{\mathrm{spacing}}^{\mathrm{detail}}$$

| 层级对 | Kitchen P50 / P90 | 有效 | Room P50 / P90 | 有效 | Garden P50 / P90 | 有效 |
|---|---:|---:|---:|---:|---:|---:|
| L0→L1 | 0.03691 / 0.04609 | 6/45 | 0.22598 / 0.30125 | 1/45 | 0.22575 / 0.37167 | 1/45 |
| L1→L2 | 0.02657 / 0.03567 | 6/45 | 0.10210 / 0.36700 | 3/45 | 0.21483 / 0.52823 | 2/45 |
| L2→L3 | 0.01878 / 0.03355 | 6/45 | 0.05429 / 0.11470 | 3/45 | 0.10984 / 0.41474 | 3/45 |
| L3→L4 | 0.01291 / 0.02357 | 6/45 | 0.04479 / 0.09450 | 3/45 | 0.06178 / 0.16151 | 6/45 |
| L4→L5 | 0.01031 / 0.01779 | 7/45 | 0.03689 / 0.07846 | 3/45 | 0.02074 / 0.07041 | 4/45 |

### 物理项

$$R_{\mathrm{phy}}^{\mathrm{detail}}$$

P50 取观测、核、间距三项 P50 的较大值，P90 同样。这一轮观测更小，表中的数由核或间距决定。

| 层级对 | Kitchen P50 / P90 | Room P50 / P90 | Garden P50 / P90 |
|---|---:|---:|---:|
| L0→L1 | 0.06566 / 0.17142 | 0.37713 / 0.45366 | 0.22575 / 0.37186 |
| L1→L2 | 0.04698 / 0.06275 | 0.13201 / 1.13241 | 0.28689 / 0.75924 |
| L2→L3 | 0.03693 / 0.09274 | 0.09821 / 0.13502 | 0.14904 / 0.70756 |
| L3→L4 | 0.02812 / 0.06259 | 0.07964 / 0.14066 | 0.07034 / 0.16151 |
| L4→L5 | 0.02338 / 0.04481 | 0.07638 / 0.09561 | 0.03986 / 0.17116 |

Garden 的 L0→L1 由间距决定 P50（0.22575），由核决定 P90（0.37186）。L3→L4 由核决定 P50（0.07034），由间距决定 P90（0.16151）。其余层级对的 P50 和 P90 都由核决定。
