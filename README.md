# PRF v2

细核组上的三个增量指标。它们描述同一次层级加细里、方向尺度变小的那批高斯，不替代 v1 的全场景 patch 统计。

\[
R_{\mathrm{phy}}^{\mathrm{detail}}
=\max(R_{\mathrm{kernel}}^{\mathrm{detail}},R_{\mathrm{spacing}}^{\mathrm{detail}})
\]

P50 与 P90 分开取较大值。这个最大值作用在软加权汇总之后的两个数上。

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
    kernel_p50, kernel_p90, spacing_p50, spacing_p90,
    delta_e, eligible, lam=0.5, planned=45,
)
print(summary["kernel_p50"], summary["spacing_p50"], summary["phy_p50"])
```

## 冻结规则

- 展开且子/父方向尺度比 \(<0.99\) 为细核；\((0.99,1.01]\) 为相近；\(>1.01\) 为变粗；未展开为保留。尺度系数是 3.77。
- 组内 P50/P90 是可见权重的离散逆 CDF，`searchsorted` 取左侧。
- 间距是切平面上 6 个互不重合中心的平均距离。近重合邻居不占名额，邻居不足时留空，不减小 k。第六近邻碰到搜索边界时，该中心不进入分位数。
- 软权重 \(g=\mathrm{clip}(\Delta E/0.02,0,1)\)，\(w=(1/N)((1-\lambda)+\lambda g)\)。没有测到的 \(\Delta E\) 不补 0。\(N\) 由 `planned` 给出。
- 可见权重、\(\Delta E\) 和有效性由调用方传入。本库不渲染。

## 测试结果

下面是 Mip-NeRF360 的 Kitchen、Room、Garden 原六层 L0–L5。主结果取 \(\lambda=0.5\)，各层对用自身有效集合做逆 CDF 加权中位数。分母固定为 45。单位是各场景自己的坐标长度，场景之间不能比较。Kitchen 五层的有效组合是同一批；Room 五层共同集合只有 1 项；Garden 五层交集为空。

### \(R_{\mathrm{kernel}}^{\mathrm{detail}}\)

| 层级对 | Kitchen P50 / P90 | 有效 | Room P50 / P90 | 有效 | Garden P50 / P90 | 有效 |
|---|---:|---:|---:|---:|---:|---:|
| L0→L1 | 0.06566 / 0.17142 | 6/45 | 0.37713 / 0.45366 | 1/45 | 0.19631 / 0.37186 | 1/45 |
| L1→L2 | 0.04698 / 0.06275 | 6/45 | 0.13201 / 1.13241 | 3/45 | 0.28689 / 0.75924 | 2/45 |
| L2→L3 | 0.03693 / 0.09274 | 6/45 | 0.09821 / 0.13502 | 3/45 | 0.14904 / 0.70756 | 3/45 |
| L3→L4 | 0.02812 / 0.06259 | 6/45 | 0.07964 / 0.14066 | 3/45 | 0.07034 / 0.15527 | 6/45 |
| L4→L5 | 0.02338 / 0.04481 | 7/45 | 0.07638 / 0.09561 | 3/45 | 0.03986 / 0.17116 | 4/45 |

### \(R_{\mathrm{spacing}}^{\mathrm{detail}}\)

| 层级对 | Kitchen P50 / P90 | 有效 | Room P50 / P90 | 有效 | Garden P50 / P90 | 有效 |
|---|---:|---:|---:|---:|---:|---:|
| L0→L1 | 0.03691 / 0.04609 | 6/45 | 0.22598 / 0.30125 | 1/45 | 0.22575 / 0.37167 | 1/45 |
| L1→L2 | 0.02657 / 0.03567 | 6/45 | 0.10210 / 0.36700 | 3/45 | 0.21483 / 0.52823 | 2/45 |
| L2→L3 | 0.01878 / 0.03355 | 6/45 | 0.05429 / 0.11470 | 3/45 | 0.10984 / 0.41474 | 3/45 |
| L3→L4 | 0.01291 / 0.02357 | 6/45 | 0.04479 / 0.09450 | 3/45 | 0.06178 / 0.16151 | 6/45 |
| L4→L5 | 0.01031 / 0.01779 | 7/45 | 0.03689 / 0.07846 | 3/45 | 0.02074 / 0.07041 | 4/45 |

### \(R_{\mathrm{phy}}^{\mathrm{detail}}\)

P50 取核与间距 P50 的较大值，P90 同样。

| 层级对 | Kitchen P50 / P90 | Room P50 / P90 | Garden P50 / P90 |
|---|---:|---:|---:|
| L0→L1 | 0.06566 / 0.17142 | 0.37713 / 0.45366 | 0.22575 / 0.37186 |
| L1→L2 | 0.04698 / 0.06275 | 0.13201 / 1.13241 | 0.28689 / 0.75924 |
| L2→L3 | 0.03693 / 0.09274 | 0.09821 / 0.13502 | 0.14904 / 0.70756 |
| L3→L4 | 0.02812 / 0.06259 | 0.07964 / 0.14066 | 0.07034 / 0.16151 |
| L4→L5 | 0.02338 / 0.04481 | 0.07638 / 0.09561 | 0.03986 / 0.17116 |

Garden 的 L0→L1 由间距决定 P50（0.22575），由核决定 P90（0.37186）。L3→L4 由核决定 P50（0.07034），由间距决定 P90（0.16151）。其余层级对的 P50 和 P90 都由核决定。
