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
