# PRF v1

旧物理分辨率场的计算库。输入是同一坐标系下的三维高斯和训练相机，输出每个表面 patch 的四项分辨率。

\[
R_{\mathrm{phys}}=\max(R_{\mathrm{obs}},R_{\mathrm{kernel}},R_{\mathrm{spacing}})
\]

系数与分位数沿用已冻结的实现，本库不改公式。MTF 在 `mtf-v1` 分支，细核组上的增量核、增量间距和 \(R_{\mathrm{phy}}^{\mathrm{detail}}\) 在 `prf-v2` 分支。

## 安装

```bash
pip install numpy scipy
pip install -e .
```

遮挡深度是可选项，需要 PyTorch、CUDA 和 gsplat 1.4.0：

```bash
pip install -r requirements-depth.txt
```

## 调用

```python
from prf import PRFConfig, build_physical_resolution_field
from prf.cameras import load_training_views
from prf.io_gs import load_gaussians_ply

scene = load_gaussians_ply("point_cloud.ply")
views = load_training_views("transforms_train.json")
field = build_physical_resolution_field(scene, views, PRFConfig())
for record in field.records:
    print(record.R_obs_m_per_equiv_pixel,
          record.R_kernel_m_per_equiv_pixel,
          record.R_spacing_m_per_equiv_pixel,
          record.R_phys_m_per_equiv_pixel,
          record.bottleneck_type)
```

命令行：

```bash
python -m prf --ply point_cloud.ply --transforms transforms_train.json --output out/frustum
```

PLY 为 binary little-endian 3DGS：scale 是对数，opacity 是 logit，四元数是 wxyz。相机默认按 OpenCV/COLMAP 的 C2W 读取。字段名里的 `_m` 是历史别名，不表示已经确认米制尺度。

## 三项定义

| 分量 | 规则 |
|---|---|
| `R_obs` | 切平面投影雅可比最小奇异值的倒数。只保留可见相机，再在夹角至少 12° 的相机对里取较差一侧的最小值。默认可见性是视锥。 |
| `R_kernel` | 切平面协方差最大标准差乘 3.77，按透明度平方和投影面积加权，取插值 CDF 的 90% 分位。 |
| `R_spacing` | 切平面上 6 近邻平均距离的 90% 分位。 |
| `R_phys` | 三项的最大值。有限项中的最大者记为瓶颈 `OBS` / `KERNEL` / `SPACING`。 |

表面块来自体素 anchor、法向聚类、空间连通拆分和去重。细节见 `DESIGN.md` 与 `SCHEMA_V2.md`。
