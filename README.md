# MTF v1

纹理传递的 MTF 计算。渲染器不在这个库里。调用方给出每个周期的余弦输入图和 90° 相移输入图，或者直接给出已经拟合的四个谐波系数。

曲线值是传递矩阵两个奇异值的均方根。MTF25 / MTF50 是增益第一次向上穿过 0.25 / 0.50 的周期，周期在对数坐标上插值。向下跌落超过 0.02 时记入状态；没有交点时不给出有限周期。

## 安装

```bash
pip install numpy
pip install -e .
```

## 调用

```python
from mtf import gain_from_images, gain_rms, mtf_periods

gains = []
for period, phase0, image0, phase90, image90 in measurements:
    gains.append(gain_from_images(phase0, image0, phase90, image90))

result = mtf_periods(periods, gains)
print(result["MTF25"]["period"], result["MTF25"]["status"])
print(result["MTF50"]["period"], result["MTF50"]["status"])

# 已有四个系数时：相位 0 的 (cos, sin)，相位 90 的 (cos, sin)
print(gain_rms(c0, s0, c90, s90))
```

相位 90 的输入是 \(\cos(\phi+\pi/2)=-\sin\phi\)。矩阵把这一列取负，再放到正正弦基上。
