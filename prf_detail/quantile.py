"""Discrete inverse-CDF quantiles used by the detail metrics.

This is not the interpolated CDF used by PRF v1 kernel weights.
"""
from __future__ import annotations

from fractions import Fraction

import numpy as np


def inverse_cdf_quantile(values, weights, quantile: float) -> float | None:
    values = np.asarray(values, dtype=np.float64)
    weights = np.asarray(weights, dtype=np.float64)
    finite = np.isfinite(values) & np.isfinite(weights) & (weights >= 0)
    values, weights = values[finite], weights[finite]
    if values.size == 0 or float(weights.sum()) <= 0.0:
        return None
    order = np.argsort(values)
    total = float(weights.sum())
    position = int(np.searchsorted(np.cumsum(weights[order]), quantile * total, side="left"))
    return float(values[order[min(position, len(order) - 1)]])


def inverse_cdf_median(items: list[tuple[float, Fraction]]) -> float | None:
    if not items:
        return None
    total = sum((weight for _, weight in items), Fraction())
    cumulative = Fraction()
    for value, weight in sorted(items, key=lambda item: item[0]):
        cumulative += weight
        if 2 * cumulative >= total:
            return float(value)
    raise AssertionError("median did not reach half the weight")
