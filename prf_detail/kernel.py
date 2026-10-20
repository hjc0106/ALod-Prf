"""Directional incremental kernel on the finer group.

A child is finer when it replaces its parent and its directional scale is
below 0.99 times the parent scale. Scale is 3.77 times the tangent standard
deviation. P50 and P90 use the visible-weight inverse CDF.
"""
from __future__ import annotations

import numpy as np

from prf_detail.quantile import inverse_cdf_quantile

BANDWIDTH = 3.77
FINER_RATIO = 0.99
COARSER_RATIO = 1.01


def covariance_world(rotations: np.ndarray, scales: np.ndarray) -> np.ndarray:
    quaternions = np.asarray(rotations, dtype=np.float64)
    quaternions = quaternions / np.maximum(np.linalg.norm(quaternions, axis=1, keepdims=True), 1e-12)
    w, x, y, z = quaternions.T
    matrices = np.stack(
        (
            1 - 2 * (y * y + z * z),
            2 * (x * y - w * z),
            2 * (x * z + w * y),
            2 * (x * y + w * z),
            1 - 2 * (x * x + z * z),
            2 * (y * z - w * x),
            2 * (x * z - w * y),
            2 * (y * z + w * x),
            1 - 2 * (x * x + y * y),
        ),
        axis=1,
    ).reshape(-1, 3, 3)
    scaled = matrices * np.asarray(scales, dtype=np.float64)[:, None, :]
    return np.matmul(scaled, scaled.transpose(0, 2, 1))


def directional_tangent(tangent_1, tangent_2, direction_deg: float) -> np.ndarray:
    angle = np.deg2rad(direction_deg)
    return np.cos(angle) * np.asarray(tangent_1, dtype=np.float64) + np.sin(angle) * np.asarray(tangent_2, dtype=np.float64)


def directional_scale(covariances: np.ndarray, tangent: np.ndarray) -> np.ndarray:
    projected = np.einsum("i,nij,j->n", tangent, covariances, tangent)
    return BANDWIDTH * np.sqrt(np.maximum(projected, 0.0))


def refinement_group(child_scale, parent_scale, expanded) -> np.ndarray:
    """0 finer, 1 similar, 2 coarser, 3 retained."""
    ratio = np.asarray(child_scale, dtype=np.float64) / np.asarray(parent_scale, dtype=np.float64)
    expanded = np.asarray(expanded, dtype=bool)
    return np.where(~expanded, 3, np.where(ratio < FINER_RATIO, 0, np.where(ratio > COARSER_RATIO, 2, 1))).astype(np.int64)


def detail_kernel_quantiles(child_scale, parent_scale, expanded, weights) -> dict[str, float | int | None]:
    group = refinement_group(child_scale, parent_scale, expanded)
    selected = group == 0
    used_values = np.asarray(child_scale, dtype=np.float64)[selected]
    used_weights = np.asarray(weights, dtype=np.float64)[selected]
    return {
        "count": int(selected.sum()),
        "p50": inverse_cdf_quantile(used_values, used_weights, 0.5),
        "p90": inverse_cdf_quantile(used_values, used_weights, 0.9),
    }
