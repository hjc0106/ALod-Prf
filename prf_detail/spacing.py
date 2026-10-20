"""2D surface spacing of detail centers.

For each query, s is the mean distance to 6 other centers. A neighbor within
``eps`` of the query, or within ``eps`` of an already accepted neighbor, does
not fill a slot. If the sixth distance reaches ``margin``, the value is kept
for diagnosis and marked not accepted. k is never reduced.
"""
from __future__ import annotations

import numpy as np
from scipy.spatial import cKDTree

from prf_detail.quantile import inverse_cdf_quantile


def neighbor_spacing(coords, queries, own_index, *, k: int = 6, eps: float, margin=None, exclude_near_query: bool = True):
    coords = np.asarray(coords, dtype=np.float64)
    queries = np.asarray(queries, dtype=np.float64)
    own_index = np.asarray(own_index, dtype=np.int64)
    count, rows = len(coords), len(queries)
    distances = np.full((rows, k), np.nan)
    if count == 0 or rows == 0:
        accepted = np.zeros(rows, dtype=bool)
        return distances.mean(axis=1), accepted
    tree = cKDTree(coords)
    pending = np.arange(rows)
    width = min(count, max(32, k + 1))
    eps2 = eps * eps
    while pending.size:
        found_d, found_i = tree.query(queries[pending], k=width, workers=1)
        found_d = np.asarray(found_d, dtype=np.float64)
        found_i = np.asarray(found_i, dtype=np.int64)
        if found_d.ndim < 2:
            found_d = found_d.reshape(len(pending), width)
            found_i = found_i.reshape(len(pending), width)
        retry = []
        for row, query_index in enumerate(pending):
            own = int(own_index[query_index])
            chosen = []
            for dist, index in zip(found_d[row], found_i[row]):
                index = int(index)
                if index < 0 or not np.isfinite(dist) or index == own:
                    continue
                if exclude_near_query and dist <= eps:
                    continue
                if any(float(np.dot(delta, delta)) <= eps2 for delta in (coords[index] - coords[item] for item in chosen)):
                    continue
                chosen.append(index)
                if len(chosen) == k:
                    break
            if len(chosen) == k:
                distances[query_index] = [float(np.linalg.norm(coords[index] - queries[query_index])) for index in chosen]
            elif width != count:
                retry.append(query_index)
        if not retry or width == count:
            break
        pending = np.asarray(retry, dtype=np.int64)
        width = min(count, width * 2)
    means = np.where(np.isfinite(distances).all(axis=1), distances.mean(axis=1), np.nan)
    if margin is None:
        accepted = np.isfinite(means)
    else:
        margin = np.asarray(margin, dtype=np.float64)
        accepted = np.isfinite(means) & (distances[:, -1] < margin)
    return means, accepted


def detail_spacing_quantiles(spacing, accepted, weights) -> dict[str, float | int | None]:
    values = np.where(np.asarray(accepted, dtype=bool), np.asarray(spacing, dtype=np.float64), np.nan)
    weights = np.asarray(weights, dtype=np.float64)
    return {
        "count": int(np.isfinite(values).sum()),
        "p50": inverse_cdf_quantile(values, weights, 0.5),
        "p90": inverse_cdf_quantile(values, weights, 0.9),
    }
