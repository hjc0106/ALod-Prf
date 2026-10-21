"""Training-view observation support on a fixed tangent plane.

Visibility is the camera frustum only. Occlusion is not tested. The finite
difference step defaults to 0.05 scene-length units, matching a unit patch
radius under the frozen Jacobian ratio and clamp.
"""
from __future__ import annotations

import numpy as np

from prf_detail.quantile import inverse_cdf_quantile

DEFAULT_DELTA = 0.05
DEFAULT_MIN_CAMERA_Z = 1.0
DEFAULT_MIN_PAIR_ANGLE_DEG = 12.0
DEFAULT_EPS = 1e-12


def _project(view, points: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    ones = np.ones((len(points), 1), dtype=np.float64)
    camera = (np.asarray(view.w2c, dtype=np.float64) @ np.concatenate([points, ones], axis=1).T).T
    depth = camera[:, 2]
    safe = np.where(np.abs(depth) > 1e-12, depth, 1.0)
    uv = np.empty((len(points), 2), dtype=np.float64)
    uv[:, 0] = float(view.fx) * camera[:, 0] / safe + float(view.cx)
    uv[:, 1] = float(view.fy) * camera[:, 1] / safe + float(view.cy)
    return uv, depth


def _one_view(points, tangent_1, tangent_2, view, delta, min_camera_z, eps):
    origin, depth0 = _project(view, points)
    plus_u, depth_u = _project(view, points + delta * tangent_1)
    plus_v, depth_v = _project(view, points + delta * tangent_2)
    jacobian = np.stack(((plus_u - origin) / delta, (plus_v - origin) / delta), axis=-1)
    gram = np.einsum("nji,njk->nik", jacobian, jacobian)
    gamma_min = np.sqrt(np.maximum(np.linalg.eigvalsh(gram), 0.0))[:, 0]
    finite = (
        (depth0 > min_camera_z)
        & (depth_u > min_camera_z)
        & (depth_v > min_camera_z)
        & np.isfinite(origin).all(axis=1)
        & np.isfinite(plus_u).all(axis=1)
        & np.isfinite(plus_v).all(axis=1)
        & (gamma_min > eps)
    )
    inside = (
        (origin[:, 0] >= 0.0)
        & (origin[:, 0] < view.width)
        & (origin[:, 1] >= 0.0)
        & (origin[:, 1] < view.height)
    )
    visible = finite & inside
    resolution = np.full(len(points), np.inf)
    resolution[visible] = 1.0 / gamma_min[visible]
    return resolution, visible


def observation_resolutions(
    points,
    tangent_1,
    tangent_2,
    views,
    *,
    delta: float = DEFAULT_DELTA,
    min_camera_z: float = DEFAULT_MIN_CAMERA_Z,
    min_pair_angle_deg: float = DEFAULT_MIN_PAIR_ANGLE_DEG,
    eps: float = DEFAULT_EPS,
    batch: int = 64,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return r_obs, visible-camera counts, and valid-pair counts.

    ``points`` lie on the tangent plane. Each view needs ``w2c``, ``center``,
    ``fx``, ``fy``, ``cx``, ``cy``, ``width`` and ``height``. A point with no
    qualifying camera pair is NaN.
    """
    points = np.asarray(points, dtype=np.float64)
    tangent_1 = np.asarray(tangent_1, dtype=np.float64)
    tangent_2 = np.asarray(tangent_2, dtype=np.float64)
    if points.ndim != 2 or points.shape[1] != 3:
        raise ValueError("points must have shape (N, 3)")
    fields = [_one_view(points, tangent_1, tangent_2, view, delta, min_camera_z, eps) for view in views]
    resolution = np.stack([item[0] for item in fields]).T
    visible = np.stack([item[1] for item in fields]).T
    cameras = np.stack([np.asarray(view.center, dtype=np.float64) for view in views])
    upper = np.triu(np.ones((len(views), len(views)), dtype=bool), 1)
    answer = np.full(len(points), np.nan)
    counts = np.zeros(len(points), dtype=np.int32)
    for start in range(0, len(points), batch):
        stop = min(start + batch, len(points))
        direction = cameras[None, :, :] - points[start:stop, None, :]
        direction /= np.maximum(np.linalg.norm(direction, axis=2, keepdims=True), 1e-12)
        cosine = np.matmul(direction, direction.transpose(0, 2, 1))
        angles = np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0)))
        valid = (angles >= min_pair_angle_deg) & upper[None, :, :] & visible[start:stop, :, None] & visible[start:stop, None, :]
        count = valid.sum(axis=(1, 2))
        counts[start:stop] = count
        pair = np.maximum(resolution[start:stop, :, None], resolution[start:stop, None, :])
        pair = np.where(valid, pair, np.inf)
        value = pair.min(axis=(1, 2))
        answer[start:stop] = np.where(count > 0, value, np.nan)
    return answer, visible.sum(axis=1).astype(np.int32), counts


def detail_observation_quantiles(points, tangent_1, tangent_2, views, weights, **kwargs) -> dict[str, float | int | None | np.ndarray]:
    """Visible-weighted P50/P90 of r_obs at the supplied plane points."""
    values, visible_views, valid_pairs = observation_resolutions(points, tangent_1, tangent_2, views, **kwargs)
    weights = np.asarray(weights, dtype=np.float64)
    return {
        "values": values,
        "visible_views": visible_views,
        "valid_pairs": valid_pairs,
        "count": int(np.isfinite(values).sum()),
        "p50": inverse_cdf_quantile(values, weights, 0.5),
        "p90": inverse_cdf_quantile(values, weights, 0.9),
    }
