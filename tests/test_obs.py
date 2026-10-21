import math
import sys
from pathlib import Path

import numpy as np

from prf_detail.obs import _one_view, _project, detail_observation_quantiles, observation_resolutions


class Camera:
    def __init__(self, center, width=400, height=200, fx=100.0, fy=100.0, cx=200.0, cy=100.0):
        rotation = np.eye(3)
        self.w2c = np.eye(4)
        self.w2c[:3, :3] = rotation
        self.w2c[:3, 3] = -rotation @ np.asarray(center, dtype=np.float64)
        self.center = np.asarray(center, dtype=np.float64)
        self.view_id = "cam"
        self.width = width
        self.height = height
        self.fx = fx
        self.fy = fy
        self.cx = cx
        self.cy = cy

    def project(self, points):
        points = np.asarray(points, dtype=np.float64)
        single = points.ndim == 1
        uv, depth = _project(self, points.reshape(1, 3) if single else points)
        return (uv[0], depth[0]) if single else (uv, depth)


def test_one_pair_uses_the_worse_camera():
    point = np.array([[0.0, 0.0, 0.0]])
    tangent_1 = np.array([1.0, 0.0, 0.0])
    tangent_2 = np.array([0.0, 1.0, 0.0])
    near = Camera((0.0, 0.0, -4.0))
    far = Camera((2.0, 0.0, -8.0))
    values, visible, pairs = observation_resolutions(point, tangent_1, tangent_2, [near, far])
    assert visible[0] == 2 and pairs[0] == 1 and math.isfinite(values[0])
    assert math.isnan(observation_resolutions(point, tangent_1, tangent_2, [near])[0][0])
    near_value = _one_view(point, tangent_1, tangent_2, near, 0.05, 1.0, 1e-12)[0][0]
    far_value = _one_view(point, tangent_1, tangent_2, far, 0.05, 1.0, 1e-12)[0][0]
    assert math.isclose(values[0], max(near_value, far_value), rel_tol=0.0, abs_tol=1e-12)


def test_narrow_pair_stays_missing():
    point = np.array([[0.0, 0.0, 0.0]])
    tangent_1 = np.array([1.0, 0.0, 0.0])
    tangent_2 = np.array([0.0, 1.0, 0.0])
    # About 4 degrees apart at the origin, below the 12 degree gate.
    values, _, pairs = observation_resolutions(
        point, tangent_1, tangent_2, [Camera((0.0, 0.0, -8.0)), Camera((0.5, 0.0, -8.0))]
    )
    assert pairs[0] == 0 and math.isnan(values[0])


def test_visible_quantile_ignores_missing_points():
    points = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, -6.0]])
    result = detail_observation_quantiles(
        points,
        np.array([1.0, 0.0, 0.0]),
        np.array([0.0, 1.0, 0.0]),
        [Camera((0.0, 0.0, -4.0)), Camera((2.0, 0.0, -8.0))],
        weights=np.array([1.0, 1.0]),
    )
    assert result["count"] == 1
    assert result["p50"] == result["p90"]
    assert math.isfinite(result["p50"])


def test_matches_frozen_prf_observation():
    library = Path("/home/hongjiacheng/codes/worktrees/Skyfall-GS-prf-v1/thirdparty/PRF")
    if not library.exists():
        return
    sys.path.insert(0, str(library))
    from prf.config import PRFConfig
    from prf.obs import compute_observation_resolution
    from prf.types import SurfacePatch

    point = np.array([0.2, -0.1, 0.0])
    tangent_1 = np.array([1.0, 0.0, 0.0])
    tangent_2 = np.array([0.0, 1.0, 0.0])
    views = [Camera((0.0, 0.0, -4.0)), Camera((3.0, 1.0, -6.0)), Camera((-2.0, 2.0, -5.0))]
    patch = SurfacePatch(
        patch_id=0, center=point, normal=np.array([0.0, 0.0, 1.0]), t1=tangent_1, t2=tangent_2,
        tangent=np.column_stack((tangent_1, tangent_2)), ids=np.zeros(0, np.int64), radius=1.0, valid=True,
    )
    expected, _ = compute_observation_resolution(patch, views, PRFConfig())
    got = observation_resolutions(point.reshape(1, 3), tangent_1, tangent_2, views)[0][0]
    assert math.isclose(got, expected, rel_tol=0.0, abs_tol=1e-8)
