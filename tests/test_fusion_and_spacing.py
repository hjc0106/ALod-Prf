"""Checks that the copied v1 formulas still satisfy their definitions."""
import math

import numpy as np

from prf.config import PRFConfig
from prf.fusion import fuse_resolutions
from prf.kernel import compute_kernel_resolution
from prf.spacing import compute_spacing_resolution
from prf.types import GaussianScene, SurfacePatch


def _patch(center, ids):
    t1 = np.array([1.0, 0.0, 0.0])
    t2 = np.array([0.0, 1.0, 0.0])
    return SurfacePatch(
        patch_id=0,
        center=np.asarray(center, float),
        normal=np.array([0.0, 0.0, 1.0]),
        t1=t1,
        t2=t2,
        tangent=np.column_stack((t1, t2)),
        ids=np.asarray(ids, np.int64),
        radius=1.0,
        valid=True,
        area_m2=1.0,
        geometry_valid=True,
    )


def test_phys_is_the_max():
    patch = _patch([0, 0, 0], [0])
    record = fuse_resolutions(patch, 0.2, 0.5, 0.3, {"quality_flags": [], "num_visible_real_views": 2, "best_camera_pair": ("a", "b"), "best_pair_angle_deg": 20})
    assert record.R_phys_m_per_equiv_pixel == 0.5
    assert record.bottleneck_type == "KERNEL"


def test_spacing_mean_neighbor_distance():
    scene = GaussianScene(
        mu=np.array([[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]]),
        scales=np.ones((2, 3)),
        rotations=np.tile([1.0, 0.0, 0.0, 0.0], (2, 1)),
        opacity=np.ones(2),
    )
    value = compute_spacing_resolution(_patch([0, 0, 0], [0, 1]), scene, PRFConfig())
    assert math.isclose(value, 3.0, rel_tol=0, abs_tol=1e-9)


def test_kernel_uses_bandwidth_and_interpolated_p90():
    scene = GaussianScene(
        mu=np.zeros((1, 3)),
        scales=np.array([[2.0, 0.5, 0.1]]),
        rotations=np.array([[1.0, 0.0, 0.0, 0.0]]),
        opacity=np.array([1.0]),
    )
    value = compute_kernel_resolution(_patch([0, 0, 0], [0]), scene, PRFConfig())
    assert math.isclose(value, 3.77 * 2.0, rel_tol=0, abs_tol=1e-9)
