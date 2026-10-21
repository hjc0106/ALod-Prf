import math

import numpy as np

from prf_detail import detail_kernel_quantiles, detail_physical, detail_physical_three, neighbor_spacing, refinement_group, summarize_detail


def test_finer_group_threshold():
    group = refinement_group([0.5, 1.0, 2.0, 0.1], [1, 1, 1, 1], [True, True, True, False])
    assert list(group) == [0, 1, 2, 3]


def test_kernel_quantile_uses_only_finer_nodes():
    result = detail_kernel_quantiles([1.0, 2.0, 9.0], [2.0, 2.0, 2.0], [True, True, True], [1, 1, 1])
    assert result["count"] == 1
    assert result["p50"] == 1.0
    assert result["p90"] == 1.0


def test_spacing_ignores_a_duplicate_neighbor():
    coords = np.array([[0.0, 0.0], [1.0, 0.0], [1.0, 0.0], [0.0, 1.0], [2.0, 0.0], [-1.0, 0.0], [0.0, -1.0], [3.0, 0.0]])
    means, accepted = neighbor_spacing(coords, coords[:1], np.array([0]), k=6, eps=1e-6, margin=np.array([10.0]))
    assert accepted[0]
    assert math.isclose(means[0], (1 + 1 + 1 + 1 + 2 + 3) / 6, abs_tol=1e-9)


def test_physical_is_componentwise_max_after_aggregation():
    summary = summarize_detail(
        obs_p50=[0.01, 0.02],
        obs_p90=[0.03, 0.04],
        kernel_p50=[0.2, 0.4],
        kernel_p90=[0.5, 0.6],
        spacing_p50=[0.3, 0.1],
        spacing_p90=[0.35, 0.9],
        delta_e=["0.04", "0.0"],
        eligible=[True, True],
        lam=0.0,
        planned=2,
    )
    assert summary["kernel_p50"] == 0.2
    assert summary["spacing_p90"] == 0.35
    physical = detail_physical(summary["obs_p50"], summary["obs_p90"], summary["kernel_p50"], summary["kernel_p90"], summary["spacing_p50"], summary["spacing_p90"])
    assert physical["p50"] == max(summary["obs_p50"], summary["kernel_p50"], summary["spacing_p50"])
    assert physical["p90"] == max(summary["obs_p90"], summary["kernel_p90"], summary["spacing_p90"])
    assert summary["phy_p50"] == physical["p50"]
    assert detail_physical_three(1.0, 1.0, 0.2, 0.3, 0.1, 0.2) == {"p50": 1.0, "p90": 1.0}
