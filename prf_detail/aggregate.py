"""Soft-weighted scene summary and R_phy^detail.

Color gain is g = clip(deltaE / 0.02, 0, 1). Each planned local has base weight
1/planned. Unmeasured deltaE is left out rather than treated as zero.
R_phy P50 is max(obs P50, kernel P50, spacing P50); P90 is the same comparison.
"""
from __future__ import annotations

import math
from fractions import Fraction

from prf_detail.quantile import inverse_cdf_median

COLOR_SCALE = Fraction(1, 50)


def _gain(delta_e: str | float | Fraction | None) -> Fraction | None:
    if delta_e is None or delta_e == "":
        return None
    value = Fraction(delta_e) if isinstance(delta_e, str) else Fraction(delta_e)
    return max(Fraction(0), min(Fraction(1), value / COLOR_SCALE))


def aggregate_quantiles(values, delta_e, eligible, *, lam: float = 0.5, planned: int | None = None) -> float | None:
    flags = list(eligible)
    gains = [_gain(item) for item in delta_e]
    count = planned if planned is not None else len(flags)
    if count <= 0:
        return None
    weight_unit = Fraction(1, count)
    mixture = Fraction(lam).limit_denominator()
    items = []
    for value, gain, keep in zip(values, gains, flags):
        if not keep or gain is None or value is None or not math.isfinite(float(value)):
            continue
        weight = weight_unit * ((1 - mixture) + mixture * gain)
        if weight > 0:
            items.append((float(value), weight))
    return inverse_cdf_median(items)


def detail_physical(obs_p50, obs_p90, kernel_p50, kernel_p90, spacing_p50, spacing_p90) -> dict[str, float | None]:
    """R_phy^detail = max(R_obs^detail, R_kernel^detail, R_spacing^detail), per quantile."""

    def larger(*values):
        if any(value is None for value in values):
            return None
        return float(max(values))

    return {
        "p50": larger(obs_p50, kernel_p50, spacing_p50),
        "p90": larger(obs_p90, kernel_p90, spacing_p90),
    }


def detail_physical_three(obs_p50, obs_p90, kernel_p50, kernel_p90, spacing_p50, spacing_p90) -> dict[str, float | None]:
    """Alias of detail_physical. Kept so earlier three-component calls stay valid."""
    return detail_physical(obs_p50, obs_p90, kernel_p50, kernel_p90, spacing_p50, spacing_p90)


def summarize_detail(obs_p50, obs_p90, kernel_p50, kernel_p90, spacing_p50, spacing_p90, delta_e, eligible, *, lam: float = 0.5, planned: int | None = None) -> dict[str, float | None]:
    o50 = aggregate_quantiles(obs_p50, delta_e, eligible, lam=lam, planned=planned)
    o90 = aggregate_quantiles(obs_p90, delta_e, eligible, lam=lam, planned=planned)
    k50 = aggregate_quantiles(kernel_p50, delta_e, eligible, lam=lam, planned=planned)
    k90 = aggregate_quantiles(kernel_p90, delta_e, eligible, lam=lam, planned=planned)
    s50 = aggregate_quantiles(spacing_p50, delta_e, eligible, lam=lam, planned=planned)
    s90 = aggregate_quantiles(spacing_p90, delta_e, eligible, lam=lam, planned=planned)
    physical = detail_physical(o50, o90, k50, k90, s50, s90)
    return {
        "obs_p50": o50,
        "obs_p90": o90,
        "kernel_p50": k50,
        "kernel_p90": k90,
        "spacing_p50": s50,
        "spacing_p90": s90,
        "phy_p50": physical["p50"],
        "phy_p90": physical["p90"],
    }
