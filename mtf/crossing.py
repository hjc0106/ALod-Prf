"""First upward threshold crossing on a period/gain curve.

Periods are interpolated in log space. A downward step larger than
``nonmonotone_drop`` is recorded on the status and, when no crossing exists,
suppresses a finite period.
"""
from __future__ import annotations

import math
from typing import Mapping, Sequence


def _max_downward_step(rows: Sequence[Mapping[str, float]], key: str) -> float:
    valid = [row for row in rows if math.isfinite(float(row[key]))]
    valid.sort(key=lambda row: float(row["period"]))
    drops = [float(first[key]) - float(second[key]) for first, second in zip(valid[:-1], valid[1:])]
    return max(0.0, max(drops, default=0.0))


def threshold_crossing(
    periods: Sequence[float],
    gains: Sequence[float],
    threshold: float,
    nonmonotone_drop: float = 0.02,
) -> dict[str, float | str]:
    rows = [
        {"period": float(period), "gain": float(gain)}
        for period, gain in zip(periods, gains)
        if math.isfinite(float(period)) and math.isfinite(float(gain))
    ]
    rows.sort(key=lambda row: row["period"])
    if len(rows) < 2:
        return {"period": math.nan, "low": math.nan, "high": math.nan, "status": "insufficient_valid_periods", "drop": 0.0}
    drop = _max_downward_step(rows, "gain")
    suffix = ""
    if drop > nonmonotone_drop:
        token = f"{nonmonotone_drop:g}".replace(".", "p")
        suffix = f";nonmonotone_drop_gt_{token}"
    for first, second in zip(rows[:-1], rows[1:]):
        y0, y1 = first["gain"], second["gain"]
        if y0 < threshold <= y1:
            low, high = first["period"], second["period"]
            estimate = math.exp(
                math.log(low) + (threshold - y0) * (math.log(high) - math.log(low)) / max(y1 - y0, 1e-12)
            )
            return {"period": estimate, "low": low, "high": high, "status": "first_upward_crossing" + suffix, "drop": drop}
    if drop > nonmonotone_drop:
        status = "nonmonotone_no_threshold_crossing"
    elif rows[0]["gain"] >= threshold:
        status = "at_or_below_tested_minimum"
    elif rows[-1]["gain"] < threshold:
        status = "above_tested_maximum"
    else:
        status = "uncrossed"
    return {"period": math.nan, "low": math.nan, "high": math.nan, "status": status, "drop": drop}


def mtf_periods(
    periods: Sequence[float],
    gains: Sequence[float],
    thresholds: Sequence[float] = (0.25, 0.50),
    nonmonotone_drop: float = 0.02,
) -> dict[str, dict[str, float | str]]:
    """MTF25 and MTF50 by default. ``thresholds`` are gain fractions, not percentages."""
    return {
        f"MTF{int(round(100 * threshold))}": threshold_crossing(periods, gains, threshold, nonmonotone_drop)
        for threshold in thresholds
    }
