"""Frozen harmonic transfer used by the published MTF curve.

phase90 injects cos(phi + pi/2) = -sin(phi). Its response column is negated
so the 2x2 matrix is in the positive-sine basis. The curve value is the RMS
of the two singular values.
"""
from __future__ import annotations

import math

import numpy as np


def fit_harmonic(phase: np.ndarray, values: np.ndarray) -> tuple[float, float, float]:
    """Return DC, cosine and sine coefficients of values on phase."""
    phase = np.asarray(phase, dtype=np.float64).ravel()
    values = np.asarray(values, dtype=np.float64).ravel()
    valid = np.isfinite(phase) & np.isfinite(values)
    phase, values = phase[valid], values[valid]
    if phase.size < 3:
        return math.nan, math.nan, math.nan
    design = np.column_stack((np.ones(phase.size), np.cos(phase), np.sin(phase)))
    coeff, *_ = np.linalg.lstsq(design, values, rcond=None)
    return float(coeff[0]), float(coeff[1]), float(coeff[2])


def transfer_matrix(c0: float, s0: float, c90: float, s90: float) -> np.ndarray:
    return np.asarray([[c0, -c90], [s0, -s90]], dtype=np.float64)


def gain_rms(c0: float, s0: float, c90: float, s90: float) -> float:
    matrix = transfer_matrix(c0, s0, c90, s90)
    if not np.isfinite(matrix).all():
        return math.nan
    singular = np.linalg.svd(matrix, compute_uv=False)
    return float(math.sqrt(float(np.mean(np.square(singular)))))


def gain_from_images(
    phase0: np.ndarray,
    image0: np.ndarray,
    phase90: np.ndarray,
    image90: np.ndarray,
) -> float:
    """RMS gain from the cosine-input image and the phase-shifted input image."""
    _, c0, s0 = fit_harmonic(phase0, image0)
    _, c90, s90 = fit_harmonic(phase90, image90)
    return gain_rms(c0, s0, c90, s90)
