import math

import numpy as np

from mtf import gain_from_images, gain_rms, mtf_periods


def test_identity_transfer_has_unit_rms():
    assert math.isclose(gain_rms(1.0, 0.0, 0.0, -1.0), 1.0, abs_tol=1e-12)


def test_image_fit_recovers_half_gain():
    phase = np.linspace(0.0, 4.0 * math.pi, 64, endpoint=False)
    image0 = 0.5 * np.cos(phase)
    image90 = 0.5 * np.cos(phase + math.pi / 2.0)
    assert math.isclose(gain_from_images(phase, image0, phase, image90), 0.5, abs_tol=1e-9)


def test_log_crossing():
    periods = np.geomspace(1.0, 8.0, 4)
    gains = [0.1, 0.2, 0.4, 0.8]
    # 0.25 lies between 0.2 at period periods[1] and 0.4 at periods[2]
    low, high = float(periods[1]), float(periods[2])
    expected = math.exp(math.log(low) + (0.25 - 0.2) * (math.log(high) - math.log(low)) / 0.2)
    result = mtf_periods(periods, gains)
    assert result["MTF25"]["status"] == "first_upward_crossing"
    assert math.isclose(float(result["MTF25"]["period"]), expected, abs_tol=1e-12)
    assert result["MTF50"]["status"] == "first_upward_crossing"
