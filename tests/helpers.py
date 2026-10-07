"""Helpers shared by several test modules (``tests`` is on the pytest path)."""

import numpy as np


def random_positions(rng, n, M=1.0, min_gap=1e-8, max_gap=1e4, pole_gap=1e-3):
    """``n`` positions, ``r - 2M`` log-uniform in ``[min_gap, max_gap] * M``, theta off the axis."""
    x = np.empty((n, 4))
    x[:, 0] = rng.normal(scale=100.0, size=n)
    x[:, 1] = 2.0 * M + M * 10.0 ** rng.uniform(np.log10(min_gap), np.log10(max_gap), size=n)
    x[:, 2] = rng.uniform(pole_gap, np.pi - pole_gap, size=n)
    x[:, 3] = rng.uniform(0.0, 2.0 * np.pi, size=n)
    return x
