"""The Schwarzschild metric in Schwarzschild (``a = 0`` Boyer-Lindquist) coordinates.

The kernels are written by hand from the tables in ``docs/theory/formulations.md`` (eqs. 2, 15
and §3.1) and checked against ``tools/derive_metric.py`` in the tests. They are written in terms
of ``r - 2M`` instead of ``f = 1 - 2M/r``: for ``M <= r <= 4M`` that subtraction is exact
(Sterbenz lemma), so every component keeps full relative accuracy down to the horizon, where
``f`` itself would carry an absolute error of one ulp of 1.
"""

from __future__ import annotations

import functools

import numpy as np
from numba import njit

from .base import Metric, MetricKernels


class Schwarzschild(Metric):
    """Schwarzschild black hole of mass ``M``, ``x = (t, r, theta, phi)``.

    Constants of motion (theory §5.1): energy ``E = -p_t``, axial angular momentum
    ``Lz = p_phi`` and total angular momentum squared ``L2 = p_theta^2 + p_phi^2 / sin^2 theta``.
    """

    invariant_names = ("E", "Lz", "L2")
    cyclic_coords = (0, 3)

    def __init__(self, M: float = 1.0):
        if not M > 0:
            raise ValueError(f"the mass must be positive, got M = {M}")
        self.M = float(M)
        super().__init__(_kernels(self.M))

    @property
    def horizon(self) -> float:
        return 2.0 * self.M

    def __repr__(self) -> str:
        return f"Schwarzschild(M={self.M})"


@functools.cache
def _kernels(M: float) -> MetricKernels:
    """Compile the kernels once per mass; ``M`` becomes a compile-time constant (ADR 0002)."""
    two_M = 2.0 * M

    @njit
    def g(x):
        r, s = x[1], np.sin(x[2])
        out = np.zeros((4, 4))
        out[0, 0] = -(r - two_M) / r
        out[1, 1] = r / (r - two_M)
        out[2, 2] = r * r
        out[3, 3] = (r * s) ** 2
        return out

    @njit
    def g_inv(x):
        r, s = x[1], np.sin(x[2])
        out = np.zeros((4, 4))
        out[0, 0] = -r / (r - two_M)
        out[1, 1] = (r - two_M) / r
        out[2, 2] = 1.0 / (r * r)
        out[3, 3] = 1.0 / (r * s) ** 2
        return out

    @njit
    def dg_inv(x):
        r, s, c = x[1], np.sin(x[2]), np.cos(x[2])
        d = r - two_M
        out = np.zeros((4, 4, 4))
        out[1, 0, 0] = two_M / (d * d)
        out[1, 1, 1] = two_M / (r * r)
        out[1, 2, 2] = -2.0 / r**3
        out[1, 3, 3] = -2.0 / (r**3 * s * s)
        out[2, 3, 3] = -2.0 * c / (r * r * s**3)
        return out

    @njit
    def christoffel(x):
        r, s, c = x[1], np.sin(x[2]), np.cos(x[2])
        d = r - two_M
        out = np.zeros((4, 4, 4))
        out[0, 0, 1] = out[0, 1, 0] = M / (r * d)
        out[1, 0, 0] = M * d / r**3
        out[1, 1, 1] = -M / (r * d)
        out[1, 2, 2] = -d
        out[1, 3, 3] = -d * s * s
        out[2, 1, 2] = out[2, 2, 1] = 1.0 / r
        out[2, 3, 3] = -s * c
        out[3, 1, 3] = out[3, 3, 1] = 1.0 / r
        out[3, 2, 3] = out[3, 3, 2] = c / s
        return out

    @njit
    def invariants(x, p):
        s = np.sin(x[2])
        out = np.empty(3)
        out[0] = -p[0]
        out[1] = p[3]
        out[2] = p[2] * p[2] + (p[3] / s) ** 2
        return out

    return MetricKernels(g, g_inv, dg_inv, christoffel, invariants)
