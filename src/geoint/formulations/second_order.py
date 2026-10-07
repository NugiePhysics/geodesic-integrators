"""Formulation (a): the geodesic equation as a first-order system in ``y = (x^mu, u^mu)``.

    dx^mu/dlambda = u^mu,    du^mu/dlambda = -Gamma^mu_{alpha beta}(x) u^alpha u^beta

(theory eq. 7). The contraction runs over all 4 x 4 index pairs, so it is correct for any metric.
"""

from __future__ import annotations

import functools

import numpy as np
from numba import njit
from numpy.typing import ArrayLike, NDArray

from ..metrics import Metric, MetricKernels
from .base import Formulation, FormulationKernels, as_array


class SecondOrder(Formulation):
    """Geodesic equation in ``(x, u)`` with ``u = dx/dlambda``.

    The mass shell is ``H = g_{mu nu} u^mu u^nu / 2``, which equals the Hamiltonian at
    ``p = g u``. The invariants ``E``, ``Lz``, ... are nonlinear functions of ``(x, u)`` here, so
    their conservation is a genuine test of the integrator (theory §5.2).
    """

    name = "second_order"

    def __init__(self, metric: Metric):
        super().__init__(metric, _kernels(metric.kernels))

    def from_xp(self, x: ArrayLike, p: ArrayLike) -> NDArray[np.float64]:
        x, p = as_array(x, (4,), "x"), as_array(p, (4,), "p")
        return np.concatenate([x, self.metric.g_inv(x) @ p])

    def to_xp(self, y: ArrayLike) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        y = as_array(y, (self.dim,), "y")
        x = y[:4].copy()
        return x, self.metric.g(x) @ y[4:]


@functools.cache
def _kernels(metric: MetricKernels) -> FormulationKernels:
    christoffel, g, invariants = metric.christoffel, metric.g, metric.invariants

    @njit
    def rhs(y):
        x, u = y[:4], y[4:]
        gamma = christoffel(x)
        dy = np.empty(8)
        for mu in range(4):
            acc = 0.0
            for a in range(4):
                for b in range(4):
                    acc += gamma[mu, a, b] * u[a] * u[b]
            dy[mu] = u[mu]
            dy[4 + mu] = -acc
        return dy

    @njit
    def point_invariants(y):
        x, u = y[:4], y[4:]
        gx = g(x)
        p = np.zeros(4)
        for mu in range(4):
            for nu in range(4):
                p[mu] += gx[mu, nu] * u[nu]
        H = 0.0
        for mu in range(4):
            H += p[mu] * u[mu]
        rest = invariants(x, p)
        out = np.empty(1 + rest.size)
        out[0] = 0.5 * H
        out[1:] = rest
        return out

    return FormulationKernels(rhs, point_invariants)
