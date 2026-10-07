"""Formulation (b): Hamilton's equations for ``H = g^{mu nu}(x) p_mu p_nu / 2`` in ``(x^mu, p_mu)``.

    dx^mu/dlambda = dH/dp_mu = g^{mu nu} p_nu,
    dp_mu/dlambda = -dH/dx^mu = -(1/2) d_mu g^{alpha beta} p_alpha p_beta

(theory eq. 12). Only ``g^{mu nu}`` and its first derivatives enter; no Christoffel symbols.
"""

from __future__ import annotations

import functools
from typing import NamedTuple

import numpy as np
from numba import njit
from numpy.typing import ArrayLike, NDArray

from ..metrics import Metric, MetricKernels
from .base import Formulation, FormulationKernels, Kernel, as_array


class HamiltonianKernels(NamedTuple):
    """The kernels of every formulation plus the two halves of the gradient of ``H``."""

    rhs: Kernel
    point_invariants: Kernel
    dH_dx: Kernel
    dH_dp: Kernel


class Hamiltonian(Formulation):
    """Canonical equations in ``(x, p)``.

    Besides ``rhs`` this exposes the two halves of the gradient, ``dH_dx(x, p)`` and
    ``dH_dp(x, p)`` (each ``(4,) -> (4,)``), which Tao's extended-phase-space method needs as
    separate functions. For every cyclic coordinate of the metric the matching component of
    ``rhs`` is an exact floating-point zero, so ``p_t`` and ``p_phi`` stay constant bit for bit
    under any Runge-Kutta method (theory §5.2).
    """

    name = "hamiltonian"

    def __init__(self, metric: Metric):
        kernels = _kernels(metric.kernels)
        super().__init__(metric, FormulationKernels(kernels.rhs, kernels.point_invariants))
        self.dH_dx = kernels.dH_dx
        self.dH_dp = kernels.dH_dp

    def from_xp(self, x: ArrayLike, p: ArrayLike) -> NDArray[np.float64]:
        return np.concatenate([as_array(x, (4,), "x"), as_array(p, (4,), "p")])

    def to_xp(self, y: ArrayLike) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        y = as_array(y, (self.dim,), "y")
        return y[:4].copy(), y[4:].copy()


@functools.cache
def _kernels(metric: MetricKernels) -> HamiltonianKernels:
    g_inv, dg_inv, invariants = metric.g_inv, metric.dg_inv, metric.invariants

    @njit
    def dH_dp(x, p):
        gi = g_inv(x)
        out = np.zeros(4)
        for mu in range(4):
            for nu in range(4):
                out[mu] += gi[mu, nu] * p[nu]
        return out

    @njit
    def dH_dx(x, p):
        dgi = dg_inv(x)
        out = np.empty(4)
        for a in range(4):
            acc = 0.0
            for mu in range(4):
                for nu in range(4):
                    acc += dgi[a, mu, nu] * p[mu] * p[nu]
            out[a] = 0.5 * acc
        return out

    @njit
    def rhs(y):
        x, p = y[:4], y[4:]
        dy = np.empty(8)
        dy[:4] = dH_dp(x, p)
        dy[4:] = -dH_dx(x, p)
        return dy

    @njit
    def point_invariants(y):
        x, p = y[:4], y[4:]
        dx = dH_dp(x, p)
        H = 0.0
        for mu in range(4):
            H += p[mu] * dx[mu]
        rest = invariants(x, p)
        out = np.empty(1 + rest.size)
        out[0] = 0.5 * H
        out[1:] = rest
        return out

    return HamiltonianKernels(rhs, point_invariants, dH_dx, dH_dp)
