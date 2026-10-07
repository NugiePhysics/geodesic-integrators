"""TC0: toy problems that verify the integrators before any general relativity is involved.

All three are canonical Hamiltonian systems in ``y = (q, p)``:

- harmonic oscillator ``H = (p^2 + q^2) / 2``, exact solution a rotation;
- Kepler problem ``H = |p|^2 / 2 - 1 / |q|`` with eccentricity ``e`` and period ``2 pi``,
  exact solution through Kepler's equation;
- Tao's non-separable example ``H = (q^2 + 1)(p^2 + 1) / 2`` (Tao 2016, §4.1), whose
  position-dependent "kinetic energy" mimics the geodesic Hamiltonian.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

import numpy as np
from numba import njit

from ..integrators import Problem


@njit
def harmonic_rhs(y):
    out = np.empty(2)
    out[0] = y[1]
    out[1] = -y[0]
    return out


@njit
def kepler_rhs(y):
    q1, q2, p1, p2 = y[0], y[1], y[2], y[3]
    r3 = (q1 * q1 + q2 * q2) ** 1.5
    out = np.empty(4)
    out[0] = p1
    out[1] = p2
    out[2] = -q1 / r3
    out[3] = -q2 / r3
    return out


@njit
def tao_example_rhs(y):
    q, p = y[0], y[1]
    out = np.empty(2)
    out[0] = (q * q + 1.0) * p
    out[1] = -(p * p + 1.0) * q
    return out


@dataclass(frozen=True)
class ToyProblem:
    name: str
    rhs: Callable
    y0: np.ndarray
    energy: Callable[[np.ndarray], np.ndarray]
    exact: Callable[[float], np.ndarray] | None = None

    def problem(self, lam_end: float, events=()) -> Problem:
        return Problem(self.rhs, self.y0, (0.0, lam_end), events, hamiltonian=True)


def harmonic(q0: float = 1.0, p0: float = 0.0) -> ToyProblem:
    def exact(t):
        c, s = np.cos(t), np.sin(t)
        return np.array([q0 * c + p0 * s, p0 * c - q0 * s])

    return ToyProblem(
        "harmonic",
        harmonic_rhs,
        np.array([q0, p0]),
        lambda y: 0.5 * (np.asarray(y)[..., 0] ** 2 + np.asarray(y)[..., 1] ** 2),
        exact,
    )


def kepler(e: float = 0.5) -> ToyProblem:
    """Start at periapsis on the x axis; semi-major axis 1, so ``H = -1/2`` and period ``2 pi``."""

    def exact(t):
        M = np.mod(t, 2 * np.pi)
        E = M + e * np.sin(M)
        for _ in range(50):  # Newton on Kepler's equation E - e sin E = M
            step = (E - e * np.sin(E) - M) / (1.0 - e * np.cos(E))
            E -= step
            if abs(step) < 1e-16:
                break
        c, s = np.cos(E), np.sin(E)
        denom = 1.0 - e * c
        w = np.sqrt(1.0 - e * e)
        return np.array([c - e, w * s, -s / denom, w * c / denom])

    def energy(y):
        y = np.asarray(y)
        return 0.5 * (y[..., 2] ** 2 + y[..., 3] ** 2) - 1.0 / np.hypot(y[..., 0], y[..., 1])

    y0 = np.array([1.0 - e, 0.0, 0.0, np.sqrt((1.0 + e) / (1.0 - e))])
    return ToyProblem(f"kepler-e{e:g}", kepler_rhs, y0, energy, exact)


def tao_example(q0: float = 1.0, p0: float = 0.5) -> ToyProblem:
    def energy(y):
        y = np.asarray(y)
        return 0.5 * (y[..., 0] ** 2 + 1.0) * (y[..., 1] ** 2 + 1.0)

    return ToyProblem("tao-example", tao_example_rhs, np.array([q0, p0]), energy)
