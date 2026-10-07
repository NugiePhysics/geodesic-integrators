"""Butcher tableaus.

The Dormand-Prince coefficients are taken from SciPy, which is the oracle our own
implementations must reproduce step for step (ADR 0005). The Gauss-Legendre tableaus
are the closed-form coefficients, together with the two derived arrays the implicit solver
needs (GNI VIII.6): ``d = b^T A^{-1}``, which gives the update directly from the stage
increments, and the matrix that extrapolates the previous step's collocation polynomial into
a starting guess.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.polynomial import Polynomial
from scipy.integrate._ivp import dop853_coefficients as _dop853
from scipy.integrate._ivp.rk import RK45 as _RK45


@dataclass(frozen=True)
class Tableau:
    A: np.ndarray
    b: np.ndarray
    c: np.ndarray
    order: int

    @property
    def stages(self) -> int:
        return self.b.size


RK4 = Tableau(
    A=np.array([[0, 0, 0, 0], [0.5, 0, 0, 0], [0, 0.5, 0, 0], [0, 0, 1.0, 0]]),
    b=np.array([1 / 6, 1 / 3, 1 / 3, 1 / 6]),
    c=np.array([0, 0.5, 0.5, 1.0]),
    order=4,
)


@dataclass(frozen=True)
class EmbeddedTableau:
    """An explicit pair in SciPy's layout (``E`` acts on the stages plus the FSAL stage)."""

    name: str
    A: np.ndarray
    B: np.ndarray
    C: np.ndarray
    E: np.ndarray  # DP5: error weights. DOP853: the 5th-order estimate E5.
    E3: np.ndarray  # DOP853 only: the 3rd-order estimate. Zeros for DP5.
    order: int
    error_estimator_order: int

    @property
    def stages(self) -> int:
        return self.B.size


DP5 = EmbeddedTableau(
    name="DP5",
    A=np.ascontiguousarray(_RK45.A, dtype=float),
    B=np.ascontiguousarray(_RK45.B, dtype=float),
    C=np.ascontiguousarray(_RK45.C, dtype=float),
    E=np.ascontiguousarray(_RK45.E, dtype=float),
    E3=np.zeros_like(_RK45.E, dtype=float),
    order=5,
    error_estimator_order=4,
)

_S = _dop853.N_STAGES
DOP853 = EmbeddedTableau(
    name="DOP853",
    A=np.ascontiguousarray(_dop853.A[:_S, :_S], dtype=float),
    B=np.ascontiguousarray(_dop853.B, dtype=float),
    C=np.ascontiguousarray(_dop853.C[:_S], dtype=float),
    E=np.ascontiguousarray(_dop853.E5, dtype=float),
    E3=np.ascontiguousarray(_dop853.E3, dtype=float),
    order=8,
    error_estimator_order=7,
)


@dataclass(frozen=True)
class GaussTableau(Tableau):
    d: np.ndarray  # b^T A^{-1}
    A_inv: np.ndarray
    extrapolation: np.ndarray  # X[i, j] = beta_j(1 + c_i) - b_j


def _gauss_coefficients(s: int):
    """Closed-form ``(A, b, c)`` of the ``s``-stage Gauss method (HNW Table II.7.3, 7.4)."""
    if s == 1:
        return np.array([[0.5]]), np.array([1.0]), np.array([0.5])
    if s == 2:
        r = np.sqrt(3.0)
        A = np.array([[1 / 4, 1 / 4 - r / 6], [1 / 4 + r / 6, 1 / 4]])
        return A, np.array([0.5, 0.5]), np.array([0.5 - r / 6, 0.5 + r / 6])
    if s == 3:
        r = np.sqrt(15.0)
        A = np.array(
            [
                [5 / 36, 2 / 9 - r / 15, 5 / 36 - r / 30],
                [5 / 36 + r / 24, 2 / 9, 5 / 36 - r / 24],
                [5 / 36 + r / 30, 2 / 9 + r / 15, 5 / 36],
            ]
        )
        return A, np.array([5 / 18, 4 / 9, 5 / 18]), np.array([0.5 - r / 10, 0.5, 0.5 + r / 10])
    raise ValueError("Gauss-Legendre is implemented for s = 1, 2, 3")


def gauss_legendre(s: int) -> GaussTableau:
    """The ``s``-stage Gauss-Legendre collocation method (order ``2s``)."""
    A, b, c = _gauss_coefficients(s)
    # beta_j(theta) = int_0^theta l_j, with l_j the Lagrange basis on the nodes. The
    # extrapolation matrix only shapes the starting guess, so its rounding does not matter.
    beta = []
    for j in range(s):
        others = np.delete(c, j)
        poly = Polynomial([1.0])
        for root in others:
            poly = poly * Polynomial([-root, 1.0]) / (c[j] - root)
        beta.append(poly.integ())
    X = np.array([[beta[j](1.0 + c[i]) - b[j] for j in range(s)] for i in range(s)])
    A_inv = np.linalg.inv(A)
    return GaussTableau(A=A, b=b, c=c, order=2 * s, d=b @ A_inv, A_inv=A_inv, extrapolation=X)


GL1, GL2, GL3 = (gauss_legendre(s) for s in (1, 2, 3))
