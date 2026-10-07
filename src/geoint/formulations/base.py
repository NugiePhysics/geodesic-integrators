"""What the two formulations of the geodesic problem have in common.

A formulation turns a :class:`~geoint.metrics.Metric` into an autonomous ODE ``dy/dlambda =
rhs(y)`` on an 8-component state ``y``, ordered ``(t, r, theta, phi)`` followed by the velocity
(formulation a) or the momentum (formulation b) in the same order (ADR 0001). Integrators only
ever see ``rhs``; everything else here is for setting up runs and measuring them.
"""

from __future__ import annotations

from collections.abc import Callable
from typing import NamedTuple

import numpy as np
from numba import njit
from numpy.typing import ArrayLike, NDArray

from ..metrics import Metric

Kernel = Callable[..., NDArray[np.float64]]


class FormulationKernels(NamedTuple):
    """Compiled functions of a formulation; ``rhs(y)`` and ``point_invariants(y)`` take ``(8,)``."""

    rhs: Kernel
    point_invariants: Kernel


@njit
def _apply_rows(fn, Y, n_out):
    """``fn`` applied to every row of the 2-D array ``Y``."""
    out = np.empty((Y.shape[0], n_out))
    for i in range(Y.shape[0]):
        out[i] = fn(Y[i])
    return out


class Formulation:
    """Geodesics of ``metric`` as a first-order system in an 8-component state.

    Attributes
    ----------
    metric : Metric
    rhs : kernel ``y -> dy/dlambda``
        The vector field; the only thing an integrator sees.
    invariant_names : tuple of str
        ``"H"`` (the mass shell, ``-eps/2`` on a geodesic) followed by the metric's invariants.
    """

    name: str = "formulation"
    dim: int = 8

    def __init__(self, metric: Metric, kernels: FormulationKernels):
        self.metric = metric
        self.kernels = kernels
        self.rhs = kernels.rhs
        self.invariant_names = ("H", *metric.invariant_names)

    def from_xp(self, x: ArrayLike, p: ArrayLike) -> NDArray[np.float64]:
        """State of this formulation at position ``x`` with covariant momentum ``p``."""
        raise NotImplementedError

    def to_xp(self, y: ArrayLike) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        """Position and covariant momentum ``(x, p)`` of a single state ``y``."""
        raise NotImplementedError

    def invariants(self, y: ArrayLike) -> dict[str, NDArray[np.float64]]:
        """Invariants of a state ``(8,)`` or a trajectory ``(N, 8)``, keyed by name."""
        Y = np.asarray(y, dtype=np.float64)
        rows = np.ascontiguousarray(Y.reshape(-1, self.dim))
        values = _apply_rows(self.kernels.point_invariants, rows, len(self.invariant_names))
        if Y.ndim == 1:
            values = values[0]
        return {name: values[..., i] for i, name in enumerate(self.invariant_names)}

    def constraint(self, y: ArrayLike, eps: float) -> NDArray[np.float64]:
        """Signed violation ``C = H + eps/2`` of the mass shell (zero on a geodesic)."""
        return self.invariants(y)["H"] + 0.5 * eps

    def relative_constraint(self, y: ArrayLike, eps: float) -> NDArray[np.float64]:
        """``delta H`` of theory §5.1: ``|H + 1/2| / (1/2)`` if timelike, ``|H| / E^2`` if null.

        The null normalization is invariant under the rescaling ``p -> alpha p`` that the free
        normalization of a null affine parameter allows.
        """
        inv = self.invariants(y)
        C = np.abs(inv["H"] + 0.5 * eps)
        return C / (0.5 * eps) if eps > 0 else C / inv["E"] ** 2

    def __repr__(self) -> str:
        return f"{type(self).__name__}({self.metric!r})"
