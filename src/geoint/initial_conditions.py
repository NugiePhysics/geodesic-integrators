"""Initial data on the mass shell ``g^{mu nu}(x) p_mu p_nu = -eps``.

Every builder works in the canonical variables: fix the position ``x`` and three components of
the covariant momentum ``p``, then solve the mass shell, a quadratic equation, for the fourth.
:meth:`Formulation.from_xp <geoint.formulations.Formulation.from_xp>` turns the result into the
state of either formulation, so both start from the same physical point.

``eps = 1`` for timelike geodesics (affine parameter = proper time) and ``eps = 0`` for null
geodesics (normalized by the caller, conventionally ``E = 1`` so that ``Lz = b``).
"""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .metrics import Metric

#: Default tolerance below which a negative discriminant is treated as rounding at a turning point.
TURNING_POINT_RTOL = 64 * np.finfo(np.float64).eps


def solve_quadratic_shell(
    G: ArrayLike,
    v: ArrayLike,
    index: int,
    eps: float,
    sign: float,
    *,
    turning_point_rtol: float = TURNING_POINT_RTOL,
) -> NDArray[np.float64]:
    """Copy of ``v`` whose component ``index`` is chosen so that ``v^T G v = -eps``.

    With ``k = index`` the condition reads ``A v_k^2 + 2 B v_k + C = 0`` where ``A = G_kk``,
    ``B = sum_{j != k} G_kj v_j`` and ``C = sum_{i, j != k} G_ij v_i v_j + eps``. Both roots
    satisfy ``(G v)_k = A v_k + B = +-sqrt(D)`` with ``D = B^2 - A C``, so ``sign`` picks the
    root by the sign of ``(G v)_k``. For ``G = g^{-1}`` that is the contravariant velocity
    ``dx^k/dlambda``: ``sign = -1`` with ``k = 1`` is an ingoing particle, ``sign = +1`` with
    ``k = 0`` a future-directed one.

    The root is evaluated in whichever of the two equivalent forms ``(-B + sign sqrt(D)) / A``
    and ``C / (-B - sign sqrt(D))`` does not subtract nearly equal numbers. That matters as soon
    as the metric has off-diagonal terms (``B != 0``, e.g. Kerr).

    At a turning point ``D`` is itself a difference of nearly equal numbers and can come out
    slightly negative. If ``-D <= turning_point_rtol * scale``, where ``scale`` bounds the terms
    of ``D``, the root is taken as the double root ``D = 0``. A clearly negative ``D`` raises.

    Raises
    ------
    ValueError
        If no real root exists (the position is forbidden for the given components), or the
        requested root is at infinity (``A = 0``).
    """
    G = np.asarray(G, dtype=np.float64)
    v = np.array(v, dtype=np.float64)
    if sign not in (-1, 1):
        raise ValueError(f"sign must be +1 or -1, got {sign}")
    k = index
    others = np.arange(v.size) != k
    v[k] = 0.0
    A = G[k, k]
    B = G[k] @ v
    terms = G[np.ix_(others, others)] * np.outer(v[others], v[others])
    C = terms.sum() + eps
    D = B * B - A * C
    if D < 0.0:
        scale = B * B + abs(A) * (np.abs(terms).sum() + abs(eps))
        if -D > turning_point_rtol * scale:
            raise ValueError(
                f"no real solution for component {k}: discriminant {D:.3e} < 0 "
                "(the position is forbidden for the given constants)"
            )
        D = 0.0
    root = sign * np.sqrt(D)
    if sign * B > 0.0:
        v[k] = C / (-B - root)
    elif A != 0.0:
        v[k] = (-B + root) / A
    else:
        raise ValueError(f"the requested root for component {k} is at infinity (G[k, k] = 0)")
    return v


def complete_momentum(
    metric: Metric,
    x: ArrayLike,
    p: ArrayLike,
    index: int,
    eps: float,
    sign: float,
    *,
    turning_point_rtol: float = TURNING_POINT_RTOL,
) -> NDArray[np.float64]:
    """Covariant momentum on the mass shell at ``x``: ``p`` with component ``index`` solved for.

    ``sign`` is the sign of the resulting velocity component ``dx^index/dlambda``. See
    :func:`solve_quadratic_shell` for the root formula and the turning-point tolerance.
    """
    x = np.asarray(x, dtype=np.float64)
    return solve_quadratic_shell(
        metric.g_inv(x), p, index, eps, sign, turning_point_rtol=turning_point_rtol
    )


def momentum_from_constants(
    metric: Metric,
    x: ArrayLike,
    E: float,
    Lz: float,
    eps: float,
    *,
    p_theta: float = 0.0,
    radial_sign: float = -1.0,
) -> NDArray[np.float64]:
    """``p = (-E, p_r, p_theta, Lz)`` with ``p_r`` from the mass shell.

    ``radial_sign = -1`` starts the particle moving inwards, ``+1`` outwards.

    Near a turning point the discriminant ``E^2 - V_eff`` is tiny and carries an absolute
    rounding error of order ``1e-16``, so ``p_r`` comes out of order ``1e-8`` where it should be
    zero. That is harmless for a generic orbit (it corresponds to a rounding-level change of the
    constants) but not for a circular one. To start exactly at a turning point, use
    :func:`turning_point_momentum`, which sets ``p_r = 0`` and solves for ``E`` instead.
    """
    p = np.array([-E, 0.0, p_theta, Lz], dtype=np.float64)
    return complete_momentum(metric, x, p, 1, eps, radial_sign)


def turning_point_momentum(
    metric: Metric, x: ArrayLike, Lz: float, eps: float, *, p_theta: float = 0.0
) -> NDArray[np.float64]:
    """``p = (p_t, 0, p_theta, Lz)``: a radial turning point, with ``p_t = -E`` from the shell.

    The energy is the future-directed root (``dt/dlambda > 0``). Solving for ``E`` rather than
    ``p_r`` involves no cancellation, so a circular orbit started this way has ``p_r = 0``
    exactly and the constraint satisfied to rounding.
    """
    p = np.array([0.0, 0.0, p_theta, Lz], dtype=np.float64)
    return complete_momentum(metric, x, p, 0, eps, +1)
