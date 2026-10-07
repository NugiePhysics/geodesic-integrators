"""Light deflection by a Schwarzschild black hole (``M = 1``), exact and approximate.

A photon with impact parameter ``b > b_c = 3 sqrt(3)`` reaches the periapsis ``r0 = 1/u2``,
where ``u1 < 0 < u2 < u3`` are the roots of ``2u^3 - u^2 + 1/b^2``. The angle swept between two
passages through ``r = r_far`` is

    Delta phi(r_far) = 2 int_{u_far}^{u2} du / sqrt(1/b^2 - u^2 + 2u^3)
                     = 4 F(psi | m) / sqrt(2 (u3 - u1)),
    m = (u2 - u1)/(u3 - u1),   sin^2 psi = (u3 - u1)(u2 - u_far) / ((u2 - u1)(u3 - u_far)).

Integrations stop at a finite ``r_far``, so they are compared with ``Delta phi(r_far)``, not with
the asymptotic angle (pitfall 18, docs/pitfalls.md). Three independent evaluations are provided: the
elliptic form above, Darwin's (1959) form in terms of ``r0``, and a direct quadrature.
"""

from __future__ import annotations

import math

import mpmath

from .schwarzschild_orbits import B_CRIT

DPS = 30
#: Constant term of Bozza's strong-deflection limit, ``log(216 (7 - 4 sqrt 3)) - pi``.
BOZZA_B = math.log(216.0 * (7.0 - 4.0 * math.sqrt(3.0))) - math.pi


def _check(b):
    if not b > B_CRIT:
        raise ValueError(f"b = {b} is captured (b_c = {B_CRIT})")


def photon_roots(b: float) -> tuple[mpmath.mpf, mpmath.mpf, mpmath.mpf]:
    """``(u1, u2, u3)``, roots of ``2u^3 - u^2 + 1/b^2``, at the current mpmath precision."""
    _check(b)
    # Trigonometric solution of the depressed cubic; all three roots are real for b > b_c.
    b = mpmath.mpf(b)
    theta = mpmath.acos(1 - 54 / b**2)
    roots = sorted((1 + 2 * mpmath.cos((theta + 2 * mpmath.pi * k) / 3)) / 6 for k in range(3))
    return roots[0], roots[1], roots[2]


def periapsis(b: float) -> float:
    """Closest approach ``r0`` of a photon with impact parameter ``b``."""
    with mpmath.workdps(DPS):
        return float(1 / photon_roots(b)[1])


def _delta_phi(b: float, r_far: float) -> mpmath.mpf:
    u1, u2, u3 = photon_roots(b)
    u_far = mpmath.mpf(0) if math.isinf(r_far) else 1 / mpmath.mpf(r_far)
    if not u_far < u2:
        raise ValueError("r_far must lie outside the periapsis")
    m = (u2 - u1) / (u3 - u1)
    sin2 = (u3 - u1) * (u2 - u_far) / ((u2 - u1) * (u3 - u_far))
    psi = mpmath.asin(mpmath.sqrt(sin2))
    return 4 * mpmath.ellipf(psi, m) / mpmath.sqrt(2 * (u3 - u1))


def delta_phi(b: float, r_far: float = math.inf) -> float:
    """Exact ``Delta phi(r_far)`` from the elliptic integral of the first kind."""
    with mpmath.workdps(DPS):
        return float(_delta_phi(b, r_far))


def deflection_angle(b: float, r_far: float = math.inf) -> float:
    """Deflection ``alpha = Delta phi - pi`` (asymptotic by default), without cancellation."""
    with mpmath.workdps(DPS):
        return float(_delta_phi(b, r_far) - mpmath.pi)


def deflection_darwin(b: float) -> float:
    """Darwin's (1959) closed form ``alpha = -pi + 4 sqrt(r0/Q) [K(k) - F(zeta, k)]``.

    In the notation of Iyer & Petters (2007), eq. (6), with the parameter ``m = k^2``:

    ``Q^2 = (r0 - 2)(r0 + 6)``, ``k^2 = (Q - r0 + 6)/(2Q)``,
    ``sin^2 zeta = (Q - r0 + 2)/(Q - r0 + 6)``.
    """
    with mpmath.workdps(DPS):
        r0 = 1 / photon_roots(b)[1]
        Q = mpmath.sqrt((r0 - 2) * (r0 + 6))
        k2 = (Q - r0 + 6) / (2 * Q)
        zeta = mpmath.asin(mpmath.sqrt((Q - r0 + 2) / (Q - r0 + 6)))
        integral = mpmath.ellipk(k2) - mpmath.ellipf(zeta, k2)
        return float(-mpmath.pi + 4 * mpmath.sqrt(r0 / Q) * integral)


def delta_phi_quadrature(b: float, r_far: float = math.inf) -> float:
    """``Delta phi(r_far)`` by direct quadrature, independent of ``u1`` and ``u3``.

    With ``u = u2 - s^2`` and ``P(u) = 1/b^2 - u^2 + 2u^3``, the integrand ``2 du / sqrt(P)``
    becomes ``4 ds / sqrt(P(u2 - s^2) / s^2)``; the quotient is the exact Taylor polynomial
    ``-P'(u2) + P''(u2) s^2 / 2 - 2 s^4``, free of the cancellation in ``P`` near ``u2``.
    """
    with mpmath.workdps(DPS + 10):
        u2 = photon_roots(b)[1]
        u_far = mpmath.mpf(0) if math.isinf(r_far) else 1 / mpmath.mpf(r_far)
        d1 = -2 * u2 + 6 * u2**2  # P'(u2)
        d2 = -2 + 12 * u2  # P''(u2)

        def integrand(s):
            s2 = s * s
            return 4 / mpmath.sqrt(-d1 + d2 / 2 * s2 - 2 * s2 * s2)

        return float(mpmath.quad(integrand, [0, mpmath.sqrt(u2 - u_far)]))


def turning_point_delta_phi(r0: float, r_exit: float) -> float:
    """Angle swept by a photon from a turning point ``r0 > 3`` outward to ``r_exit > r0``.

    The photon has ``b^2 = r0^3 / (r0 - 2)``. Its turning point ``u2 = 1/r0`` is known exactly,
    so the other two roots follow without solving the cubic:
    ``u1, u3 = [(1/2 - u2) -+ sqrt((1/2 - u2)(3 u2 + 1/2))] / 2``. This stays accurate as
    ``r0 -> 3``, where ``b - b_c = O((r0 - 3)^2)`` is below float64 resolution (test case TC2).
    """
    if not 3.0 < r0 < r_exit:
        raise ValueError("need 3 < r0 < r_exit")
    with mpmath.workdps(2 * DPS):
        u2 = 1 / mpmath.mpf(r0)
        half = mpmath.mpf(1) / 2
        root = mpmath.sqrt((half - u2) * (3 * u2 + half))
        u1, u3 = ((half - u2) - root) / 2, ((half - u2) + root) / 2
        u_exit = 1 / mpmath.mpf(r_exit)
        m = (u2 - u1) / (u3 - u1)
        sin2 = (u3 - u1) * (u2 - u_exit) / ((u2 - u1) * (u3 - u_exit))
        psi = mpmath.asin(mpmath.sqrt(sin2))
        return float(2 * mpmath.ellipf(psi, m) / mpmath.sqrt(2 * (u3 - u1)))


def weak_field_deflection(b: float, order: int = 4) -> float:
    """Series ``4/b + (15 pi/4)/b^2 + (128/3)/b^3 + (3465 pi/64)/b^4``, first ``order`` terms."""
    if not 1 <= order <= 4:
        raise ValueError("order must be between 1 and 4")
    coefficients = [4.0, 15.0 * math.pi / 4.0, 128.0 / 3.0, 3465.0 * math.pi / 64.0]
    return sum(c / b ** (n + 1) for n, c in enumerate(coefficients[:order]))


def bozza_deflection(b: float) -> float:
    """Strong-deflection limit (Bozza 2002): ``-log(b/b_c - 1) + log(216 (7 - 4 sqrt 3)) - pi``.

    ``b/b_c - 1`` is formed in extended precision: in float64 it would lose
    ``log10(b_c / (b - b_c))`` digits right where the limit applies.
    """
    _check(b)
    with mpmath.workdps(DPS):
        x = mpmath.mpf(b) / (3 * mpmath.sqrt(3)) - 1
        return float(-mpmath.log(x) + mpmath.log(216 * (7 - 4 * mpmath.sqrt(3))) - mpmath.pi)


def windings(b: float, r_far: float = math.inf) -> float:
    """Number of revolutions ``Delta phi / (2 pi)`` around the hole."""
    return delta_phi(b, r_far) / (2.0 * math.pi)
