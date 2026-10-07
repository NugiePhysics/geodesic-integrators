"""Closed-form and quadrature references for timelike Schwarzschild orbits (``M = 1``).

Bound orbits are labelled by the semi-latus rectum ``p`` and eccentricity ``e`` (Cutler,
Kennefick & Poisson 1994): ``r_p = p/(1+e)``, ``r_a = p/(1-e)``. With ``u = 1/r`` the orbit
equation is ``(du/dphi)^2 = 2 (u - u1)(u - u2)(u - u3)`` with ``u1 = (1-e)/p``,
``u2 = (1+e)/p`` and ``u3 = 1/2 - 2/p``, which gives

- the angle per radial period ``Phi = 4 K(m) / sqrt(2 (u3 - u1))``, ``m = (u2-u1)/(u3-u1)``
  (SciPy's parameter convention ``m = k^2``, roadmap pitfall 17);
- the orbit ``u(phi) = u1 + (u2 - u1) cd^2(sqrt((u3 - u1)/2) phi | m)`` with ``phi`` measured
  from periapsis.

Periods in proper and coordinate time, and ``Phi`` a second time, come from an independent
mpmath quadrature in the angle ``chi`` of ``r = (r_a + r_p)/2 - (r_a - r_p)/2 cos chi``, which
removes the square-root singularities at both turning points.
"""

from __future__ import annotations

import functools
import math
from dataclasses import dataclass

import mpmath
import numpy as np
from scipy.special import ellipj, ellipk

ISCO = 6.0
PHOTON_SPHERE = 3.0
B_CRIT = 3.0 * math.sqrt(3.0)
#: Lyapunov exponent of the photon sphere in coordinate time; equal to its orbital frequency.
LYAPUNOV_PHOTON_SPHERE = 1.0 / B_CRIT


def circular_orbit_constants(r: float) -> tuple[float, float]:
    """``(E, L)`` of the timelike circular orbit of radius ``r > 3``."""
    if not r > 3.0:
        raise ValueError("timelike circular orbits need r > 3M")
    s = math.sqrt(1.0 - 3.0 / r)
    return (1.0 - 2.0 / r) / s, math.sqrt(r) / s


def circular_orbit_frequency(r: float) -> float:
    """Coordinate angular velocity ``dphi/dt = r^(-3/2)``."""
    return r**-1.5


def circular_orbit_period(r: float) -> tuple[float, float]:
    """``(T_tau, T_t)``: one revolution in proper and in coordinate time."""
    _, L = circular_orbit_constants(r)
    return 2.0 * math.pi * r * r / L, 2.0 * math.pi * r**1.5


def pe_to_EL(p: float, e: float) -> tuple[float, float]:
    """Energy and angular momentum of the orbit ``(p, e)`` (Cutler et al. 1994)."""
    if not p > 6.0 + 2.0 * e:
        raise ValueError(f"(p, e) = ({p}, {e}) is not a stable bound orbit (p <= 6 + 2e)")
    denom = p - 3.0 - e * e
    return math.sqrt(((p - 2.0) ** 2 - 4.0 * e * e) / (p * denom)), p / math.sqrt(denom)


def EL_to_pe(E: float, L: float) -> tuple[float, float]:
    """Inverse of :func:`pe_to_EL`, from the two largest roots of the radial cubic."""
    roots = np.sort(np.roots([E * E - 1.0, 2.0, -L * L, 2.0 * L * L]).real)
    r_p, r_a = roots[1], roots[2]
    return 2.0 * r_a * r_p / (r_a + r_p), (r_a - r_p) / (r_a + r_p)


def separatrix(e: float) -> float:
    """Smallest ``p`` of a bound orbit with eccentricity ``e``."""
    return 6.0 + 2.0 * e


@dataclass(frozen=True)
class EccentricOrbit:
    """Exact description of the bound timelike orbit ``(p, e)``, started at periapsis."""

    p: float
    e: float

    @property
    def r_p(self) -> float:
        return self.p / (1.0 + self.e)

    @property
    def r_a(self) -> float:
        return self.p / (1.0 - self.e)

    @property
    def r3(self) -> float:
        """Third root of the radial cubic, inside the periapsis."""
        return 2.0 * self.p / (self.p - 4.0)

    @property
    def EL(self) -> tuple[float, float]:
        return pe_to_EL(self.p, self.e)

    def _u(self):
        u1, u2 = (1.0 - self.e) / self.p, (1.0 + self.e) / self.p
        return u1, u2, 0.5 - 2.0 / self.p

    @property
    def m(self) -> float:
        u1, u2, u3 = self._u()
        return (u2 - u1) / (u3 - u1)

    @property
    def Phi(self) -> float:
        """Angle swept per radial period (closed form)."""
        u1, _, u3 = self._u()
        return 4.0 * ellipk(self.m) / math.sqrt(2.0 * (u3 - u1))

    @property
    def precession(self) -> float:
        """Periapsis advance per radial period, ``Phi - 2 pi``."""
        return self.Phi - 2.0 * math.pi

    def r_of_phi(self, phi):
        """Radius at angle ``phi`` from periapsis (closed form, Jacobi elliptic functions)."""
        u1, u2, u3 = self._u()
        sn, cn, dn, _ = ellipj(math.sqrt(0.5 * (u3 - u1)) * np.asarray(phi, dtype=float), self.m)
        return 1.0 / (u1 + (u2 - u1) * (cn / dn) ** 2)

    @functools.cached_property
    def periods(self) -> dict[str, float]:
        """``Phi``, ``T_tau`` and ``T_t`` per radial period, by 30-digit quadrature."""
        return _periods(self.p, self.e)


@functools.cache
def _periods(p: float, e: float) -> dict[str, float]:
    with mpmath.workdps(30):
        p_, e_ = mpmath.mpf(p), mpmath.mpf(e)
        r_p, r_a = p_ / (1 + e_), p_ / (1 - e_)
        r3 = 2 * p_ / (p_ - 4)
        one_minus_E2 = (1 - e_**2) * (p_ - 4) / (p_ * (p_ - 3 - e_**2))
        E = mpmath.sqrt(1 - one_minus_E2)
        L = p_ / mpmath.sqrt(p_ - 3 - e_**2)

        def r(chi):
            return (r_a + r_p) / 2 - (r_a - r_p) / 2 * mpmath.cos(chi)

        def dtau(chi):  # d tau / d chi
            rr = r(chi)
            return rr**1.5 / mpmath.sqrt(one_minus_E2 * (rr - r3))

        T_tau = 2 * mpmath.quad(dtau, [0, mpmath.pi])
        Phi = 2 * mpmath.quad(lambda c: L / r(c) ** 2 * dtau(c), [0, mpmath.pi])
        T_t = 2 * mpmath.quad(lambda c: E / (1 - 2 / r(c)) * dtau(c), [0, mpmath.pi])
        return {"Phi": float(Phi), "T_tau": float(T_tau), "T_t": float(T_t)}
