"""The test-case matrix TC1-TC6 of the roadmap (§1.3) for the Schwarzschild metric, ``M = 1``.

Errors are always measured on physical quantities (``r``, ``phi``, ``t`` at an event or at a
given proper time), never on the raw state, because the states of (a) and (b) differ.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from .. import analytic
from ..initial_conditions import momentum_from_constants, turning_point_momentum
from ..integrators import Event
from .base import (
    PHI_INDEX,
    R_INDEX,
    T_INDEX,
    THETA_INDEX,
    GeodesicCase,
    capture_event,
    periapsis_event,
)

EQUATOR = math.pi / 2


def _failed(solution, expected_status: str = "completed") -> dict | None:
    if solution.status != expected_status:
        return {"error": math.inf, "status": solution.status}
    return None


class _Deflection(GeodesicCase):
    """A photon from ``r_far`` past the hole and back out to ``r_far`` (``E = 1``, ``L = b``)."""

    eps = 0
    b: float
    r_far: float

    def initial_xp(self, metric):
        x = np.array([0.0, self.r_far, EQUATOR, 0.0])
        return x, momentum_from_constants(metric, x, 1.0, self.b, 0, radial_sign=-1.0)

    @property
    def lam_end(self) -> float:
        return 3.0 * self.r_far + 1000.0

    def events(self, metric):
        return (capture_event(metric), Event("escaped", R_INDEX, self.r_far, +1, stop_after=1))

    def reference(self) -> dict:
        delta_phi = analytic.delta_phi(self.b, self.r_far)
        return {
            "delta_phi": delta_phi,
            "windings": delta_phi / (2 * math.pi),
            "alpha_inf": analytic.deflection_angle(self.b),
        }

    def errors(self, solution, formulation) -> dict:
        failed = _failed(solution, "escaped")
        if failed:
            return failed
        ref = self.reference()
        delta_phi = solution.y_end[PHI_INDEX] - solution.y[0, PHI_INDEX]
        return {
            "error": abs(delta_phi - ref["delta_phi"]),
            "delta_phi": delta_phi,
            "windings": delta_phi / (2 * math.pi),
            **self.constraint_errors(solution, formulation),
        }


@dataclass(frozen=True)
class Deflection(_Deflection):
    """TC1: light deflection with impact parameter ``b``."""

    b: float
    r_far: float = 1000.0
    family = "TC1"


@dataclass(frozen=True)
class NearCritical(_Deflection):
    """TC1b: ``b = b_c (1 + 10^-k)``, winding ``~ -log(b/b_c - 1) / (2 pi)`` times."""

    k: float
    r_far: float = 1000.0
    family = "TC1b"

    @property
    def b(self) -> float:
        return analytic.B_CRIT * (1.0 + 10.0**-self.k)

    def reference(self) -> dict:
        ref = super().reference()
        ref["alpha_bozza"] = analytic.bozza_deflection(self.b)
        return ref


def in_plane_angle(solution, formulation) -> np.ndarray:
    """Angle swept in the orbital plane at every saved point, about the initial ``L`` vector.

    For a prograde orbit (``L_z > 0``, inclination below 90 degrees) the in-plane angle and
    ``phi`` advance together and never differ by more than ``pi/2``. Taking the number of whole
    turns from ``phi``, which is continuous in the state, makes the result independent of how
    far apart the saved points are (``np.unwrap`` would lose turns on large adaptive steps).
    """
    L0 = angular_momentum_vector(*formulation.to_xp(solution.y[0]))
    n = unit_position(solution.y)
    axis = L0 / np.linalg.norm(L0)
    raw = np.arctan2(np.cross(n[0], n) @ axis, n @ n[0])
    phi = solution.y[:, PHI_INDEX] - solution.y[0, PHI_INDEX]
    return phi + np.angle(np.exp(1j * (raw - phi)))


def angular_momentum_vector(x, p) -> np.ndarray:
    """Cartesian ``(L_x, L_y, L_z)`` from the rotation Killing vectors (theory eq. 18)."""
    _, _, theta, phi = x
    cot = math.cos(theta) / math.sin(theta)
    return np.array(
        [
            -math.sin(phi) * p[2] - cot * math.cos(phi) * p[3],
            math.cos(phi) * p[2] - cot * math.sin(phi) * p[3],
            p[3],
        ]
    )


def unit_position(y) -> np.ndarray:
    theta, phi = y[..., THETA_INDEX], y[..., PHI_INDEX]
    return np.stack(
        [np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta)], axis=-1
    )


@dataclass(frozen=True)
class PhotonSphere(GeodesicCase):
    """TC2: photon at a turning point ``r0 = 3 + 10^-k`` with ``b = b_c`` (``k = inf``: r0 = 3).

    The unstable circular orbit amplifies an offset like ``cosh(lambda_L t)``: the linearized
    prediction for leaving ``|r - 3| < 0.1`` is ``arccosh(0.1 / delta0) / (2 pi)`` orbits. The
    exact number follows from the elliptic integral from the turning point to ``r = 3.1``.

    In the equatorial plane the circular orbit ``r = 3`` is an exact fixed point of every
    Runge-Kutta map: the ``(r, p_r)`` components of the vector field vanish there and ``t``,
    ``phi`` are linear in ``lam``, so no truncation error can perturb it. Tilting the orbital
    plane by ``inclination`` degrees makes ``theta`` oscillate; the truncation error in ``L^2``
    then perturbs ``r``, which is what makes the orbit's lifetime depend on the integrator.
    """

    k: float = math.inf
    inclination: float = 0.0
    family = "TC2"
    eps = 0
    exit_offset = 0.1

    @property
    def delta0(self) -> float:
        """The offset actually represented: ``10^-k`` rounded to the grid of doubles near 3.

        Below half an ulp of 3 (``k >= 16``) ``r0`` rounds to 3 and the offset is 0.
        """
        return 0.0 if math.isinf(self.k) else (3.0 + 10.0**-self.k) - 3.0

    def initial_xp(self, metric):
        i = math.radians(self.inclination)
        b = analytic.B_CRIT
        x = np.array([0.0, 3.0 + self.delta0, EQUATOR, 0.0])
        return x, turning_point_momentum(metric, x, b * math.cos(i), 0, p_theta=-b * math.sin(i))

    @property
    def lam_end(self) -> float:
        # About 55 orbits of the photon sphere; every perturbed orbit leaves well before.
        return 600.0

    def events(self, metric):
        d = self.exit_offset
        return (
            Event("left_outward", R_INDEX, 3.0 + d, +1, stop_after=1),
            Event("left_inward", R_INDEX, 3.0 - d, -1, stop_after=1),
        )

    def reference(self) -> dict:
        exact = linear = math.nan
        if self.delta0 > 0:
            r_exit = 3.0 + self.exit_offset
            exact = analytic.turning_point_delta_phi(3.0 + self.delta0, r_exit) / (2 * math.pi)
            linear = math.acosh(self.exit_offset / self.delta0) / (2 * math.pi)
        return {
            "orbits": exact,
            "orbits_linear": linear,
            "lyapunov": analytic.LYAPUNOV_PHOTON_SPHERE,
        }

    def errors(self, solution, formulation) -> dict:
        if solution.status.startswith("failed"):
            return {"error": math.inf, "status": solution.status}
        y = solution.y
        orbits = in_plane_angle(solution, formulation)[-1] / (2 * math.pi)
        offset = np.abs(y[:, R_INDEX] - 3.0)
        window = (offset > max(10 * self.delta0, 1e-13)) & (offset < 0.03)
        rate = math.nan
        if window.sum() >= 3:
            rate = float(np.polyfit(y[window, T_INDEX], np.log(offset[window]), 1)[0])
        predicted = self.reference()["orbits"]
        return {
            "error": abs(orbits - predicted) if self.delta0 > 0 else math.nan,
            "orbits": orbits,
            "t_exit": float(solution.y_end[T_INDEX]),
            "exit": solution.status,
            "growth_rate": rate,
            **self.constraint_errors(solution, formulation),
        }


@dataclass(frozen=True)
class Circular(GeodesicCase):
    """TC3: timelike circular orbit of radius ``r_c`` for ``n_orbits`` revolutions."""

    r_c: float
    n_orbits: float = 5.0
    family = "TC3"

    @property
    def constants(self):
        return analytic.circular_orbit_constants(self.r_c)

    def initial_xp(self, metric):
        x = np.array([0.0, self.r_c, EQUATOR, 0.0])
        return x, turning_point_momentum(metric, x, self.constants[1], 1)

    @property
    def lam_end(self) -> float:
        return self.n_orbits * analytic.circular_orbit_period(self.r_c)[0]

    def reference(self) -> dict:
        E, L = self.constants
        return {
            "E": E,
            "L": L,
            "omega": analytic.circular_orbit_frequency(self.r_c),
            "phi_end": L / self.r_c**2 * self.lam_end,
        }

    def errors(self, solution, formulation) -> dict:
        failed = _failed(solution)
        if failed:
            return failed
        ref = self.reference()
        phi, t = solution.y_end[PHI_INDEX], solution.y_end[T_INDEX]
        return {
            "error": abs(phi - ref["phi_end"]),
            "dr_max": float(np.max(np.abs(solution.y[:, R_INDEX] - self.r_c))),
            "omega": phi / t,
            "omega_rel_error": abs(phi / t - ref["omega"]) / ref["omega"],
            **self.constraint_errors(solution, formulation),
        }


@dataclass(frozen=True)
class MarginalCircular(GeodesicCase):
    """ISCO (r = 6) circular orbit in an inclined plane, followed until ``|r - 6| > 1``.

    The ISCO is marginally stable (an inflection point of the effective potential), so a
    perturbation grows algebraically rather than exponentially until the orbit plunges. In the
    equatorial plane it is an exact fixed point of every Runge-Kutta map (like the photon
    sphere, TC2); the inclination lets truncation errors in ``L^2`` perturb ``r``. Used for
    the time-to-plunge experiment (F11).
    """

    inclination: float = 45.0
    n_orbits: float = 2000.0
    family = "TC3m"
    exit_offset = 1.0

    def initial_xp(self, metric):
        _, L = analytic.circular_orbit_constants(analytic.ISCO)
        i = math.radians(self.inclination)
        x = np.array([0.0, analytic.ISCO, EQUATOR, 0.0])
        p = turning_point_momentum(metric, x, L * math.cos(i), 1, p_theta=-L * math.sin(i))
        return x, p

    @property
    def period(self) -> float:
        return analytic.circular_orbit_period(analytic.ISCO)[0]

    @property
    def lam_end(self) -> float:
        return self.n_orbits * self.period

    def events(self, metric):
        d = self.exit_offset
        return (
            capture_event(metric),
            Event("plunged", R_INDEX, analytic.ISCO - d, -1, stop_after=1),
            Event("outward", R_INDEX, analytic.ISCO + d, +1, stop_after=1),
        )

    def errors(self, solution, formulation) -> dict:
        if solution.status.startswith("failed"):
            return {"error": math.inf, "status": solution.status}
        return {
            "error": math.nan,
            "exit": solution.status,
            "lam_exit": solution.lam_end,
            "orbits": solution.lam_end / self.period,
            **self.constraint_errors(solution, formulation),
        }


@dataclass(frozen=True)
class Eccentric(GeodesicCase):
    """TC4: bound orbit ``(p, e)`` from periapsis, for ``n_orbits`` radial periods."""

    p: float
    e: float
    n_orbits: float = 1.0
    family = "TC4"

    @property
    def orbit(self) -> analytic.EccentricOrbit:
        return analytic.EccentricOrbit(self.p, self.e)

    def initial_xp(self, metric):
        o = self.orbit
        x = np.array([0.0, o.r_p, EQUATOR, 0.0])
        return x, turning_point_momentum(metric, x, o.EL[1], 1)

    @property
    def lam_end(self) -> float:
        return self.n_orbits * self.orbit.periods["T_tau"]

    def events(self, metric):
        return (capture_event(metric), periapsis_event())

    def reference(self) -> dict:
        o = self.orbit
        periods = o.periods
        return {
            "Phi": o.Phi,
            "precession": o.precession,
            "precession_weak_field": 6 * math.pi / self.p,
            "T_tau": periods["T_tau"],
            "T_t": periods["T_t"],
        }

    def errors(self, solution, formulation) -> dict:
        failed = _failed(solution)
        if failed:
            return failed
        ref = self.reference()
        y_end = solution.y_end
        dphi = y_end[PHI_INDEX] - self.n_orbits * ref["Phi"]
        shape = np.abs(solution.y[:, R_INDEX] - self.orbit.r_of_phi(solution.y[:, PHI_INDEX]))
        out = {
            "error": abs(dphi),
            "dphi": dphi,
            "dr": y_end[R_INDEX] - self.orbit.r_p,
            "dt": y_end[T_INDEX] - self.n_orbits * ref["T_t"],
            "shape_error": float(np.max(shape)),
            **self.constraint_errors(solution, formulation),
        }
        peri = solution.events.get("periapsis")
        if peri is not None and len(peri.lam) > 0:
            out["precession"] = peri.y[0, PHI_INDEX] - 2 * math.pi
        return out


@dataclass(frozen=True)
class Inclined(GeodesicCase):
    """TC6: the orbit ``(p, e)`` in a plane inclined by ``inclination`` degrees."""

    p: float
    e: float
    inclination: float
    n_orbits: float = 1.0
    family = "TC6"

    @property
    def orbit(self) -> analytic.EccentricOrbit:
        return analytic.EccentricOrbit(self.p, self.e)

    def initial_xp(self, metric):
        o = self.orbit
        L = o.EL[1]
        i = math.radians(self.inclination)
        x = np.array([0.0, o.r_p, EQUATOR, 0.0])
        p = turning_point_momentum(metric, x, L * math.cos(i), 1, p_theta=-L * math.sin(i))
        return x, p

    @property
    def lam_end(self) -> float:
        return self.n_orbits * self.orbit.periods["T_tau"]

    def reference(self) -> dict:
        return {"Phi": self.orbit.Phi, "L2": self.orbit.EL[1] ** 2}

    def errors(self, solution, formulation) -> dict:
        failed = _failed(solution)
        if failed:
            return failed
        ref = self.reference()
        L0 = angular_momentum_vector(*formulation.to_xp(solution.y[0]))
        L1 = angular_momentum_vector(*formulation.to_xp(solution.y_end))
        tilt = math.atan2(np.linalg.norm(np.cross(L0, L1)), float(L0 @ L1))
        psi = in_plane_angle(solution, formulation)
        L2 = formulation.invariants(solution.y)["L2"]
        return {
            "error": abs(psi[-1] - self.n_orbits * ref["Phi"]),
            "dpsi": psi[-1] - self.n_orbits * ref["Phi"],
            "tilt": tilt,
            "dL2_max": float(np.max(np.abs(L2 - ref["L2"])) / ref["L2"]),
            **self.constraint_errors(solution, formulation),
        }


FAMILIES = {
    cls.family: cls
    for cls in (
        Deflection,
        NearCritical,
        PhotonSphere,
        Circular,
        MarginalCircular,
        Eccentric,
        Inclined,
    )
}

#: The parameter grid of roadmap §1.3.
TC1_IMPACT_PARAMETERS = (5.3, 6.0, 8.0, 10.0, 20.0, 50.0, 100.0, 1000.0)
TC3_RADII = (6.0, 7.0, 10.0, 20.0, 100.0)
TC4_ORBITS = ((100.0, 0.5), (20.0, 0.5), (7.5, 0.5))
TC6_INCLINATIONS = (30.0, 60.0, 85.0)


def make_case(family: str, **params) -> GeodesicCase:
    return FAMILIES[family](**params)


def standard_cases() -> list[GeodesicCase]:
    """One instance of every test case of the matrix, short integrations."""
    cases: list[GeodesicCase] = [Deflection(b) for b in TC1_IMPACT_PARAMETERS]
    cases += [Deflection(b, 1e4) for b in TC1_IMPACT_PARAMETERS]
    cases += [NearCritical(k) for k in range(1, 11)]
    cases += [PhotonSphere(k) for k in (2, 4, 6, 8)]
    cases += [Circular(r) for r in TC3_RADII]
    cases += [Eccentric(p, e, 2.0) for p, e in TC4_ORBITS]
    cases += [Inclined(20.0, 0.5, i, 2.0) for i in TC6_INCLINATIONS]
    return cases
