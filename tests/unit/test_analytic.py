"""The analytic references, each checked along two independent paths and in its limits."""

import math

import numpy as np
import pytest

from geoint import analytic as an

IMPACT = np.geomspace(5.2, 1e4, 15)


@pytest.mark.parametrize("b", IMPACT)
@pytest.mark.parametrize("r_far", [math.inf, 1e3, 1e4])
def test_deflection_closed_form_matches_quadrature(b, r_far):
    # Phase 4 exit criterion: <= 1e-12 relative for b in [5.2, 1e4].
    if an.periapsis(b) >= r_far:
        pytest.skip("periapsis outside r_far")
    exact = an.delta_phi(b, r_far)
    assert an.delta_phi_quadrature(b, r_far) == pytest.approx(exact, rel=1e-12)


@pytest.mark.parametrize("b", IMPACT)
def test_darwin_form_matches(b):
    assert an.deflection_darwin(b) == pytest.approx(an.deflection_angle(b), rel=1e-12)


def test_weak_field_limit():
    # The fourth-order series is off by the fifth-order term, O(b^-5).
    for b in (1e2, 1e3, 1e4):
        alpha = an.deflection_angle(b)
        assert abs(alpha - an.weak_field_deflection(b, 4)) < 1e3 / b**5
        assert alpha == pytest.approx(4 / b, rel=5 / b)


def test_strong_deflection_limit():
    # Bozza (2002): the error of the logarithmic limit vanishes like b - b_c.
    errors = []
    for k in (4, 6, 8, 10):
        b = an.B_CRIT * (1 + 10.0**-k)
        errors.append(abs(an.deflection_angle(b) - an.bozza_deflection(b)))
    ratios = np.array(errors[:-1]) / np.array(errors[1:])
    assert np.all((ratios > 50) & (ratios < 200))
    assert an.BOZZA_B == pytest.approx(-0.4002, abs=1e-4)


def test_photon_capture_threshold():
    with pytest.raises(ValueError, match="captured"):
        an.delta_phi(5.19)
    assert an.periapsis(5.2) == pytest.approx(3.0, abs=0.1)


@pytest.mark.parametrize(("p", "e"), [(100, 0.5), (20, 0.5), (7.5, 0.5), (13.0, 0.1)])
def test_precession_formula_matches_quadrature_and_orbit_shape(p, e):
    orbit = an.EccentricOrbit(p, e)
    # Phi from the elliptic integral and from the independent 30-digit quadrature.
    assert orbit.Phi == pytest.approx(orbit.periods["Phi"], rel=1e-14)
    # r(phi) from Jacobi functions: periapsis at 0 and Phi, apoapsis at Phi/2, and in between
    # the orbit equation holds.
    assert orbit.r_of_phi(0.0) == pytest.approx(orbit.r_p, rel=1e-15)
    assert orbit.r_of_phi(orbit.Phi) == pytest.approx(orbit.r_p, rel=1e-13)
    assert orbit.r_of_phi(orbit.Phi / 2) == pytest.approx(orbit.r_a, rel=1e-13)
    E, L = orbit.EL
    phi = np.linspace(0.1, orbit.Phi - 0.1, 7)
    h = 1e-5
    u = 1 / orbit.r_of_phi(phi)
    du = (1 / orbit.r_of_phi(phi + h) - 1 / orbit.r_of_phi(phi - h)) / (2 * h)
    rhs = (E**2 - (1 - 2 * u) * (1 + L**2 * u**2)) / L**2
    np.testing.assert_allclose(du**2, rhs, rtol=1e-8, atol=1e-12)


def test_weak_field_precession_limit():
    orbit = an.EccentricOrbit(1e4, 0.5)
    assert orbit.precession == pytest.approx(6 * math.pi / 1e4, rel=1e-3)


def test_pe_round_trip_and_separatrix():
    for p, e in [(20, 0.5), (7.5, 0.5), (100, 0.1)]:
        assert an.EL_to_pe(*an.pe_to_EL(p, e)) == pytest.approx((p, e), rel=1e-10)
    with pytest.raises(ValueError, match="bound"):
        an.pe_to_EL(6.9, 0.5)
    assert an.separatrix(0.5) == 7.0


def test_circular_orbits():
    assert an.circular_orbit_constants(6.0) == pytest.approx((math.sqrt(8 / 9), 2 * math.sqrt(3)))
    assert an.circular_orbit_frequency(10.0) == pytest.approx(10**-1.5)
    T_tau, T_t = an.circular_orbit_period(10.0)
    assert T_t / T_tau == pytest.approx(an.circular_orbit_constants(10.0)[0] / 0.8)


def test_photon_sphere_turning_point_integral():
    # The one-way angle from a turning point is half the symmetric deflection integral.
    r0, r_exit = 4.0, 50.0
    b = math.sqrt(r0**3 / (r0 - 2))
    assert an.turning_point_delta_phi(r0, r_exit) == pytest.approx(
        an.delta_phi(b, r_exit) / 2, rel=1e-12
    )
    assert an.LYAPUNOV_PHOTON_SPHERE == pytest.approx(1 / (3 * math.sqrt(3)))
