"""Formulations (a) and (b) trace the same geodesics (Phase 1 exit criterion)."""

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from geoint.initial_conditions import momentum_from_constants, turning_point_momentum

EQUATOR = np.pi / 2
INCLINATION = np.radians(60)

# id: (x0, builder, keyword arguments, eps, lambda_end). The timelike orbit is the r_apo = 20,
# L = 4.2 case of the predecessor's validation table (about three radial periods); the photon
# passes at r_min = 4.45, well clear of the photon sphere, so the comparison is well conditioned.
ORBITS = {
    "eccentric-equatorial": (
        [0.0, 20.0, EQUATOR, 0.0],
        turning_point_momentum,
        dict(Lz=4.2, eps=1),
        1500.0,
    ),
    "eccentric-inclined-60deg": (
        [0.0, 20.0, EQUATOR, 0.0],
        turning_point_momentum,
        dict(Lz=4.2 * np.cos(INCLINATION), eps=1, p_theta=4.2 * np.sin(INCLINATION)),
        1500.0,
    ),
    "photon-b6": (
        [0.0, 1000.0, EQUATOR, 0.0],
        momentum_from_constants,
        dict(E=1.0, Lz=6.0, eps=0),
        1950.0,
    ),
}


@pytest.mark.parametrize("orbit", ORBITS)
def test_same_trajectory_in_both_formulations(schwarzschild, second_order, hamiltonian, orbit):
    x0, build, kwargs, lam_end = ORBITS[orbit]
    p0 = build(schwarzschild, x0, **kwargs)
    lam = np.linspace(0.0, lam_end, 2001)
    paths = []
    for form in (second_order, hamiltonian):
        sol = solve_ivp(
            lambda _, y, rhs=form.rhs: rhs(y),
            (0.0, lam_end),
            form.from_xp(x0, p0),
            method="DOP853",
            rtol=1e-13,
            atol=1e-13,
            dense_output=True,
        )
        assert sol.success
        paths.append(sol.sol(lam))
    a, b = paths
    for index in (1, 2, 3):  # r, theta, phi
        assert np.max(np.abs(a[index] - b[index])) <= 1e-9, index
