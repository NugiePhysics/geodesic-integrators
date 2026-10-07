"""The validation table of the predecessor project ``schwarzschild-geodesics`` (README).

| Check | Numerical (old) | Analytical (old) |
| Circular orbit r_c = 10: dphi/dt | 0.0316227766 | 0.0316227766 |
| Periapsis advance r_apo = 20, L = 4.2 | 2.127094008 rad | 2.127093985 rad |
| Photon capture threshold | b = 5.19 captured, b = 5.20 escapes | b_c = 5.196 |
| Light deflection at b = 40 | 0.108104 rad | 0.108096 rad (4th-order series) |
| Global order of RK4 | 3.996 | 4 |

The old "analytical" precession, 2.127093985, was the elliptic formula evaluated at the
periapsis found as the minimum over the *sampled* numerical trajectory. The exact value is
2.12709400813 (closed form and 40-digit quadrature agree), which the old numerical value matched.
"""

import math

import numpy as np
import pytest

from geoint import analytic as an
from geoint.initial_conditions import turning_point_momentum
from geoint.integrators import Problem, get
from geoint.testcases.base import periapsis_event
from geoint.testcases.schwarzschild import Circular, Deflection, Eccentric


def test_circular_orbit_angular_velocity(hamiltonian):
    case = Circular(10.0, 5.0)
    sol = get("RK4").solve(case.problem(hamiltonian), h=1.0)
    assert case.errors(sol, hamiltonian)["omega"] == pytest.approx(0.0316227766, abs=1e-10)


def test_periapsis_advance(schwarzschild, hamiltonian):
    x0 = [0.0, 20.0, math.pi / 2, 0.0]
    y0 = hamiltonian.from_xp(x0, turning_point_momentum(schwarzschild, x0, 4.2, 1))
    problem = Problem(hamiltonian.rhs, y0, (0.0, 2000.0), (periapsis_event(),))
    sol = get("DOP853").solve(problem, rtol=1e-12, atol=1e-12)
    phi = sol.events["periapsis"].y[:, 3]
    advance = phi[1] - phi[0] - 2 * math.pi
    assert advance == pytest.approx(2.127094008, abs=1e-9)  # the old numerical value
    E = math.sqrt(0.9 * (1 + 4.2**2 / 400))
    exact = an.EccentricOrbit(*an.EL_to_pe(E, 4.2)).precession
    assert exact == pytest.approx(2.12709400813, abs=1e-11)
    assert advance == pytest.approx(exact, abs=1e-10)


@pytest.mark.parametrize(("b", "status"), [(5.19, "captured"), (5.20, "escaped")])
def test_capture_threshold(hamiltonian, b, status):
    case = Deflection(b)
    sol = get("DOP853").solve(case.problem(hamiltonian), rtol=1e-10, atol=1e-10)
    assert sol.status == status


def test_deflection_at_b40(hamiltonian):
    # The integration stops at r_far, so it is compared with Delta phi(r_far) (pitfall 18); the
    # asymptotic angle, which the old table quotes, then follows from the exact formula.
    case = Deflection(40.0, 1e4)
    sol = get("DOP853").solve(case.problem(hamiltonian), rtol=1e-12, atol=1e-12)
    assert case.errors(sol, hamiltonian)["error"] < 1e-11
    alpha = an.deflection_angle(40.0)
    assert alpha == pytest.approx(0.108104, abs=1e-6)  # the old numerical value
    assert alpha == pytest.approx(0.10810407654, abs=1e-11)
    assert an.weak_field_deflection(40.0) == pytest.approx(0.108096, abs=1e-6)


@pytest.mark.parametrize("formulation_name", ["second_order", "hamiltonian"])
def test_rk4_global_order(request, formulation_name):
    # Phase 2 exit criterion, on a strong-field orbit where RK4 is already asymptotic.
    form = request.getfixturevalue(formulation_name)
    case = Eccentric(7.5, 0.5, 1.0)
    steps = np.array([800, 1600, 3200])
    errors = [
        case.errors(get("RK4").solve(case.problem(form), n_steps=n, save_every=n), form)["error"]
        for n in steps
    ]
    order = -np.polyfit(np.log(steps), np.log(errors), 1)[0]
    assert order == pytest.approx(4.0, abs=0.1)
