"""Our compiled DP5 and DOP853 against SciPy (Phase 2 exit criterion: nfev within 1 %)."""

import numpy as np
import pytest

from geoint.initial_conditions import turning_point_momentum
from geoint.integrators import Problem, get
from geoint.testcases.toy import kepler

PAIRS = [("DP5", "scipy-RK45"), ("DOP853", "scipy-DOP853")]
TOLS = [1e-5, 1e-8, 1e-11]


def compare(problem, own, oracle, tol, atol=None):
    atol = tol if atol is None else atol
    a = get(own).solve(problem, rtol=tol, atol=atol)
    b = get(oracle).solve(problem, rtol=tol, atol=atol)
    assert a.status == b.status == "completed"
    assert abs(a.nfev - b.nfev) <= 0.01 * b.nfev
    assert a.n_steps == b.n_steps
    return a, b


@pytest.mark.parametrize(("own", "oracle"), PAIRS)
@pytest.mark.parametrize("tol", TOLS)
def test_matches_scipy_on_kepler(own, oracle, tol):
    a, b = compare(kepler(0.5).problem(20 * np.pi), own, oracle, tol)
    # Same steps, so the two solutions differ only by rounding (amplified along the orbit).
    np.testing.assert_allclose(a.y_end, b.y_end, rtol=0, atol=1e-9)


@pytest.mark.parametrize(("own", "oracle"), PAIRS)
@pytest.mark.parametrize("tol", TOLS)
def test_matches_scipy_on_a_geodesic(schwarzschild, hamiltonian, own, oracle, tol):
    # The r_apo = 20, L = 4.2 orbit of the predecessor's validation table, three radial periods.
    x0 = [0.0, 20.0, np.pi / 2, 0.0]
    y0 = hamiltonian.from_xp(x0, turning_point_momentum(schwarzschild, x0, 4.2, 1))
    a, b = compare(Problem(hamiltonian.rhs, y0, (0.0, 1500.0)), own, oracle, tol)
    np.testing.assert_allclose(a.y_end[1:4], b.y_end[1:4], rtol=0, atol=max(tol, 1e-9) * 10)
