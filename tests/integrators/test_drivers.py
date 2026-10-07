"""Event location and the mechanics of the drivers (all integrators, toy problems)."""

import numpy as np
import pytest

from geoint.integrators import Event, Problem, get
from geoint.testcases.toy import harmonic, kepler

OSC = harmonic(1.0, 0.0)  # q = cos t, p = -sin t
# q crosses 0 upward at t = 3 pi/2 + 2 pi k, downward at pi/2 + 2 pi k.
UP = Event("q_up", 0, 0.0, +1)
DOWN = Event("q_down", 0, 0.0, -1)
METHODS = {
    "RK4": {"h": 0.01},
    "GL2": {"h": 0.05},
    "GL3": {"h": 0.1},
    "Tao4": {"h": 0.002, "omega": 50.0},
    "DP5": {"rtol": 1e-10, "atol": 1e-12},
    "DOP853": {"rtol": 1e-12, "atol": 1e-14},
    "scipy-DOP853": {"rtol": 1e-12, "atol": 1e-14},
}


@pytest.mark.parametrize("name", METHODS)
def test_event_times_and_directions(name):
    sol = get(name).solve(OSC.problem(20.0, events=(UP, DOWN)), **METHODS[name])
    up, down = sol.events["q_up"].lam, sol.events["q_down"].lam
    np.testing.assert_allclose(up, 1.5 * np.pi + 2 * np.pi * np.arange(3), atol=2e-8)
    np.testing.assert_allclose(down, 0.5 * np.pi + 2 * np.pi * np.arange(3), atol=2e-8)
    np.testing.assert_allclose(sol.events["q_up"].y[:, 0], 0.0, atol=1e-12)


@pytest.mark.parametrize("name", ["RK4", "GL3", "DOP853"])
def test_event_location_has_the_order_of_the_method(name):
    # Halving the step must shrink the event-time error like the global error, not like O(h^2).
    errors = []
    for k in (1, 2):
        options = dict(METHODS[name])
        if "h" in options:
            options["h"] *= 4 / k
        else:
            options = {"rtol": 1e-6 / 10**k, "atol": 1e-8 / 10**k}
        sol = get(name).solve(OSC.problem(5.0, events=(UP,)), **options)
        errors.append(abs(sol.events["q_up"].lam[0] - 1.5 * np.pi))
    if name == "DOP853":
        assert errors[1] < errors[0]
    else:
        assert errors[0] / errors[1] > 0.7 * 2 ** get(name).order


@pytest.mark.parametrize("name", ["RK4", "GL2", "DP5", "scipy-RK45"])
def test_terminal_event_and_count(name):
    stop_second = Event("q_up", 0, 0.0, +1, stop_after=2)
    options = METHODS.get(name, {"rtol": 1e-10, "atol": 1e-12})
    sol = get(name).solve(OSC.problem(100.0, events=(stop_second,)), **options)
    assert sol.status == "q_up"
    assert sol.lam_end == pytest.approx(1.5 * np.pi + 2 * np.pi, abs=1e-7)
    assert len(sol.events["q_up"].lam) == 2
    np.testing.assert_allclose(sol.y_end, sol.events["q_up"].y[-1])


def test_start_on_the_surface_is_not_an_event():
    # p starts at exactly 0 (a turning point); the first upward crossing of p is at t = pi.
    sol = get("RK4").solve(OSC.problem(4.0, events=(Event("p_up", 1, 0.0, +1),)), h=0.01)
    np.testing.assert_allclose(sol.events["p_up"].lam, [np.pi], atol=1e-8)


def test_fixed_step_parameter_is_n_times_h():
    sol = get("RK4").solve(OSC.problem(1.0), n_steps=7)
    h = 1.0 / 7
    assert sol.options["h"] == h
    np.testing.assert_array_equal(sol.lam, np.arange(8) * h)


def test_h_is_shrunk_to_divide_the_interval():
    sol = get("RK4").solve(OSC.problem(1.0), h=0.3)
    assert sol.n_steps == 4 and sol.options["h"] == 0.25


def test_save_every_keeps_the_end_point():
    sol = get("GL2").solve(OSC.problem(1.0), n_steps=10, save_every=4)
    np.testing.assert_allclose(sol.lam, [0.0, 0.4, 0.8, 1.0])
    full = get("GL2").solve(OSC.problem(1.0), n_steps=10)
    np.testing.assert_array_equal(sol.y_end, full.y_end)


def test_compensated_summation_reduces_round_off():
    # 10^6 tiny steps of a high-order method: the error is pure round-off accumulation.
    problem = OSC.problem(10.0)
    exact = OSC.exact(10.0)
    errors = {}
    for comp in (False, True):
        sol = get("GL3").solve(problem, n_steps=10**6, save_every=10**6, compensated=comp)
        errors[comp] = np.max(np.abs(sol.y_end - exact))
    assert errors[True] < 0.1 * errors[False]


def test_iteration_tolerance_stops_earlier():
    problem = kepler(0.5).problem(2 * np.pi)
    strict = get("GL2").solve(problem, n_steps=200)
    loose = get("GL2").solve(problem, n_steps=200, iter_tol=1e-6)
    assert loose.n_iter < strict.n_iter
    assert strict.n_iter / strict.n_steps < 15


def test_divergent_implicit_iteration_is_reported():
    # Midpoint rule on the oscillator: the fixed-point map has spectral radius h/2 = 2.5.
    sol = get("GL1").solve(OSC.problem(10.0), n_steps=2)
    assert sol.status == "failed: implicit iteration did not converge"
    assert not sol.success


def test_adaptive_failure_is_reported():
    sol = get("DOP853").solve(kepler(0.5).problem(2 * np.pi), rtol=1e-10, atol=1e-10, max_steps=5)
    assert sol.status == "failed: maximum number of steps"


def test_problem_validation():
    with pytest.raises(ValueError, match="forward"):
        Problem(OSC.rhs, OSC.y0, (1.0, 0.0))
    with pytest.raises(ValueError, match="component"):
        Problem(OSC.rhs, OSC.y0, (0.0, 1.0), events=(Event("bad", 5, 0.0),))
    with pytest.raises(ValueError, match="exactly one"):
        get("RK4").solve(OSC.problem(1.0))
    with pytest.raises(ValueError, match="Hamiltonian"):
        get("Tao2").solve(Problem(OSC.rhs, OSC.y0, (0.0, 1.0)), h=0.1, omega=1.0)


def test_tao_reports_the_first_copy_and_counts_three_evaluations_per_step():
    sol = get("Tao2").solve(OSC.problem(1.0), n_steps=100, omega=50.0)
    assert sol.y.shape[1] == 2 and sol.extended.shape[1] == 4
    np.testing.assert_array_equal(sol.y, sol.extended[:, :2])
    # 4 gradient flows per Strang step, but the closing phi_A of one step and the opening phi_A
    # of the next are evaluated at the same (q, y): 3 per step plus the very first one.
    assert sol.nfev == 3 * 100 + 1
    assert get("Tao4").solve(OSC.problem(1.0), n_steps=100, omega=50.0).nfev == 9 * 100 + 1
