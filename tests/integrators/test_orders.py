"""Measured convergence orders on TC0 (Phase 2 and 3 exit criteria: theory +- 0.1)."""

import numpy as np
import pytest

from geoint.integrators import get
from geoint.testcases.toy import kepler

KEPLER = kepler(0.5)
PERIOD = 2 * np.pi

# Step counts in the asymptotic range of each method: above the round-off floor and, for Tao,
# below the step sizes where omega h is too large.
CASES = {
    "RK4": (4, [512, 1024, 2048]),
    "GL1": (2, [256, 512, 1024]),
    "GL2": (4, [128, 256, 512]),
    "GL3": (6, [32, 64, 128]),
    "Tao2": (2, [512, 1024, 2048]),
    "Tao4": (4, [1024, 2048, 4096]),
}


def measured_order(name, steps):
    errors = []
    for n in steps:
        kwargs = {"omega": 20.0} if name.startswith("Tao") else {}
        sol = get(name).solve(KEPLER.problem(PERIOD), n_steps=n, save_every=n, **kwargs)
        errors.append(np.max(np.abs(sol.y_end - KEPLER.exact(PERIOD))))
    slope, _ = np.polyfit(np.log(PERIOD / np.array(steps)), np.log(errors), 1)
    return slope


@pytest.mark.parametrize("name", CASES)
def test_order_on_kepler(name):
    order, steps = CASES[name]
    assert get(name).order == order
    assert measured_order(name, steps) == pytest.approx(order, abs=0.1)


@pytest.mark.parametrize("name", ["DP5", "DOP853"])
def test_adaptive_error_is_tolerance_proportional(name):
    # The global error of a well-tuned pair scales roughly like the tolerance.
    errors = []
    for tol in (1e-7, 1e-9, 1e-11):
        sol = get(name).solve(KEPLER.problem(PERIOD), rtol=tol, atol=tol)
        errors.append(np.max(np.abs(sol.y_end - KEPLER.exact(PERIOD))))
    slope = np.polyfit(np.log([1e-7, 1e-9, 1e-11]), np.log(errors), 1)[0]
    assert 0.7 < slope < 1.2
