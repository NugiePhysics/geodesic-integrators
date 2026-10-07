"""Fits and interpolations used by the figure scripts."""

import numpy as np
import pytest

from geoint.experiments.analysis import (
    classify_growth,
    cost_at_error,
    envelope,
    error_at_cost,
    fit_slope,
    pareto_front,
)


def test_fit_slope_ignores_the_floor_and_the_pre_asymptotic_regime():
    h = np.geomspace(1e-3, 1, 13)
    err = np.maximum(3.0 * h**4, 1e-13)  # round-off floor at small h
    err[-2:] = 1e-1  # pre-asymptotic points at large h
    slope, se, n = fit_slope(h, err, lo=1e-12, hi=1e-3)
    assert slope == pytest.approx(4.0, abs=1e-10) and n >= 3 and se < 1e-8


def test_work_precision_interpolation():
    cost = np.array([10.0, 100.0, 1000.0, 50.0])
    error = np.array([1e-2, 1e-4, 1e-6, 1e-1])  # the last point is dominated
    c, e = pareto_front(cost, error)
    np.testing.assert_array_equal(c, [10.0, 100.0, 1000.0])
    assert cost_at_error(cost, error, 1e-5) == pytest.approx(10**2.5)
    assert np.isnan(cost_at_error(cost, error, 1e-8))
    assert error_at_cost(cost, error, 10**1.5) == pytest.approx(1e-3)


def test_envelope_drops_a_short_trailing_window():
    x = np.arange(0, 10.01, 0.1)  # ten windows of width 1 plus the single end point
    lam, env = envelope(x, np.sin(x) ** 2, 1.0)
    assert lam.size == 10 and np.all(env <= 1)


@pytest.mark.parametrize(
    ("p", "law"), [(0.05, "bounded"), (0.45, "random walk"), (1.1, "linear"), (1.8, "quadratic")]
)
def test_growth_classification(p, law):
    assert classify_growth(p) == law
