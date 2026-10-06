"""Smoke tests: the package imports and the Numba pattern the design relies on compiles."""

import importlib

import numpy as np
import pytest
from numba import njit

import geoint


def test_version_is_exposed():
    assert isinstance(geoint.__version__, str)


@pytest.mark.parametrize(
    "module", ["numpy", "scipy.integrate", "scipy.special", "matplotlib", "mpmath", "numba"]
)
def test_runtime_dependency_imports(module):
    importlib.import_module(module)


def make_decay_rhs(k):
    """Factory returning a jitted RHS that captures its parameter, as formulations will."""

    @njit
    def rhs(y):
        return -k * y

    return rhs


@njit
def euler(rhs, y0, h, n):
    y = y0.copy()
    for _ in range(n):
        y = y + h * rhs(y)
    return y


def test_jitted_loop_accepts_jitted_rhs_closure():
    # docs/decisions/0002: integrator loops are njit functions taking a jitted RHS as argument.
    k, h, n = 0.5, 0.01, 1000
    y0 = np.array([1.0, -2.0])
    y = euler(make_decay_rhs(k), y0, h, n)
    # Explicit Euler on y' = -k y is exactly y_n = (1 - h k)^n y_0.
    np.testing.assert_allclose(y, (1.0 - h * k) ** n * y0, rtol=1e-12)
