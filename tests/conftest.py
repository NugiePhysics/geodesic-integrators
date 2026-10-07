"""Shared fixtures. Metrics and formulations are session-scoped so their kernels compile once."""

import pytest

from geoint.metrics import Schwarzschild


@pytest.fixture(scope="session")
def schwarzschild():
    return Schwarzschild()
