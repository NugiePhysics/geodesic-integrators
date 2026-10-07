"""Shared fixtures. Metrics and formulations are session-scoped so their kernels compile once."""

import pytest

from geoint.formulations import Hamiltonian, SecondOrder
from geoint.metrics import Schwarzschild


@pytest.fixture(scope="session")
def schwarzschild():
    return Schwarzschild()


@pytest.fixture(scope="session")
def second_order(schwarzschild):
    return SecondOrder(schwarzschild)


@pytest.fixture(scope="session")
def hamiltonian(schwarzschild):
    return Hamiltonian(schwarzschild)


@pytest.fixture(scope="session", params=["second_order", "hamiltonian"])
def formulation(request):
    return request.getfixturevalue(request.param)
