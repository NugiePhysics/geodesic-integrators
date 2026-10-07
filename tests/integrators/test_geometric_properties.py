"""Symplecticity, symmetry, reversibility and long-time energy behaviour on TC0 (Phase 3)."""

import numpy as np
import pytest

from geoint.integrators import get
from geoint.integrators.tableaus import GL1, GL2, GL3
from geoint.integrators.tao import extend
from geoint.testcases.toy import harmonic, kepler, tao_example

KEPLER = kepler(0.5)
Y = np.array([0.6, 0.3, -0.4, 1.1])
H = 0.3
SYMPLECTIC = ["GL1", "GL2", "GL3", "Tao2", "Tao4"]


def J(n):
    m = n // 2
    return np.block([[np.zeros((m, m)), np.eye(m)], [-np.eye(m), np.zeros((m, m))]])


def jacobian(f, y, eps=2.5e-4):
    """Fourth-order central differences."""
    out = np.empty((y.size, y.size))
    for k in range(y.size):
        e = np.zeros(y.size)
        e[k] = eps
        out[:, k] = (-f(y + 2 * e) + 8 * f(y + e) - 8 * f(y - e) + f(y - 2 * e)) / (12 * eps)
    return out


def one_step(name):
    """The one-step map and its symplectic structure; Tao acts on (q, p, x, y)."""
    method = get(name)
    if name.startswith("Tao"):
        structure = np.block([[J(4), np.zeros((4, 4))], [np.zeros((4, 4)), J(4)]])
        return (lambda z, h=H: method.step(KEPLER.rhs, z, h, omega=20.0)), extend(Y), structure
    return (lambda y, h=H: method.step(KEPLER.rhs, y, h)), Y, J(4)


@pytest.mark.parametrize("name", SYMPLECTIC)
def test_symplectic(name):
    step, y, structure = one_step(name)
    P = jacobian(step, y)
    assert np.max(np.abs(P.T @ structure @ P - structure)) <= 1e-10


def test_rk4_is_not_symplectic():
    step, y, structure = one_step("RK4")
    P = jacobian(step, y)
    assert np.max(np.abs(P.T @ structure @ P - structure)) > 1e-4


@pytest.mark.parametrize("name", SYMPLECTIC)
def test_symmetric(name):
    step, y, _ = one_step(name)
    assert np.max(np.abs(step(step(y), -H) - y)) <= 1e-13


def test_rk4_is_not_symmetric():
    step, y, _ = one_step("RK4")
    assert np.max(np.abs(step(step(y), -H) - y)) > 1e-6


@pytest.mark.parametrize("name", ["GL2", "Tao4"])
def test_reversible_under_momentum_flip(name):
    # Integrate forward, flip the momenta, integrate forward again: back to the start.
    kwargs = {"omega": 20.0} if name.startswith("Tao") else {}
    method = get(name)
    forward = method.solve(KEPLER.problem(10.0), n_steps=200, **kwargs)
    flipped = forward.y_end * np.array([1, 1, -1, -1])
    if name.startswith("Tao"):
        # The extended state must be flipped as a whole: (q, -p, x, -y).
        z = forward.extended[-1] * np.array([1, 1, -1, -1, 1, 1, -1, -1])
        back = np.array(z)
        for _ in range(200):
            back = method.step(KEPLER.rhs, back, 10.0 / 200, omega=20.0)
        end = back[:4]
    else:
        end = method.solve(
            type(KEPLER)(KEPLER.name, KEPLER.rhs, flipped, KEPLER.energy).problem(10.0),
            n_steps=200,
        ).y_end
    np.testing.assert_allclose(end * np.array([1, 1, -1, -1]), KEPLER.y0, atol=1e-12)


@pytest.mark.parametrize(
    ("tableau", "s"), [(GL1, 1), (GL2, 2), (GL3, 3)], ids=["GL1", "GL2", "GL3"]
)
def test_gauss_tableau_is_symplectic_and_consistent(tableau, s):
    A, b, c = tableau.A, tableau.b, tableau.c
    M = b[:, None] * A + (b[:, None] * A).T - np.outer(b, b)  # GNI Thm VI.4.3
    assert np.max(np.abs(M)) <= 4e-17
    np.testing.assert_allclose(A.sum(axis=1), c, atol=1e-16)
    for k in range(1, 2 * s + 1):  # quadrature order 2s: sum b c^(k-1) = 1/k
        assert b @ c ** (k - 1) == pytest.approx(1 / k, abs=1e-15)


@pytest.mark.parametrize("name", ["GL1", "GL2", "Tao2", "Tao4"])
def test_energy_error_is_bounded(name):
    # Tao's non-separable example: the symplectic methods keep |dH| bounded (no drift).
    problem = tao_example()
    kwargs = {"omega": 50.0} if name.startswith("Tao") else {}
    sol = get(name).solve(problem.problem(1e3), n_steps=10**5, save_every=100, **kwargs)
    dH = np.abs(problem.energy(sol.y) - problem.energy(problem.y0))
    half = dH.size // 2
    assert dH[half:].max() <= 1.2 * dH[:half].max()


def test_rk4_energy_drifts():
    problem = tao_example()
    sol = get("RK4").solve(problem.problem(1e3), n_steps=10**5, save_every=100)
    dH = np.abs(problem.energy(sol.y) - problem.energy(problem.y0))
    half = dH.size // 2
    assert dH[half:].max() >= 1.8 * dH[:half].max()


@pytest.mark.slow
@pytest.mark.parametrize("name", SYMPLECTIC)
def test_energy_error_is_bounded_for_a_million_steps(name):
    # Phase 3 exit criterion.
    problem = tao_example()
    kwargs = {"omega": 50.0} if name.startswith("Tao") else {}
    sol = get(name).solve(problem.problem(1e4), n_steps=10**6, save_every=1000, **kwargs)
    dH = np.abs(problem.energy(sol.y) - problem.energy(problem.y0))
    half = dH.size // 2
    assert dH[half:].max() <= 1.2 * dH[:half].max()


def test_tao_conserves_the_extended_hamiltonian_structure():
    # Harmonic oscillator, omega large: the two copies stay together to O(h^2).
    sol = get("Tao2").solve(harmonic().problem(100.0), n_steps=10**4, omega=50.0)
    q, p, x, y = sol.extended.T
    assert np.max(np.abs(q - x)) < 1e-3 and np.max(np.abs(p - y)) < 1e-3
