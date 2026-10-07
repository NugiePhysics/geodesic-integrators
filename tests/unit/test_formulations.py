"""Right-hand sides and invariants of formulations (a) and (b), and how they relate."""

import mpmath
import numpy as np
import pytest

import derive_metric
from helpers import random_positions

RHO = np.array([1, 1, 1, 1, -1, -1, -1, -1], dtype=float)  # flips velocity or momentum


@pytest.fixture(scope="module")
def dg_lower():
    """``d_a g_{mu nu}`` from SymPy, independent of the kernels under test.

    Evaluated at 30 digits: in float64 the SymPy form ``1 - 2M/r`` would itself lose digits near
    the horizon, which is exactly what the hand-written kernels avoid.
    """
    sym = derive_metric.schwarzschild()
    fn = derive_metric.lambdify(sym, derive_metric.derive(sym).dg)

    def dg(x):
        with mpmath.workdps(30):
            return np.array(fn([mpmath.mpf(v) for v in x]), dtype=float)

    return dg


def random_states(rng, n, **kwargs):
    """Positions as in ``random_positions`` with an O(1) velocity/momentum, not on the shell."""
    return np.hstack([random_positions(rng, n, **kwargs), rng.normal(size=(n, 4))])


def test_second_order_rhs_matches_explicit_equations(second_order):
    # Theory eq. (8): the generic Gamma contraction against the hand-expanded geodesic equation.
    for y in random_states(np.random.default_rng(10), 100, min_gap=1e-3):
        r, th = y[1], y[2]
        ut, ur, uth, uph = y[4:]
        f, df = 1 - 2 / r, 2 / r**2
        expected = [
            -df / f * ut * ur,
            -0.5 * f * df * ut**2
            + df / (2 * f) * ur**2
            + r * f * (uth**2 + np.sin(th) ** 2 * uph**2),
            -2 / r * ur * uth + np.sin(th) * np.cos(th) * uph**2,
            -2 / r * ur * uph - 2 / np.tan(th) * uth * uph,
        ]
        dy = second_order.rhs(y)
        assert np.array_equal(dy[:4], y[4:])
        np.testing.assert_allclose(dy[4:], expected, rtol=1e-12, atol=1e-14)


def test_formulations_agree_through_p_equals_g_u(
    schwarzschild, second_order, hamiltonian, dg_lower
):
    # Theory eq. (17): for any (x, u), on or off the shell, with p = g u,
    # dx/dlambda of (b) is u, and dp/dlambda of (b) is d_rho g_{mu nu} u^rho u^nu + g_{mu nu} du^nu.
    for y in random_states(np.random.default_rng(11), 200, min_gap=1e-6):
        x, u = y[:4], y[4:]
        g = schwarzschild.g(x)
        dyb = hamiltonian.rhs(hamiltonian.from_xp(x, g @ u))
        du = second_order.rhs(y)[4:]
        dg = dg_lower(x)  # [rho, mu, nu]
        first = np.einsum("rmn,r,n->rmn", dg, u, u)  # terms of d_rho g_{mu nu} u^rho u^nu, per mu
        second = g * du  # terms of g_{mu nu} du^nu
        expected = first.sum(axis=(0, 2)) + second.sum(axis=1)
        scale = np.abs(first).sum(axis=(0, 2)) + np.abs(second).sum(axis=1)
        np.testing.assert_allclose(dyb[:4], u, rtol=4 * np.finfo(float).eps, atol=0)
        assert np.all(np.abs(dyb[4:] - expected) <= 1e-14 * scale), (y, dyb[4:], expected)


def test_cyclic_momenta_have_exactly_zero_rate(hamiltonian):
    # Theory §5.2: the p_t and p_phi increments of every RK stage are floating-point zeros.
    for y in random_states(np.random.default_rng(12), 200):
        dy = hamiltonian.rhs(y)
        assert dy[4] == 0.0 and dy[7] == 0.0
        dH_dx = hamiltonian.dH_dx(y[:4], y[4:])
        assert dH_dx[0] == 0.0 and dH_dx[3] == 0.0


def test_hamiltonian_rhs_is_the_split_gradient(hamiltonian):
    for y in random_states(np.random.default_rng(13), 50):
        x, p = y[:4], y[4:]
        expected = np.concatenate([hamiltonian.dH_dp(x, p), -hamiltonian.dH_dx(x, p)])
        assert np.array_equal(hamiltonian.rhs(y), expected)


def test_hamiltonian_value(hamiltonian):
    # Theory eq. (14).
    x = np.array([0.0, 7.0, 1.1, 0.3])
    p = np.array([-0.95, 0.2, 1.5, 3.0])
    f = 1 - 2 / 7.0
    expected = 0.5 * (-(0.95**2) / f + f * 0.2**2 + 1.5**2 / 49 + 3.0**2 / (49 * np.sin(1.1) ** 2))
    assert hamiltonian.invariants(hamiltonian.from_xp(x, p))["H"] == pytest.approx(expected, 1e-15)


def test_reversibility_is_exact(formulation):
    # Theory eq. (23): rho F(y) = -F(rho y), with rho flipping the velocity or momentum.
    for y in random_states(np.random.default_rng(14), 100):
        assert np.array_equal(RHO * formulation.rhs(y), -formulation.rhs(RHO * y))


def test_xp_round_trip(formulation):
    for y in random_states(np.random.default_rng(15), 100, min_gap=1e-3):
        x, p = y[:4], y[4:]
        x2, p2 = formulation.to_xp(formulation.from_xp(x, p))
        assert np.array_equal(x2, x)
        np.testing.assert_allclose(p2, p, rtol=4 * np.finfo(float).eps)


def test_invariants_agree_between_formulations(second_order, hamiltonian):
    for y in random_states(np.random.default_rng(16), 100, min_gap=1e-3):
        x, p = y[:4], y[4:]
        a = second_order.invariants(second_order.from_xp(x, p))
        b = hamiltonian.invariants(hamiltonian.from_xp(x, p))
        assert a.keys() == b.keys() == {"H", "E", "Lz", "L2"}
        for name in ("E", "Lz", "L2"):
            assert a[name] == pytest.approx(b[name], rel=1e-15, abs=1e-300)
        # H is a sum of terms of both signs; compare on the scale of the largest one.
        scale = np.abs(np.diag(hamiltonian.metric.g_inv(x)) * p**2).max()
        assert abs(a["H"] - b["H"]) <= 1e-15 * scale


def test_invariants_of_a_trajectory(formulation):
    Y = random_states(np.random.default_rng(17), 5)
    batch = formulation.invariants(Y)
    assert set(batch) == set(formulation.invariant_names)
    for i, y in enumerate(Y):
        single = formulation.invariants(y)
        for name in formulation.invariant_names:
            assert batch[name].shape == (5,)
            assert np.ndim(single[name]) == 0
            assert single[name] == batch[name][i]


def test_relative_constraint_normalizations(hamiltonian):
    x = np.array([0.0, 10.0, np.pi / 2, 0.0])
    f = 0.8
    # Timelike, H = -1/2 + d  ->  delta H = 2 d.  Null with E = 2: delta H = |H| / 4.
    for eps, E, H_target, expected in [(1, 1.0, -0.5 + 1e-6, 2e-6), (0, 2.0, 3e-6, 3e-6 / 4)]:
        p_r = np.sqrt((2 * H_target + E**2 / f) / f)
        y = hamiltonian.from_xp(x, [-E, p_r, 0.0, 0.0])
        assert hamiltonian.relative_constraint(y, eps) == pytest.approx(expected, rel=1e-8)
        assert hamiltonian.constraint(y, eps) == pytest.approx(H_target + eps / 2, rel=1e-8)
