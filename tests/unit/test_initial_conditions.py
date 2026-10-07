"""Initial data on the mass shell: canonical orbits, random orbits, and the root solver itself."""

import mpmath
import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from geoint.initial_conditions import (
    momentum_from_constants,
    solve_quadratic_shell,
    turning_point_momentum,
)

EQUATOR = np.pi / 2
ULP = np.finfo(float).eps


def circular(rc):
    """E and L of the circular timelike orbit at ``rc`` (theory §9)."""
    return (1 - 2 / rc) / np.sqrt(1 - 3 / rc), np.sqrt(rc) / np.sqrt(1 - 3 / rc)


def eccentric(p, e):
    """E and L for semi-latus rectum ``p`` and eccentricity ``e`` (Cutler et al. 1994)."""
    E = np.sqrt(((p - 2) ** 2 - 4 * e**2) / (p * (p - 3 - e**2)))
    return E, p / np.sqrt(p - 3 - e**2)


def canonical_cases():
    """(id, x, builder, keyword arguments) for the orbits of the test-case matrix.

    The matrix is listed in docs/testcases.md.
    """
    cases = []
    for rc in (6, 7, 10, 20, 100):
        x = [0.0, rc, EQUATOR, 0.0]
        cases.append((f"circular-{rc}", x, turning_point_momentum, dict(Lz=circular(rc)[1], eps=1)))
    for p, e in ((100, 0.5), (20, 0.5), (7.5, 0.5)):
        E, L = eccentric(p, e)
        for r in (p / (1 + e), p, p / (1 - e)):
            x = [0.0, r, EQUATOR, 0.0]
            kwargs = dict(E=E, Lz=L, eps=1)
            cases.append((f"eccentric-{p}-{e}-r{r:g}", x, momentum_from_constants, kwargs))
        for inclination in (30, 60, 85):
            i = np.radians(inclination)
            x = [0.0, p, EQUATOR, 0.0]
            kwargs = dict(E=E, Lz=L * np.cos(i), eps=1, p_theta=L * np.sin(i))
            cases.append((f"inclined-{p}-{inclination}deg", x, momentum_from_constants, kwargs))
    for b in (5.3, 6, 8, 10, 20, 50, 100, 1e3):
        for r_far in (1e3, 1e4):
            x = [0.0, r_far, EQUATOR, 0.0]
            kwargs = dict(E=1.0, Lz=b, eps=0)
            cases.append((f"null-b{b:g}-r{r_far:g}", x, momentum_from_constants, kwargs))
    return cases


CASES = canonical_cases()


@pytest.mark.parametrize(
    ("x", "build", "kwargs"), [c[1:] for c in CASES], ids=[c[0] for c in CASES]
)
def test_canonical_initial_constraint_below_1e15(schwarzschild, formulation, x, build, kwargs):
    # Phase 1 exit criterion: relative constraint (theory §5.1) <= 1e-15 in both formulations.
    p = build(schwarzschild, x, **kwargs)
    y = formulation.from_xp(x, p)
    assert formulation.relative_constraint(y, kwargs["eps"]) <= 1e-15


@pytest.mark.parametrize("rc", [6, 7, 10, 20, 100])
def test_circular_orbit_from_turning_point(schwarzschild, rc):
    E, L = circular(rc)
    p = turning_point_momentum(schwarzschild, [0.0, rc, EQUATOR, 0.0], L, 1)
    assert p[1] == 0.0 and p[3] == L
    assert -p[0] == pytest.approx(E, rel=2 * ULP)


def test_isco_constants(schwarzschild):
    p = turning_point_momentum(schwarzschild, [0.0, 6.0, EQUATOR, 0.0], 2 * np.sqrt(3), 1)
    assert -p[0] == pytest.approx(np.sqrt(8 / 9), rel=2 * ULP)


def test_radial_sign_selects_direction(schwarzschild, hamiltonian):
    x = np.array([0.0, 30.0, EQUATOR, 0.0])
    for sign in (-1, 1):
        p = momentum_from_constants(schwarzschild, x, 1.0, 6.0, 0, radial_sign=sign)
        assert np.sign(hamiltonian.rhs(hamiltonian.from_xp(x, p))[1]) == sign


def test_forbidden_region_raises(schwarzschild):
    # E = 0.9 with L = 4 is bound well inside r = 100: E^2 = 0.81 < V_eff(100) ~ 0.98.
    with pytest.raises(ValueError, match="forbidden"):
        momentum_from_constants(schwarzschild, [0.0, 100.0, EQUATOR, 0.0], 0.9, 4.0, 1)


@settings(max_examples=300, deadline=None)
@given(
    gap=st.floats(-6, 4),  # log10(r - 2)
    theta=st.floats(0.05, np.pi - 0.05),
    L=st.floats(0, 50),
    p_theta_share=st.floats(0, 1),
    excess=st.floats(-8, 2),  # log10 of (E^2 - V_eff) / V_eff
    eps=st.sampled_from([0, 1]),
    sign=st.sampled_from([-1, 1]),
)
def test_random_initial_data_is_on_the_shell(
    schwarzschild, second_order, hamiltonian, gap, theta, L, p_theta_share, excess, eps, sign
):
    # |C| is bounded by a few units of rounding in the largest term of H, wherever we start.
    r = 2.0 + 10.0**gap
    if eps == 0:
        L = max(L, 1e-3)  # a null geodesic needs some angular momentum to have a finite E here
    f = (r - 2) / r
    V = f * (eps + L**2 / r**2)
    E = np.sqrt(V * (1 + 10.0**excess))
    Lz, p_theta = L * np.sqrt(1 - p_theta_share) * np.sin(theta), L * np.sqrt(p_theta_share)
    x = np.array([0.0, r, theta, 0.0])
    p = momentum_from_constants(schwarzschild, x, E, Lz, eps, p_theta=p_theta, radial_sign=sign)
    scale = 0.5 * np.abs(np.diag(schwarzschild.g_inv(x)) * p**2).sum()
    for form in (second_order, hamiltonian):
        y = form.from_xp(x, p)
        assert abs(form.constraint(y, eps)) <= 16 * ULP * scale
        inv = form.invariants(y)
        assert inv["E"] == pytest.approx(E, rel=4 * ULP)
        assert inv["Lz"] == pytest.approx(Lz, rel=4 * ULP, abs=1e-300)
    # In the canonical variables the conserved momenta are set, not computed.
    assert (p[0], p[3]) == (-E, Lz)


class TestQuadraticSolver:
    """``solve_quadratic_shell`` on synthetic forms, including off-diagonal ones (Kerr-like)."""

    @staticmethod
    def random_lorentzian(rng):
        boost = np.eye(4) + 0.3 * rng.normal(size=(4, 4))
        return boost.T @ np.diag([-1.0, 1.0, 1.0, 1.0]) @ boost

    def test_both_roots_on_off_diagonal_forms(self):
        rng = np.random.default_rng(20)
        solved = 0
        for _ in range(200):
            G, v, eps = self.random_lorentzian(rng), rng.normal(size=4), rng.choice([0.0, 1.0])
            k = int(rng.integers(4))
            for sign in (-1, 1):
                try:
                    w = solve_quadratic_shell(G, v, k, eps, sign)
                except ValueError:
                    continue
                solved += 1
                assert np.array_equal(np.delete(w, k), np.delete(v, k))
                assert np.sign((G @ w)[k]) == sign
                terms = np.abs(G * np.outer(w, w)).sum()
                assert abs(w @ G @ w + eps) <= 32 * ULP * terms
        assert solved >= 100  # about half the random forms have real roots; not vacuous

    def test_small_root_without_cancellation(self):
        # A v^2 + 2 B v + C = 0 with A << B: the small root is about -C / (2B), and the textbook
        # formula (-B + sqrt(B^2 - AC)) / A would lose about eight digits to cancellation.
        A, B, C = 1e-8, 1.0, 1.0
        G = np.array([[A, B], [B, 0.0]])
        w = solve_quadratic_shell(G, [0.0, 1.0], 0, C, +1)
        with mpmath.workdps(50):
            exact = (-B + mpmath.sqrt(B * B - mpmath.mpf(A) * C)) / A
        assert w[0] == pytest.approx(float(exact), rel=2 * ULP)

    def test_linear_case(self):
        # G_kk = 0 (null coordinate): the equation is linear and only the root with sign(B) exists.
        G = np.array([[0.0, 1.0], [1.0, 0.0]])
        assert solve_quadratic_shell(G, [0.0, 2.0], 0, 1.0, +1)[0] == -0.25
        with pytest.raises(ValueError, match="infinity"):
            solve_quadratic_shell(G, [0.0, 2.0], 0, 1.0, -1)

    def test_turning_point_rounding_is_clamped(self):
        G = np.diag([-1.0, 1.0])
        # Exactly at the turning point up to one ulp: D = -2.2e-16 is rounding, so p_r = 0.
        assert solve_quadratic_shell(G, [1.0, 0.0], 1, 1.0 + ULP, -1)[1] == 0.0
        with pytest.raises(ValueError, match="forbidden"):
            solve_quadratic_shell(G, [1.0, 0.0], 1, 1.0 + 1e-6, -1)

    def test_rejects_bad_sign(self):
        with pytest.raises(ValueError, match="sign"):
            solve_quadratic_shell(np.eye(2), [0.0, 1.0], 0, 0.0, 0)
