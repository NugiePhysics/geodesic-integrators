"""Hand-written Schwarzschild kernels against SymPy, finite differences and their own identities."""

import mpmath
import numpy as np
import pytest

import derive_metric
from geoint.metrics import Schwarzschild, StopSurface
from helpers import random_positions

KERNELS = ("g", "g_inv", "dg_inv", "christoffel")


@pytest.fixture(scope="module", params=[1.0, 2.5], ids=lambda M: f"M={M}")
def metric_and_reference(request):
    M = request.param
    sym = derive_metric.schwarzschild(M)
    derived = derive_metric.derive(sym).arrays()
    reference = {name: derive_metric.lambdify(sym, derived[name]) for name in KERNELS}
    return Schwarzschild(M), reference


def test_kernels_match_sympy_to_1e14(metric_and_reference):
    # Phase 1 exit criterion. The reference is evaluated at 50 digits on the same float64 inputs,
    # so the comparison measures the kernels alone, down to r - 2M = 1e-8 M.
    metric, reference = metric_and_reference
    for x in random_positions(np.random.default_rng(7), 100, M=metric.M):
        with mpmath.workdps(50):
            mx = [mpmath.mpf(v) for v in x]
            for name in KERNELS:
                ref = np.array(reference[name](mx), dtype=float)
                num = getattr(metric, name)(x)
                zero = ref == 0
                assert np.all(num[zero] == 0.0), name
                rel = np.abs(num[~zero] - ref[~zero]) / np.abs(ref[~zero])
                assert rel.max() <= 1e-14, (name, x, rel.max())


def test_inverse(schwarzschild):
    for x in random_positions(np.random.default_rng(1), 200):
        product = schwarzschild.g(x) @ schwarzschild.g_inv(x)
        np.testing.assert_allclose(product, np.eye(4), rtol=0, atol=4 * np.finfo(float).eps)


def test_symmetries_are_exact(schwarzschild):
    for x in random_positions(np.random.default_rng(2), 50):
        for name in ("g", "g_inv"):
            a = getattr(schwarzschild, name)(x)
            assert np.array_equal(a, a.T)
        dgi = schwarzschild.dg_inv(x)
        assert np.array_equal(dgi, dgi.transpose(0, 2, 1))
        gamma = schwarzschild.christoffel(x)
        assert np.array_equal(gamma, gamma.transpose(0, 2, 1))


def test_cyclic_coordinates_have_exactly_zero_derivatives(schwarzschild):
    # Theory §5.2: exact zeros here are what make p_t and p_phi conserved bit for bit in (b).
    assert schwarzschild.cyclic_coords == (0, 3)
    for x in random_positions(np.random.default_rng(3), 50):
        dgi = schwarzschild.dg_inv(x)
        for mu in schwarzschild.cyclic_coords:
            assert not np.any(dgi[mu])


def test_dg_inv_matches_central_differences(schwarzschild):
    rng = np.random.default_rng(4)
    for x in random_positions(rng, 50, min_gap=0.5, max_gap=1e2, pole_gap=0.1):
        exact = schwarzschild.dg_inv(x)
        for a in range(4):
            h = 1e-5 * max(1.0, abs(x[a]))
            e = np.zeros(4)
            e[a] = h
            fd = (
                -schwarzschild.g_inv(x + 2 * e)
                + 8 * schwarzschild.g_inv(x + e)
                - 8 * schwarzschild.g_inv(x - e)
                + schwarzschild.g_inv(x - 2 * e)
            ) / (12 * h)
            scale = np.abs(schwarzschild.g_inv(x)).max() / max(1.0, abs(x[a]))
            np.testing.assert_allclose(exact[a], fd, rtol=0, atol=1e-8 * scale)


def test_invariants_in_canonical_variables(schwarzschild):
    x = np.array([0.0, 7.0, 1.1, 0.3])
    p = np.array([-0.95, 0.2, 1.5, 3.0])
    E, Lz, L2 = schwarzschild.invariants(x, p)
    assert (E, Lz) == (0.95, 3.0)
    assert L2 == pytest.approx(1.5**2 + 3.0**2 / np.sin(1.1) ** 2, rel=1e-15)
    assert schwarzschild.invariant_names == ("E", "Lz", "L2")


def test_horizon_and_stop_surface():
    metric = Schwarzschild(M=1.5)
    assert metric.horizon == 3.0
    (surface,) = metric.stop_surfaces(horizon_margin=1e-2)
    assert surface == StopSurface("captured", index=1, value=3.0 * 1.01, direction=-1)


def test_kernels_are_cached_per_mass():
    assert Schwarzschild(1.0).kernels is Schwarzschild(1.0).kernels
    assert Schwarzschild(1.0).kernels is not Schwarzschild(2.0).kernels


@pytest.mark.parametrize("M", [0.0, -1.0, float("nan")])
def test_rejects_non_positive_mass(M):
    with pytest.raises(ValueError, match="positive"):
        Schwarzschild(M)
