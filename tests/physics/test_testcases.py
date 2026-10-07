"""Every test case against its exact reference with a tight DOP853 (Phase 4 exit criterion)."""

import math

import numpy as np
import pytest

from geoint.integrators import get
from geoint.testcases.schwarzschild import (
    Circular,
    Deflection,
    Eccentric,
    Inclined,
    NearCritical,
    PhotonSphere,
    standard_cases,
)

TIGHT = {"rtol": 1e-13, "atol": 1e-13}

# Well-conditioned cases: the error against the exact reference must be <= 1e-10.
WELL_CONDITIONED = [
    *(Circular(r) for r in (6.0, 7.0, 10.0, 20.0, 100.0)),
    *(Eccentric(p, e, 2.0) for p, e in ((100.0, 0.5), (20.0, 0.5), (7.5, 0.5))),
    *(Inclined(20.0, 0.5, i, 2.0) for i in (30.0, 60.0, 85.0)),
    *(Deflection(b) for b in (5.3, 6.0, 8.0, 10.0, 20.0, 50.0, 100.0, 1000.0)),
]


@pytest.mark.parametrize("case", WELL_CONDITIONED, ids=lambda c: c.id)
def test_hamiltonian_formulation_meets_the_reference(case, hamiltonian):
    sol = get("DOP853").solve(case.problem(hamiltonian), **TIGHT)
    assert case.errors(sol, hamiltonian)["error"] <= 1e-10


@pytest.mark.parametrize(
    "case", [c for c in WELL_CONDITIONED if c.family != "TC1"], ids=lambda c: c.id
)
def test_second_order_formulation_meets_the_reference(case, second_order):
    sol = get("DOP853").solve(case.problem(second_order), **TIGHT)
    assert case.errors(sol, second_order)["error"] <= 1e-10


def test_second_order_loses_accuracy_on_far_rays(second_order, hamiltonian):
    # In (a), L = r^2 u^phi is rebuilt from a component that falls off like 1/r^2, so its
    # relative accuracy is not controlled far from the hole; in (b), L = p_phi is exact.
    # The impact parameter drifts and Delta phi inherits |d Delta phi / d b| times that drift.
    case = Deflection(5.3, 1e4)
    err = {}
    for form in (second_order, hamiltonian):
        sol = get("DOP853").solve(case.problem(form), **TIGHT)
        err[form.name] = case.errors(sol, form)["error"]
        if form is second_order:
            inv = form.invariants(sol.y)
            drift_L = abs(inv["Lz"][-1] - inv["Lz"][0]) / inv["Lz"][0]
    assert err["hamiltonian"] < 1e-11
    assert err["second_order"] > 100 * err["hamiltonian"]
    assert drift_L > 1e-11


@pytest.mark.parametrize("k", [2, 4, 6])
def test_photon_sphere_exit_matches_the_elliptic_integral(k, hamiltonian):
    case = PhotonSphere(k)
    sol = get("DOP853").solve(case.problem(hamiltonian), **TIGHT)
    e = case.errors(sol, hamiltonian)
    # The error is amplified by e^(2 pi) per orbit spent near r = 3.
    assert e["error"] <= 1e-12 * math.exp(2 * math.pi * e["orbits"])
    if k >= 4:
        assert e["growth_rate"] == pytest.approx(1 / (3 * math.sqrt(3)), rel=2e-3)
    # Linear theory misses the exit by a few thousandths of an orbit (nonlinear terms).
    assert abs(case.reference()["orbits_linear"] - e["orbits"]) < 0.01


@pytest.mark.parametrize("k", [2, 4, 6])
def test_near_critical_windings(k, hamiltonian):
    case = NearCritical(k)
    sol = get("DOP853").solve(case.problem(hamiltonian), **TIGHT)
    e = case.errors(sol, hamiltonian)
    assert e["windings"] == pytest.approx(case.reference()["windings"], abs=1e-6)


def test_every_standard_case_runs(hamiltonian):
    for case in standard_cases():
        sol = get("DOP853").solve(case.problem(hamiltonian), rtol=1e-10, atol=1e-10)
        assert sol.success, case.id
        assert np.isfinite(case.errors(sol, hamiltonian)["error"]) or case.family == "TC2"


def test_case_ids_are_unique_and_stable():
    ids = [c.id for c in standard_cases()]
    assert len(ids) == len(set(ids))
    assert Eccentric(20.0, 0.5, 2.0).id == "TC4-p20-e0.5-n_orbits2"
