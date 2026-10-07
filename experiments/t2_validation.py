"""T2: numerical against analytical values for every test case (DOP853, tol 1e-13, (b))."""

from __future__ import annotations

import math

import pandas as pd
from _common import SUMMARY, markdown_table, save_table

from geoint import analytic
from geoint.experiments import RunSpec, solve
from geoint.testcases.schwarzschild import (
    Circular,
    Deflection,
    Eccentric,
    Inclined,
    NearCritical,
    PhotonSphere,
)

TIGHT = {"rtol": 1e-13, "atol": 1e-13}


def run(case, form="b"):
    sol, f = solve(RunSpec.make(case, form, "DOP853", **TIGHT))
    return sol, f, case.errors(sol, f)


def compute() -> pd.DataFrame:
    rows = []

    def add(check, numerical, analytical):
        rows.append(
            {
                "check": check,
                "numerical": numerical,
                "analytical": analytical,
                "difference": abs(numerical - analytical),
            }
        )

    for b, r_far in ((6.0, 1e3), (40.0, 1e4), (1000.0, 1e4)):
        case = Deflection(b, r_far)
        _, _, e = run(case)
        add(
            f"TC1 deflection Δφ(r_far), b = {b:g}, r_far = {r_far:g}",
            e["delta_phi"],
            case.reference()["delta_phi"],
        )
    case = NearCritical(6.0)
    _, _, e = run(case)
    add("TC1b windings, b = b_c(1 + 1e-6)", e["windings"], case.reference()["windings"])
    case = PhotonSphere(6.0)
    _, _, e = run(case)
    add("TC2 orbits before |r - 3| > 0.1, r0 = 3 + 1e-6", e["orbits"], case.reference()["orbits"])
    add(
        "TC2 growth rate of |r - 3| in t (Lyapunov exponent)",
        e["growth_rate"],
        analytic.LYAPUNOV_PHOTON_SPHERE,
    )
    for r_c in (10.0, 6.0):
        case = Circular(r_c, 5.0)
        _, _, e = run(case)
        add(
            f"TC3 circular orbit r_c = {r_c:g}: dφ/dt",
            e["omega"],
            analytic.circular_orbit_frequency(r_c),
        )
    for p, e_ in ((20.0, 0.5), (7.5, 0.5)):
        case = Eccentric(p, e_, 1.2)
        sol, _, e = run(case)
        add(
            f"TC4 periapsis advance, (p, e) = ({p:g}, {e_:g})",
            e["precession"],
            case.reference()["precession"],
        )
        add(
            f"TC4 radial period in proper time, (p, e) = ({p:g}, {e_:g})",
            sol.events["periapsis"].lam[0],
            case.reference()["T_tau"],
        )
    case = Inclined(20.0, 0.5, 60.0, 2.0)
    _, _, e = run(case)
    add(
        "TC6 inclined 60°: in-plane angle after 2 periods",
        e["dpsi"] + 2 * case.reference()["Phi"],
        2 * case.reference()["Phi"],
    )
    add("TC6 inclined 60°: tilt of the orbital plane (rad)", e["tilt"], 0.0)
    add("TC6 inclined 60°: max relative change of L²", e["dL2_max"], 0.0)
    add(
        "Old README: precession r_apo = 20, L = 4.2",
        analytic.EccentricOrbit(
            *analytic.EL_to_pe(math.sqrt(0.9 * (1 + 4.2**2 / 400)), 4.2)
        ).precession,
        2.127093985,
    )
    return pd.DataFrame(rows)


def main():
    df = compute()
    save_table(df, "t2_validation")
    (SUMMARY / "t2_validation.md").write_text(markdown_table(df, ".12g") + "\n")
    return df


if __name__ == "__main__":
    print(main().to_string())
