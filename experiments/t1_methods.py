"""T1: properties of every integrator of the study.

Order, structure and nominal cost come from the integrators themselves. The measured cost of
the implicit methods (fixed-point iterations to round-off stagnation) is read from T4: the
evaluations per step over 10^4 periods of (20, 0.5) at 1000 steps per period, so this script
runs after f07.
"""

from __future__ import annotations

import pandas as pd
from _common import SUMMARY, markdown_table, save_table

from geoint.experiments.sweeps import ALL, TAO_OMEGA
from geoint.integrators import get

STEPS_PER_PERIOD = 1000  # of the T4 runs
PARAMETERS = {
    "RK4": "h",
    "DP5": "rtol, atol",
    "DOP853": "rtol, atol",
    "GL1": "h; iteration to stagnation",
    "GL2": "h; iteration to stagnation",
    "GL3": "h; iteration to stagnation",
    "Tao2": f"h; ω = {TAO_OMEGA:g}, (r, θ) coupled",
    "Tao4": f"h; ω = {TAO_OMEGA:g}, (r, θ) coupled",
}
STRUCTURE = {
    "RK4": "explicit Runge-Kutta",
    "DP5": "explicit embedded pair 5(4)",
    "DOP853": "explicit embedded pair 8(5,3)",
    "GL1": "implicit midpoint (Gauss, s = 1)",
    "GL2": "Gauss collocation, s = 2",
    "GL3": "Gauss collocation, s = 3",
    "Tao2": "extended phase space, Strang",
    "Tao4": "extended phase space, triple jump",
}


def table() -> pd.DataFrame:
    t4 = pd.read_csv(SUMMARY / "t4_growth.csv").set_index(["method", "formulation"])
    rows = []
    for name in ALL:
        method = get(name)
        p = method.properties()
        row = {
            "method": name,
            "structure": STRUCTURE[name],
            "order": p["order"],
            "explicit": "no" if p["implicit"] else "yes",
            "symplectic": "yes" if p["symplectic"] else "no",
            "symmetric": "yes" if p["symmetric"] else "no",
            "adaptive": "yes" if p["adaptive"] else "no",
            "evals_per_step": str(p["evals_per_step"]),
        }
        if p["implicit"]:
            a, b = (t4.loc[(name, f), "nfev_per_orbit"] / STEPS_PER_PERIOD for f in "ab")
            row["evals_per_step"] += f" per iteration; measured {b:.1f} (b), {a:.1f} (a)"
        row["parameters"] = PARAMETERS[name]
        row["formulations"] = "(b) only" if method.requires_hamiltonian else "(a), (b)"
        rows.append(row)
    return pd.DataFrame(rows)


def main():
    t = table()
    save_table(t, "t1_methods")
    (SUMMARY / "t1_methods.md").write_text(markdown_table(t) + "\n")
    return t


if __name__ == "__main__":
    print(main().to_string())
