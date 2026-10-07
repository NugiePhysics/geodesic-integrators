"""F20 and T5b: work-precision over 1000 radial periods of the orbit (20, 0.5) (RQ3).

F6 and T5 answer RQ3 for short integrations (a few periods), where DOP853 wins everywhere.
This repeats the sweep over 1000 periods, so that errors which grow with time (linear drift
in H, quadratic phase error) have time to matter. Two error measures, against the number of
vector-field evaluations: the phase error |phi - N Phi| at lambda = N T_tau, and the largest
relative mass-shell error |delta H| along the run. T5b interpolates the cost of reaching
phase errors 1e-4 and 1e-8 and a mass-shell error 1e-10, per period, for comparison with T5.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import pandas as pd
from _common import SUMMARY, markdown_table, save_table

from geoint.experiments import RunSpec, run_many
from geoint.experiments.analysis import cost_at_error, pareto_front
from geoint.experiments.sweeps import ALL, TAO_OMEGA, formulations_for, step_grid, tolerance_grid
from geoint.integrators import get
from geoint.plotting import style
from geoint.testcases.schwarzschild import Eccentric

N_ORBITS = 1000
CASE = Eccentric(20.0, 0.5, float(N_ORBITS))
STEPS_PER_ORBIT = step_grid(30, 4000, 9)
TOLERANCES = tolerance_grid(1e-13, 1e-4, 10)
PHASE_TARGETS = (1e-4, 1e-8)
DH_TARGET = 1e-10
SAVED_POINTS = 20_000  # per run; dH_max is taken over the saved points


def all_specs() -> list[RunSpec]:
    out = []
    for method in ALL:
        for form in formulations_for(method):
            if get(method).is_adaptive:
                for tol in TOLERANCES:
                    out.append(RunSpec.make(CASE, form, method, rtol=tol, atol=tol))
                continue
            for spo in STEPS_PER_ORBIT:
                n = N_ORBITS * spo
                opts = {"n_steps": n, "save_every": max(1, n // SAVED_POINTS)}
                if get(method).requires_hamiltonian:
                    opts["omega"] = TAO_OMEGA
                out.append(RunSpec.make(CASE, form, method, **opts))
    return out


def compute() -> pd.DataFrame:
    df = run_many(all_specs())
    for col in ("error", "dH_max"):
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df["nfev_per_orbit"] = df["nfev"] / N_ORBITS
    return df


def cost_table(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (form, method), g in df.groupby(["formulation", "method"]):
        row = {"formulation": form, "method": method}
        for target in PHASE_TARGETS:
            row[f"phase {target:g}"] = cost_at_error(g["nfev_per_orbit"], g["error"], target)
        row[f"dH {DH_TARGET:g}"] = cost_at_error(g["nfev_per_orbit"], g["dH_max"], DH_TARGET)
        rows.append(row)
    return pd.DataFrame(rows).sort_values(["formulation", "method"], ignore_index=True)


def plot(df: pd.DataFrame) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(8.4, 6.6), sharex=True, constrained_layout=True)
    for col, form in enumerate(("a", "b")):
        sub = df[df["formulation"] == form]
        for row, (key, label) in enumerate(
            (("error", r"phase error $|\phi - N\Phi|$ (rad)"), ("dH_max", r"max $|\delta H|$"))
        ):
            ax = axes[row, col]
            for method in ALL:
                g = sub[sub["method"] == method]
                if g.empty:
                    continue
                c, e = pareto_front(g["nfev_per_orbit"], g[key])
                ax.plot(c, e, **style.method_style(method, form, markersize=3.5), label=method)
            ax.set(xscale="log", yscale="log", ylim=(1e-13, 1e3) if row == 0 else (1e-16, 1))
            if row == 0:
                ax.set_title(style.FORMULATIONS[form]["label"])
            if row == 1:
                ax.set_xlabel("vector-field evaluations per radial period")
            if col == 0:
                ax.set_ylabel(label)
    handles, labels = axes[0, 1].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside lower center", ncol=8)
    fig.suptitle("Work-precision over 1000 radial periods of (20, 0.5)", fontsize=10)
    style.save(fig, "F20_long_work_precision")


def main():
    df = compute()
    save_table(
        df,
        "f20_long_work_precision",
        [
            "formulation",
            "method",
            "opt_n_steps",
            "opt_rtol",
            "error",
            "dH_max",
            "nfev",
            "n_steps",
            "n_rejected",
            "n_iter",
            "status",
        ],
    )
    table = cost_table(df)
    save_table(table, "t5b_cost_long")
    (SUMMARY / "t5b_cost_long.md").write_text(markdown_table(table) + "\n")
    plot(df)
    return df, table


if __name__ == "__main__":
    _, table = main()
    with pd.option_context("display.width", 200):
        print(table.to_string())
