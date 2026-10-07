"""F6 and T5: work-precision diagrams and the cost of reaching a target error.

For three test cases, every method is run over a range of step sizes (fixed step) or
tolerances (adaptive). F6 plots the error against the number of vector-field evaluations and
against wall time (median of repeats after compilation, our compiled implementations only;
ADR 0003), one figure per formulation. T5 interpolates each method's work-precision front at
the target errors 1e-6 and 1e-10.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import pandas as pd
from _common import SUMMARY, markdown_table, save_table

from geoint.experiments import run_many
from geoint.experiments.analysis import cost_at_error, pareto_front
from geoint.experiments.sweeps import ALL, specs, step_grid, tolerance_grid
from geoint.plotting import style
from geoint.testcases.schwarzschild import Deflection, Eccentric, Inclined

CASES = {
    "deflection b = 6": (Deflection(6.0), step_grid(300, 1_000_000, 12)),
    "orbit (20, 0.5), 3 periods": (Eccentric(20.0, 0.5, 3.0), step_grid(150, 300_000, 12)),
    "inclined orbit, 60°, 2 periods": (Inclined(20.0, 0.5, 60.0, 2.0), step_grid(100, 200_000, 12)),
}
TOLERANCES = tolerance_grid(1e-13, 1e-3, 11)
TARGETS = (1e-6, 1e-10)


def compute() -> pd.DataFrame:
    frames = []
    for label, (case, steps) in CASES.items():
        df = run_many(specs(case, ALL, steps, TOLERANCES), repeats=3)
        df["case_label"] = label
        frames.append(df)
    df = pd.concat(frames, ignore_index=True)
    df["error"] = pd.to_numeric(df["error"], errors="coerce")
    return df


def cost_table(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (label, form, method), g in df.groupby(["case_label", "formulation", "method"]):
        row = {"case": label, "formulation": form, "method": method}
        for target in TARGETS:
            row[f"nfev@{target:g}"] = cost_at_error(g["nfev"], g["error"], target)
            row[f"time_ms@{target:g}"] = 1e3 * cost_at_error(g["wall_time"], g["error"], target)
        rows.append(row)
    return pd.DataFrame(rows)


def plot(df: pd.DataFrame, form: str) -> None:
    fig, axes = plt.subplots(len(CASES), 2, figsize=(8.4, 9.0), constrained_layout=True)
    for row, label in enumerate(CASES):
        sub = df[(df["case_label"] == label) & (df["formulation"] == form)]
        for col, cost in enumerate(("nfev", "wall_time")):
            ax = axes[row, col]
            for method in ALL:
                g = sub[sub["method"] == method]
                if g.empty:
                    continue
                c, e = pareto_front(g[cost], g["error"])
                if cost == "wall_time":
                    c = c * 1e3
                ax.plot(c, e, **style.method_style(method, form, markersize=3.5), label=method)
            ax.set(xscale="log", yscale="log", ylim=(1e-15, 1))
            if row == len(CASES) - 1:
                ax.set_xlabel("vector-field evaluations" if cost == "nfev" else "wall time (ms)")
            if col == 0:
                ax.set_ylabel(f"{label}\nerror (rad)")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside lower center", ncol=8)
    fig.suptitle(f"Work-precision, formulation {style.FORMULATIONS[form]['label']}", fontsize=10)
    style.save(fig, f"F06_work_precision_{form}")


def main():
    df = compute()
    save_table(
        df,
        "f06_work_precision",
        [
            "case_label",
            "formulation",
            "method",
            "opt_n_steps",
            "opt_rtol",
            "error",
            "nfev",
            "n_steps",
            "n_rejected",
            "n_iter",
            "wall_time",
            "status",
        ],
    )
    table = cost_table(df)
    save_table(table, "t5_cost_to_target")
    (SUMMARY / "t5_cost_to_target.md").write_text(markdown_table(table) + "\n")
    for form in ("a", "b"):
        plot(df, form)
    return df, table


if __name__ == "__main__":
    _, table = main()
    print(table.to_string())
