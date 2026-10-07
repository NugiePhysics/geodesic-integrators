"""F3: deflection error against impact parameter for every method at a fixed cost.

For each b (r_far = 1000) every method runs a small sweep of step sizes or tolerances; the
error at a budget of 2x10^4 vector-field evaluations is read off its work-precision front.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import pandas as pd
from _common import save_table

from geoint.experiments import run_many
from geoint.experiments.analysis import error_at_cost
from geoint.experiments.sweeps import ALL, specs, step_grid, tolerance_grid
from geoint.plotting import style
from geoint.testcases.schwarzschild import TC1_IMPACT_PARAMETERS, Deflection

BUDGET = 2e4


def compute() -> pd.DataFrame:
    specs_ = []
    for b in TC1_IMPACT_PARAMETERS:
        specs_ += specs(
            Deflection(b), ALL, step_grid(500, 200_000, 7), tolerance_grid(1e-12, 1e-4, 7)
        )
    df = run_many(specs_)
    df["error"] = pd.to_numeric(df["error"], errors="coerce")
    rows = []
    for (b, form, method), g in df.groupby(["case_b", "formulation", "method"]):
        rows.append(
            {
                "b": b,
                "formulation": form,
                "method": method,
                "error_at_budget": error_at_cost(g["nfev"], g["error"], BUDGET),
            }
        )
    return pd.DataFrame(rows)


def plot(df: pd.DataFrame) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(8.6, 3.6), sharey=True, constrained_layout=True)
    for ax, form in zip(axes, ("a", "b"), strict=True):
        for method in ALL:
            g = df[(df["formulation"] == form) & (df["method"] == method)].sort_values("b")
            if g.empty:
                continue
            ax.plot(g["b"], g["error_at_budget"], **style.method_style(method, form), label=method)
        ax.set(
            xscale="log",
            yscale="log",
            xlabel="impact parameter b (M)",
            title=style.FORMULATIONS[form]["label"],
        )
    axes[0].set_ylabel(rf"error in $\Delta\phi$ at {BUDGET:.0e} evaluations".replace("e+04", "e4"))
    handles, labels = axes[1].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside lower center", ncol=8)
    fig.suptitle(r"Light deflection, $r_{far} = 1000$, fixed cost", fontsize=10)
    style.save(fig, "F03_deflection_cost")


def main():
    df = compute()
    save_table(df, "f03_deflection_cost")
    plot(df)
    return df


if __name__ == "__main__":
    main()
