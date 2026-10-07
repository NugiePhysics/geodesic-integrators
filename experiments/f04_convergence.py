"""F4 and T3: convergence of the fixed-step methods in both formulations.

Error against the exact reference as a function of the step size h, for the strong-field orbit
(p, e) = (7.5, 0.5) over one radial period (phase error at the final periapsis) and for the
deflection of a photon with b = 6 (error in Delta phi at r_far = 1000). T3 lists the slopes,
fitted over the asymptotic range between the round-off floor and the pre-asymptotic regime.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import pandas as pd
from _common import SUMMARY, markdown_table, save_table

from geoint.experiments import run_many
from geoint.experiments.analysis import fit_slope
from geoint.experiments.sweeps import FIXED, specs, step_grid
from geoint.integrators import get
from geoint.plotting import style
from geoint.testcases.schwarzschild import Deflection, Eccentric

CASES = {
    "orbit": (Eccentric(7.5, 0.5, 1.0), step_grid(40, 100_000, 14)),
    "deflection": (Deflection(6.0), step_grid(200, 1_600_000, 16)),
}


def compute() -> pd.DataFrame:
    frames = []
    for label, (case, steps) in CASES.items():
        df = run_many(specs(case, FIXED, steps, ()))
        df["case_label"] = label
        frames.append(df)
    df = pd.concat(frames, ignore_index=True)
    df["error"] = pd.to_numeric(df["error"], errors="coerce")
    df["iter_per_step"] = df["n_iter"] / df["n_steps"]
    return df


def orders(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (label, form, method), g in df.groupby(["case_label", "formulation", "method"]):
        # Errors between 1e-12 and 1e-5: above the round-off floor, below the regime where
        # the photon's near-critical passage makes the error pre-asymptotic.
        slope, se, n = fit_slope(g["h"], g["error"], lo=1e-12, hi=1e-5)
        rows.append(
            {
                "case": label,
                "formulation": form,
                "method": method,
                "theory": get(method).order,
                "measured": slope,
                "fit_se": se,
                "points": n,
            }
        )
    return pd.DataFrame(rows)


def plot(df: pd.DataFrame) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(8.6, 3.5), sharey=True, constrained_layout=True)
    sub = df[df["case_label"] == "orbit"]
    for ax, form in zip(axes, ("a", "b"), strict=True):
        for method in FIXED:
            g = sub[(sub["formulation"] == form) & (sub["method"] == method)].sort_values("h")
            if g.empty:
                continue
            ax.plot(g["h"], g["error"], **style.method_style(method, form), label=method)
        # Each guide runs one decade parallel to (and a decade below) a method of that order.
        for p, method in ((2, "GL1"), (4, "GL2"), (6, "GL3")):
            g = sub[(sub["formulation"] == form) & (sub["method"] == method)].sort_values("h")
            g = g[(g["error"] > 1e-11) & (g["error"] < 1e-3)]
            h0 = g["h"].iloc[0]
            style.reference_slope(
                ax, h0, 0.1 * g["error"].iloc[0], p, length=1.0, label=rf"$h^{p}$"
            )
        ax.set(
            xscale="log",
            yscale="log",
            xlabel="step size h",
            title=f"{style.FORMULATIONS[form]['label']}",
            ylim=(1e-15, 10),
        )
    axes[0].set_ylabel("phase error after one radial period (rad)")
    axes[1].legend(loc="lower right", ncol=2)
    fig.suptitle("Convergence on the orbit (p, e) = (7.5, 0.5)", fontsize=10)
    style.save(fig, "F04_convergence")


def main():
    df = compute()
    save_table(
        df,
        "f04_convergence",
        [
            "case_label",
            "formulation",
            "method",
            "opt_n_steps",
            "h",
            "error",
            "nfev",
            "n_iter",
            "iter_per_step",
            "status",
        ],
    )
    table = orders(df)
    save_table(table, "t3_orders")
    (SUMMARY / "t3_orders.md").write_text(markdown_table(table) + "\n")
    plot(df)
    return df, table


if __name__ == "__main__":
    _, table = main()
    print(table.to_string())
