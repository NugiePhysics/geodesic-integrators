"""F5: tolerance proportionality of the adaptive pairs, ours and SciPy's.

Global error against rtol (= atol) for DP5/RK45 and DOP853, on the orbit (20, 0.5) over three
radial periods and on the deflection with b = 6, formulation (b). A well-behaved pair has a
global error roughly proportional to the tolerance (slope 1), down to the round-off floor.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import pandas as pd
from _common import save_table

from geoint.experiments import run_many
from geoint.experiments.analysis import fit_slope
from geoint.experiments.sweeps import specs, tolerance_grid
from geoint.plotting import style
from geoint.testcases.schwarzschild import Deflection, Eccentric

CASES = {
    "orbit (20, 0.5), 3 periods": Eccentric(20.0, 0.5, 3.0),
    "deflection b = 6": Deflection(6.0),
}
METHODS = ("DP5", "scipy-RK45", "DOP853", "scipy-DOP853")
TOLERANCES = tolerance_grid(1e-13, 1e-3, 11)


def compute() -> pd.DataFrame:
    frames = []
    for label, case in CASES.items():
        df = run_many(specs(case, METHODS, (), TOLERANCES, formulations=("b",)))
        df["case_label"] = label
        frames.append(df)
    df = pd.concat(frames, ignore_index=True)
    df["error"] = pd.to_numeric(df["error"], errors="coerce")
    return df


def plot(df: pd.DataFrame) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.4), sharey=True, constrained_layout=True)
    for ax, label in zip(axes, CASES, strict=True):
        sub = df[df["case_label"] == label]
        for method in METHODS:
            g = sub[sub["method"] == method].sort_values("opt_rtol")
            slope, _, _ = fit_slope(g["opt_rtol"], g["error"], lo=1e-12, hi=1e-2)
            ax.plot(
                g["opt_rtol"],
                g["error"],
                **style.method_style(method, "b", linestyle="-" if "scipy" not in method else ":"),
                label=f"{method} (slope {slope:.2f})",
            )
        ax.plot([1e-13, 1e-3], [1e-13, 1e-3], color=style.MUTED, linewidth=0.8, linestyle="--")
        ax.annotate(
            "error = tol",
            (1e-6, 1e-6),
            color=style.INK_2,
            fontsize=7,
            xytext=(4, -10),
            textcoords="offset points",
        )
        ax.set(xscale="log", yscale="log", xlabel="rtol = atol", title=label)
        ax.legend(loc="upper left")
    axes[0].set_ylabel("error (rad)")
    style.save(fig, "F05_tolerance")


def main():
    df = compute()
    save_table(
        df,
        "f05_tolerance",
        ["case_label", "method", "opt_rtol", "error", "nfev", "n_steps", "status"],
    )
    plot(df)
    return df


if __name__ == "__main__":
    main()
