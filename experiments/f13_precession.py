"""F13: periapsis advance per radial period against p (e = 0.5), near and far from the separatrix.

Numerical advance from the first periapsis event (DOP853, tol 1e-13, both formulations) against
the exact elliptic formula Phi - 2 pi and the weak-field estimate 6 pi M / p. Bottom: error of
the integration, and relative error of the weak-field estimate.
"""

from __future__ import annotations

import math

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from _common import save_table

from geoint import analytic
from geoint.experiments import RunSpec, run_many
from geoint.plotting import style
from geoint.testcases.schwarzschild import Eccentric

E = 0.5
P_GRID = 7.0 + np.geomspace(0.02, 200, 16)


def compute() -> pd.DataFrame:
    specs = [
        RunSpec.make(Eccentric(float(p), E, 1.2), form, "DOP853", rtol=1e-13, atol=1e-13)
        for form in ("a", "b")
        for p in P_GRID
    ]
    df = run_many(specs)
    df["p"] = df["case_p"]
    df["exact"] = [analytic.EccentricOrbit(p, E).precession for p in df["p"]]
    df["weak_field"] = 6 * math.pi / df["p"]
    df["abs_error"] = np.abs(df["precession"] - df["exact"])
    return df


def plot(df: pd.DataFrame) -> None:
    fig, (top, bottom) = plt.subplots(
        2,
        1,
        figsize=(6.4, 5.4),
        sharex=True,
        gridspec_kw={"height_ratios": [3, 2]},
        constrained_layout=True,
    )
    df = df.assign(x=df["p"] - analytic.separatrix(E))
    b = df[df["formulation"] == "b"].sort_values("p")
    top.plot(b["x"], b["exact"], color=style.INK, linewidth=1.6, label=r"exact, $\Phi - 2\pi$")
    top.plot(
        b["p"],
        b["weak_field"],
        color="#2a78d6",
        linestyle="--",
        linewidth=1.3,
        label=r"weak field $6\pi M/p$",
    )
    for form in ("a", "b"):
        g = df[df["formulation"] == form].sort_values("p")
        st = style.method_style("DOP853", form, markersize=3.5)
        top.plot(
            g["p"],
            g["precession"],
            linestyle="none",
            marker=st["marker"],
            color=st["color"],
            fillstyle="none" if form == "a" else "full",
            label=f"DOP853, {style.FORMULATIONS[form]['label']}",
        )
        bottom.plot(g["x"], g["abs_error"], **st, label=f"integration, {form}")
    bottom.plot(
        b["p"],
        np.abs(b["weak_field"] - b["exact"]) / b["exact"],
        color="#2a78d6",
        linestyle="--",
        linewidth=1.3,
        label="weak field (relative)",
    )
    top.axvline(analytic.separatrix(E), color=style.MUTED, linewidth=0.8)
    top.annotate(
        "separatrix p = 6 + 2e",
        (analytic.separatrix(E), 0.02),
        xycoords=("data", "axes fraction"),
        xytext=(3, 0),
        textcoords="offset points",
        color=style.INK_2,
        fontsize=7,
    )
    top.set(xscale="log", yscale="log", ylabel="periapsis advance per period (rad)")
    top.legend(loc="upper right")
    bottom.set(
        xscale="log",
        yscale="log",
        ylabel="error",
        xlabel=r"distance from the separatrix, $p - (6 + 2e)$ (M)",
    )
    bottom.legend(loc="center right")
    style.save(fig, "F13_precession")


def main():
    df = compute()
    save_table(
        df, "f13_precession", ["formulation", "p", "precession", "exact", "weak_field", "abs_error"]
    )
    plot(df)
    return df


if __name__ == "__main__":
    main()
