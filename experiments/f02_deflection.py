"""F2: deflection angle against impact parameter, numerical and every approximation.

Top: the exact asymptotic angle (Darwin), the weak-field 4M/b and fourth-order series, Bozza's
strong-deflection limit, and the numerical Delta phi(r_far) - pi with r_far = 10^4 (DOP853,
tol 1e-13, formulation (b)), against b - b_c. Bottom: relative error of each approximation
with respect to the exact angle, and of the integration with respect to the exact
Delta phi(r_far) at the same r_far (pitfall 18).
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
from geoint.testcases.schwarzschild import Deflection

B_GRID = analytic.B_CRIT + np.geomspace(1e-3, 2e3, 26)  # b < r_far: a photon must reach r_far
R_FAR = 1e4


def compute() -> pd.DataFrame:
    cases = [Deflection(float(b), R_FAR) for b in B_GRID]
    df = run_many([RunSpec.make(c, "b", "DOP853", rtol=1e-13, atol=1e-13) for c in cases])
    df["b"] = df["case_b"]
    df["alpha_exact"] = [analytic.deflection_angle(b) for b in df["b"]]
    df["alpha_numeric"] = df["delta_phi"] - math.pi
    df["alpha_exact_rfar"] = [analytic.deflection_angle(b, R_FAR) for b in df["b"]]
    df["numeric_rel_error"] = pd.to_numeric(df["error"]) / df["delta_phi"]
    df["weak1"] = 4 / df["b"]
    df["weak4"] = [analytic.weak_field_deflection(b, 4) for b in df["b"]]
    df["bozza"] = [analytic.bozza_deflection(b) for b in df["b"]]
    return df


def plot(df: pd.DataFrame) -> None:
    fig, (top, bottom) = plt.subplots(
        2,
        1,
        figsize=(6.4, 5.6),
        sharex=True,
        gridspec_kw={"height_ratios": [3, 2]},
        constrained_layout=True,
    )
    x = df["b"] - analytic.B_CRIT
    approx = {
        "weak field 4M/b": ("weak1", "#2a78d6", "--"),
        "weak field, 4th order": ("weak4", "#1baf7a", "-."),
        "Bozza (strong deflection)": ("bozza", "#eb6834", ":"),
    }
    top.plot(
        x, df["alpha_exact"], color=style.INK, linewidth=1.6, label=r"exact, $r_{far} = \infty$"
    )
    top.plot(
        x,
        df["alpha_exact_rfar"],
        color=style.MUTED,
        linewidth=1.0,
        label=r"exact, $r_{far} = 10^4$",
    )
    for label, (col, color, ls) in approx.items():
        top.plot(x, df[col], color=color, linestyle=ls, linewidth=1.3, label=label)
        bottom.plot(
            x,
            np.abs(df[col] - df["alpha_exact"]) / df["alpha_exact"],
            color=color,
            linestyle=ls,
            linewidth=1.3,
        )
    top.plot(
        x,
        df["alpha_numeric"],
        linestyle="none",
        marker="o",
        markersize=3.5,
        color=style.INK,
        fillstyle="none",
        label=r"DOP853, $\Delta\phi(r_{far}) - \pi$",
    )
    bottom.plot(
        x,
        df["numeric_rel_error"],
        linestyle="none",
        marker="o",
        markersize=3.5,
        color=style.INK,
        fillstyle="none",
    )
    top.set(xscale="log", yscale="log", ylabel=r"deflection angle $\alpha$ (rad)")
    top.legend(loc="lower left")
    bottom.set(
        xscale="log",
        yscale="log",
        xlabel=r"$b - b_c$ (M)",
        ylabel="relative error",
        ylim=(1e-17, 10),
    )
    bottom.annotate(
        r"integration vs exact $\Delta\phi(r_{far})$", (3e-3, 3e-15), color=style.INK_2, fontsize=7
    )
    style.save(fig, "F02_deflection")


def main():
    df = compute()
    save_table(
        df,
        "f02_deflection",
        [
            "b",
            "alpha_exact",
            "alpha_numeric",
            "delta_phi",
            "numeric_rel_error",
            "weak1",
            "weak4",
            "bozza",
        ],
    )
    plot(df)
    return df


if __name__ == "__main__":
    main()
