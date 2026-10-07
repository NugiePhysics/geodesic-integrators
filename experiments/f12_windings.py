"""F12: windings of near-critical photons, b = b_c (1 + 10^-k), against -ln(b/b_c - 1).

Numerical windings Delta phi(r_far) / 2 pi (DOP853, tol 1e-13, formulation (b), r_far = 1000)
against the exact elliptic integral and Bozza's logarithmic prediction, whose slope is 1/2pi
windings per e-fold of b - b_c. Bottom: error of the integration, which grows like the
condition number d(Delta phi)/db ~ 1/(b - b_c) (pitfall 19).
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
from geoint.testcases.schwarzschild import NearCritical

K_GRID = np.linspace(0.5, 10, 20)


def compute() -> pd.DataFrame:
    rows = []
    for form in ("a", "b"):
        specs = [
            RunSpec.make(NearCritical(float(k)), form, "DOP853", rtol=1e-13, atol=1e-13)
            for k in K_GRID
        ]
        rows.append(run_many(specs))
    df = pd.concat(rows, ignore_index=True)
    df["b"] = analytic.B_CRIT * (1 + 10.0 ** -df["case_k"])
    df["x"] = -np.log(df["b"] / analytic.B_CRIT - 1)
    df["windings_exact"] = [analytic.windings(b, 1000.0) for b in df["b"]]
    df["windings_bozza"] = [
        (analytic.bozza_deflection(b) + math.pi) / (2 * math.pi) for b in df["b"]
    ]
    df["error"] = pd.to_numeric(df["error"], errors="coerce")
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
    b = df[df["formulation"] == "b"].sort_values("x")
    top.plot(
        b["x"],
        b["windings_exact"],
        color=style.INK,
        linewidth=1.6,
        label=r"exact, $r_{far} = 1000$",
    )
    top.plot(
        b["x"],
        b["windings_bozza"],
        color="#eb6834",
        linestyle=":",
        linewidth=1.4,
        label=r"Bozza, slope $1/2\pi$ (asymptotic)",
    )
    for form in ("a", "b"):
        g = df[df["formulation"] == form].sort_values("x")
        st = style.method_style("DOP853", form, markersize=3.5)
        top.plot(
            g["x"],
            g["windings"],
            linestyle="none",
            marker=st["marker"],
            color=st["color"],
            fillstyle="none" if form == "a" else "full",
            label=f"DOP853, {style.FORMULATIONS[form]['label']}",
        )
        bottom.plot(
            g["x"], g["error"] / (2 * math.pi), **st, label=style.FORMULATIONS[form]["label"]
        )
    kappa = 1 / (b["b"] - analytic.B_CRIT)
    bottom.plot(
        b["x"],
        1e-13 * kappa,
        color=style.MUTED,
        linestyle="--",
        linewidth=0.9,
        label=r"$10^{-13} / (b - b_c)$",
    )
    top.set(ylabel="windings  $\\Delta\\phi / 2\\pi$")
    top.legend(loc="upper left")
    bottom.set(yscale="log", xlabel=r"$-\ln(b/b_c - 1)$", ylabel="error (windings)")
    bottom.legend(loc="upper left")
    style.save(fig, "F12_windings")


def main():
    df = compute()
    save_table(
        df,
        "f12_windings",
        [
            "formulation",
            "case_k",
            "b",
            "x",
            "windings",
            "windings_exact",
            "windings_bozza",
            "error",
        ],
    )
    plot(df)
    return df


if __name__ == "__main__":
    main()
