"""F11: time until an inclined ISCO orbit leaves 5 < r < 7, against tolerance and step size.

The ISCO is marginally stable: the effective potential has an inflection point there. With
r = 6 + delta and L^2 = 12, V(r) = f (1 + L^2/r^2) = 8/9 + delta^3/1944 + O(delta^4), so the
radial equation is delta'' = -delta^2/1296 in proper time. An integrator adds a small force
eps (proportional to h^p, or to the tolerance) and the start is no longer at rest:

    delta'' = eps - delta^2 / 1296.

For eps > 0 there is a stable equilibrium at delta = 36 sqrt(eps) and the orbit oscillates
about it without plunging. For eps < 0 the orbit drifts inward under constant acceleration
until delta ~ -36 sqrt|eps|, then runs away in finite time: tau ~ |eps|^(-1/4). The growth is
algebraic, not exponential as on the photon sphere, and the sign of the error decides between
survival and plunge. With eps ~ h^p the time to plunge scales as (steps per orbit)^(p/4).

Equatorial ISCO orbits are exact fixed points of every Runge-Kutta map, so the orbit is
inclined by 45 degrees (TC3m). Runs are capped at 2000 orbits. Tao4 is run with the global
omega = 1e-3 and with omega = 0.1: at 1e-3 the two copies drift apart along the marginal
direction (by O(1) in r), and the method fails here.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from _common import save_table
from matplotlib.lines import Line2D

from geoint.experiments import run_many
from geoint.experiments.analysis import fit_slope
from geoint.experiments.runner import RunSpec
from geoint.experiments.sweeps import TAO_OMEGA, tolerance_grid
from geoint.plotting import style
from geoint.testcases.schwarzschild import MarginalCircular

CASE = MarginalCircular(45.0, 2000.0)
TOLS = tolerance_grid(1e-13, 1e-3, 11)
STEPS_PER_ORBIT = [int(n) for n in np.geomspace(10, 1000, 7)]
FIXED = ("RK4", "GL2", "GL3", "Tao4")
ORDER = {"RK4": 4, "GL2": 4, "GL3": 6, "Tao4": 4}
TAO_STRONG = 0.1  # omega that holds Tao's copies together at the ISCO


def compute() -> pd.DataFrame:
    specs = [
        RunSpec.make(CASE, form, m, rtol=t, atol=t, save_every=10)
        for m in ("DP5", "DOP853")
        for form in ("a", "b")
        for t in TOLS
    ]
    for m in FIXED:
        for form in ("b",) if m.startswith("Tao") else ("a", "b"):
            for omega in (TAO_OMEGA, TAO_STRONG) if m.startswith("Tao") else (None,):
                for spo in STEPS_PER_ORBIT:
                    opts = {"n_steps": int(CASE.n_orbits * spo), "save_every": spo}
                    if omega is not None:
                        opts["omega"] = omega
                    specs.append(RunSpec.make(CASE, form, m, **opts))
    df = run_many(specs)
    df["steps_per_orbit"] = df["opt_n_steps"] / CASE.n_orbits
    df["plunged"] = df["exit"] == "plunged"
    return df


def series(df):
    """(label, method, formulation, x column, rows) for every curve."""
    out = []
    for m in ("DP5", "DOP853"):
        for form in ("a", "b"):
            g = df[(df["method"] == m) & (df["formulation"] == form)].sort_values("opt_rtol")
            out.append((f"{m} ({form})", m, form, "opt_rtol", g))
    for m in FIXED:
        for form in ("b",) if m.startswith("Tao") else ("a", "b"):
            g = df[(df["method"] == m) & (df["formulation"] == form)]
            if m.startswith("Tao"):
                for omega in (TAO_OMEGA, TAO_STRONG):
                    h = g[g["opt_omega"] == omega].sort_values("steps_per_orbit")
                    out.append((f"{m} ({form}), ω = {omega:g}", m, form, "steps_per_orbit", h))
            else:
                g = g.sort_values("steps_per_orbit")
                out.append((f"{m} ({form})", m, form, "steps_per_orbit", g))
    return out


def scaling_table(df) -> pd.DataFrame:
    """Slope of log(orbits to plunge) against log(steps per orbit) or log(tol), plunges only."""
    rows = []
    for label, m, _, x, g in series(df):
        p = g[g["plunged"]]
        slope, se, _ = fit_slope(p[x], p["orbits"], lo=0, hi=np.inf)
        predicted = -0.25 if x == "opt_rtol" else ORDER[m] / 4
        rows.append(
            {
                "series": label,
                "plunged": int(g["plunged"].sum()),
                "runs": len(g),
                "slope": slope,
                "slope_se": se,
                "predicted": predicted,
            }
        )
    return pd.DataFrame(rows)


def plot(df: pd.DataFrame, table: pd.DataFrame) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.4), sharey=True, constrained_layout=True)
    slopes = table.set_index("series")
    handles = ([], [])
    for label, m, form, x, g in series(df):
        col = 0 if x == "opt_rtol" else 1
        ax = axes[col]
        st = style.method_style(m, form)
        if "ω = 0.1" in label:
            st["linestyle"] = ":"
        ax.plot(g[x], g["orbits"], **{**st, "marker": None})
        for plunged, fill in ((True, "full"), (False, "none")):
            sel = g[g["plunged"] == plunged]
            ax.plot(sel[x], sel["orbits"], **{**st, "linestyle": "none", "fillstyle": fill})
        row = slopes.loc[label]
        if row["plunged"] >= 3:
            label = f"{label}: slope {row['slope']:.2f} ({row['predicted']:g})"
        handles[col].append(Line2D([], [], **st, label=label))
    for ax in axes:
        ax.axhline(CASE.n_orbits, color=style.MUTED, linewidth=0.8)
    axes[0].set(
        xscale="log",
        yscale="log",
        xlabel="rtol = atol",
        ylabel="orbits before |r - 6| > 1",
        title="adaptive pairs",
    )
    axes[1].set(xscale="log", xlabel="steps per ISCO orbit n", title="fixed-step methods")
    axes[0].legend(handles=handles[0], loc="lower left")
    axes[1].legend(handles=handles[1], loc="lower right", fontsize=6.5)
    axes[1].annotate("cap: 2000 orbits", (12, 2300), color=style.INK_2, fontsize=7)
    fig.suptitle(
        "Inclined ISCO orbit. Filled: plunged (r < 5); hollow: left outward (r > 7) or still "
        "bound at the cap.\nSlopes: fit over the plunges (prediction: p/4 against n, "
        "-1/4 against tol)",
        fontsize=8,
        color=style.INK_2,
    )
    style.save(fig, "F11_isco_plunge")


def main():
    df = compute()
    save_table(
        df,
        "f11_isco",
        [
            "formulation",
            "method",
            "opt_rtol",
            "steps_per_orbit",
            "opt_omega",
            "orbits",
            "exit",
            "nfev",
            "dH_max",
        ],
    )
    table = scaling_table(df)
    save_table(table, "f11_isco_scaling")
    plot(df, table)
    return df, table


if __name__ == "__main__":
    _, table = main()
    print(table.to_string())
