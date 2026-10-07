"""F10 and T6: how long the photon-sphere orbit survives, against tolerance and step size (RQ4).

A photon starts exactly on the photon sphere r = 3 (b = b_c) in a plane inclined by 45 degrees;
in the equatorial plane the orbit is an exact fixed point of every Runge-Kutta map (TC2). The
integration error delta0 is amplified by exp(lambda_L t), lambda_L = 1/(3 sqrt 3), until
|r - 3| > 0.1, after about ln(0.1 / delta0) / (2 pi) orbits. If delta0 is proportional to the
tolerance, every decade of tolerance buys ln(10)/(2 pi) = 0.37 orbits, and double precision
caps the gain near ln(0.1 / 1e-16) / (2 pi) = 5.5 orbits whatever the method.

Panels: (a) log|r - 3| against coordinate time for DOP853 at several tolerances; (b) orbits
survived against the tolerance (adaptive pairs); (c) against the step size (fixed-step
methods); (d) equatorial orbits with an explicit offset r0 - 3 = 10^-k against the exact
elliptic-integral prediction, down to the round-off floor.
"""

from __future__ import annotations

import math

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from _common import SUMMARY, markdown_table, save_table

from geoint import analytic
from geoint.experiments import RunSpec, run_many, solve
from geoint.experiments.sweeps import TAO_OMEGA, tolerance_grid
from geoint.plotting import style
from geoint.testcases.schwarzschild import PhotonSphere

CASE = PhotonSphere(math.inf, 45.0)
TOLS = tolerance_grid(1e-13, 1e-3, 11)
ORBIT = 2 * math.pi * 9 / analytic.B_CRIT  # affine length of one photon-sphere orbit
STEPS_PER_ORBIT = [int(n) for n in np.geomspace(8, 4000, 10)]
FIXED = ("RK4", "GL2", "GL3", "Tao4")
PER_DECADE = math.log(10) / (2 * math.pi)
#: Orbits survived from an offset of one machine epsilon (linear theory).
ROUNDOFF_ORBITS = math.acosh(0.1 / np.finfo(float).eps) / (2 * math.pi)


def compute():
    adaptive = run_many(
        [
            RunSpec.make(CASE, form, m, rtol=t, atol=t)
            for m in ("DP5", "DOP853")
            for form in ("a", "b")
            for t in TOLS
        ]
    )
    fixed_specs = []
    for m in FIXED:
        for form in ("b",) if m.startswith("Tao") else ("a", "b"):
            for spo in STEPS_PER_ORBIT:
                n = int(round(CASE.lam_end / ORBIT * spo))
                opts = {"n_steps": n, "save_every": max(1, spo // 8)}
                if m.startswith("Tao"):
                    opts["omega"] = TAO_OMEGA
                fixed_specs.append(RunSpec.make(CASE, form, m, **opts))
    fixed = run_many(fixed_specs)
    fixed["steps_per_orbit"] = fixed["opt_n_steps"] / (CASE.lam_end / ORBIT)
    offsets = run_many(
        [
            RunSpec.make(PhotonSphere(float(k)), "b", "DOP853", rtol=1e-13, atol=1e-13)
            for k in range(2, 16)
        ]
    )
    offsets["exact"] = [PhotonSphere(float(k)).reference()["orbits"] for k in offsets["case_k"]]
    traces = {}
    for tol in (1e-5, 1e-8, 1e-11, 1e-13):
        sol, _ = solve(RunSpec.make(CASE, "b", "DOP853", rtol=tol, atol=tol))
        traces[tol] = (sol.y[:, 0], np.abs(sol.y[:, 1] - 3.0))
    return adaptive, fixed, offsets, traces


def survival_table(adaptive, fixed) -> pd.DataFrame:
    rows = []
    for (m, form), g in adaptive.groupby(["method", "formulation"]):
        g = g.sort_values("opt_rtol")
        k = np.polyfit(np.log10(g["opt_rtol"]), g["orbits"], 1)
        rows.append(
            {
                "method": m,
                "formulation": form,
                "setting": "tolerance",
                "orbits_tightest": g["orbits"].iloc[0],
                "orbits_loosest": g["orbits"].iloc[-1],
                "orbits_per_decade": -k[0],
                "fit_points": len(g),
                "predicted_per_decade": PER_DECADE,
            }
        )
    for (m, form), g in fixed.groupby(["method", "formulation"]):
        g = g.sort_values("opt_n_steps")
        # Leave out runs near the round-off limit, where the step size no longer matters.
        ok = g["exit"].isin(["left_outward", "left_inward"]) & (g["orbits"] < ROUNDOFF_ORBITS - 0.5)
        k = np.polyfit(np.log10(g.loc[ok, "steps_per_orbit"]), g.loc[ok, "orbits"], 1)
        order = {"RK4": 4, "GL2": 4, "GL3": 6, "Tao4": 4}[m]
        rows.append(
            {
                "method": m,
                "formulation": form,
                "setting": "step size",
                "orbits_tightest": g["orbits"].iloc[-1],
                "orbits_loosest": g["orbits"].iloc[0],
                "orbits_per_decade": k[0],
                "fit_points": int(ok.sum()),
                "predicted_per_decade": order * PER_DECADE,
            }
        )
    return pd.DataFrame(rows)


def roundoff_line(ax, offset=0.08) -> None:
    ax.axhline(ROUNDOFF_ORBITS, color=style.MUTED, linewidth=0.8)
    ax.annotate(
        rf"$\delta_0 = \epsilon_\mathrm{{mach}}$: {ROUNDOFF_ORBITS:.1f} orbits",
        (0.02, ROUNDOFF_ORBITS + offset),
        xycoords=("axes fraction", "data"),
        color=style.INK_2,
        fontsize=7,
    )


def plot(adaptive, fixed, offsets, traces) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(9.6, 7.0), constrained_layout=True)
    ax = axes[0, 0]
    shades = ["#86b6ef", "#3987e5", "#1c5cab", "#0d366b"]
    for (tol, (t, dr)), color in zip(traces.items(), shades, strict=True):
        ax.plot(t, np.maximum(dr, 1e-17), color=color, linewidth=1.1, label=f"tol {tol:g}")
    t = np.linspace(0, 200, 50)
    ax.plot(
        t,
        1e-15 * np.exp(t / analytic.B_CRIT),
        color=style.MUTED,
        linestyle=":",
        label=r"$\propto e^{t/(3\sqrt{3})}$",
    )
    ax.axhline(0.1, color=style.MUTED, linewidth=0.8)
    ax.set(
        yscale="log",
        xlabel="coordinate time t (M)",
        ylabel=r"$|r - 3|$",
        title="(a) DOP853, inclined 45°, starting at r = 3",
        ylim=(1e-17, 1),
    )
    ax.legend(loc="lower right")

    ax = axes[0, 1]
    for m in ("DP5", "DOP853"):
        for form in ("a", "b"):
            g = adaptive[(adaptive["method"] == m) & (adaptive["formulation"] == form)]
            g = g.sort_values("opt_rtol")
            ax.plot(
                g["opt_rtol"], g["orbits"], **style.method_style(m, form), label=f"{m} ({form})"
            )
    # The predicted slope, with the intercept fitted to DOP853 (b).
    g = adaptive[(adaptive["method"] == "DOP853") & (adaptive["formulation"] == "b")]
    intercept = np.mean(g["orbits"] + PER_DECADE * np.log10(g["opt_rtol"]))
    tol = np.array([1e-13, 1e-3])
    ax.plot(
        tol,
        intercept - PER_DECADE * np.log10(tol),
        color=style.MUTED,
        linestyle=":",
        label=r"slope $\ln 10 / 2\pi$ = 0.37 per decade",
    )
    roundoff_line(ax)
    ax.set(
        xscale="log",
        xlabel="rtol = atol",
        ylabel="orbits before |r - 3| > 0.1",
        title="(b) adaptive pairs",
    )
    ax.legend(loc="lower left", ncol=2)

    ax = axes[1, 0]
    for m in FIXED:
        for form in ("b",) if m.startswith("Tao") else ("a", "b"):
            g = fixed[(fixed["method"] == m) & (fixed["formulation"] == form)]
            g = g.sort_values("steps_per_orbit")
            ax.plot(
                g["steps_per_orbit"],
                g["orbits"],
                **style.method_style(m, form),
                label=f"{m} ({form})",
            )
    roundoff_line(ax, offset=-0.3)
    for order in (4, 6):
        x = np.array([8, 80])
        y = 3.4 + order * PER_DECADE * np.log10(x / x[0])
        ax.plot(x, y, color=style.MUTED, linestyle=":")
        ax.annotate(f"{order} × 0.37", (x[1], y[1]), color=style.INK_2, fontsize=7, va="center")
    ax.set(
        xscale="log",
        xlabel="steps per photon-sphere orbit",
        ylabel="orbits before |r - 3| > 0.1",
        title="(c) fixed-step methods",
    )
    ax.legend(loc="lower right", ncol=2)

    ax = axes[1, 1]
    o = offsets.sort_values("case_k")
    ax.plot(o["case_k"], o["exact"], color=style.INK, linewidth=1.5, label="exact (elliptic)")
    ax.plot(
        o["case_k"], o["orbits"], **style.method_style("DOP853", "b"), label="DOP853 (b), tol 1e-13"
    )
    ax.set(
        xlabel=r"$k$, with $r_0 = 3 + 10^{-k}$ (equatorial)",
        ylabel="orbits before |r - 3| > 0.1",
        title="(d) explicit offset",
    )
    ax.legend(loc="upper left")
    style.save(fig, "F10_photon_sphere")


def main():
    adaptive, fixed, offsets, traces = compute()
    table = survival_table(adaptive, fixed)
    save_table(table, "t6_photon_sphere")
    (SUMMARY / "t6_photon_sphere.md").write_text(markdown_table(table) + "\n")
    save_table(
        pd.concat([adaptive, fixed], ignore_index=True),
        "f10_photon_sphere_runs",
        ["formulation", "method", "opt_rtol", "opt_n_steps", "orbits", "exit", "nfev"],
    )
    save_table(offsets, "f10_photon_sphere_offsets", ["case_k", "orbits", "exact", "exit"])
    plot(adaptive, fixed, offsets, traces)
    return table


if __name__ == "__main__":
    print(main().to_string())
