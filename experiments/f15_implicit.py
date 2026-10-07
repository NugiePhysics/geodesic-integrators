"""F15: cost and accuracy of the implicit iteration of the Gauss-Legendre methods.

(i) Fixed-point iterations per step against h (orbit (7.5, 0.5), both formulations, from F4).
(ii) Iterations per step along the orbit, against r (h fixed), formulation (b).
(iii) Relative energy error |dH| over 1000 radial periods of the orbit (20, 0.5) for GL2 with
the iteration stopped at a tolerance, against iteration to round-off stagnation (pitfall 10).
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from _common import SUMMARY, save_table

from geoint.experiments import RunSpec, solve
from geoint.experiments.analysis import envelope
from geoint.plotting import style
from geoint.testcases.schwarzschild import Eccentric

ORBIT = Eccentric(7.5, 0.5, 1.0)
LONG = Eccentric(20.0, 0.5, 1000.0)
ITER_TOLS = (1e-6, 1e-8, 1e-10, 1e-12, 0.0)
STEPS_PER_ORBIT = 400


def along_orbit():
    out = {}
    for method in ("GL1", "GL2", "GL3"):
        sol, _ = solve(RunSpec.make(ORBIT, "b", method, n_steps=2000))
        out[method] = (sol.y[1:, 1], np.diff(sol.counts[:, 1]))
    return out


def drift():
    T = LONG.reference()["T_tau"]
    out = {}
    for tol in ITER_TOLS:
        spec = RunSpec.make(
            LONG,
            "b",
            "GL2",
            n_steps=int(1000 * STEPS_PER_ORBIT),
            iter_tol=tol,
            save_every=STEPS_PER_ORBIT // 10,
        )
        sol, form = solve(spec)
        dH = form.relative_constraint(sol.y, 1)
        lam, env = envelope(sol.lam, dH, T)
        out[tol] = (lam / T, env, sol.n_iter / sol.n_steps)
    return out


def plot(conv: pd.DataFrame, orbit, long) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(10.8, 3.4), constrained_layout=True)
    ax = axes[0]
    for method in ("GL1", "GL2", "GL3"):
        for form in ("a", "b"):
            g = conv[(conv["method"] == method) & (conv["formulation"] == form)].sort_values("h")
            g = g[g["status"] == "completed"]
            ax.plot(
                g["h"],
                g["iter_per_step"],
                **style.method_style(method, form),
                label=f"{method} ({form})",
            )
    ax.set(
        xscale="log",
        xlabel="step size h",
        ylabel="fixed-point iterations per step",
        title="(i) orbit (7.5, 0.5)",
    )
    ax.legend(loc="upper left", ncol=2)
    ax = axes[1]
    for method, (r, its) in orbit.items():
        st = style.method_style(method, "b")
        ax.plot(
            r,
            its,
            linestyle="none",
            marker=st["marker"],
            color=st["color"],
            markersize=2.5,
            alpha=0.6,
            label=method,
        )
    ax.set(
        xlabel="r (M)",
        ylabel="iterations in that step",
        title="(ii) along the orbit, 2000 steps/period",
    )
    ax.legend(loc="upper right")
    ax = axes[2]
    shades = ["#b7d3f6", "#86b6ef", "#5598e7", "#256abf", style.INK]
    for (tol, (orbits, env, its)), color in zip(long.items(), shades, strict=True):
        label = "stagnation" if tol == 0 else f"tol {tol:g}"
        ax.plot(orbits, env, color=color, linewidth=1.2, label=f"{label} ({its:.1f} it/step)")
    ax.set(
        xscale="log",
        yscale="log",
        xlabel="radial periods",
        ylabel=r"max $|\delta H|$ per period",
        title="(iii) GL2, orbit (20, 0.5), 400 steps/period",
    )
    ax.legend(loc="upper left")
    style.save(fig, "F15_implicit")


def main():
    conv = pd.read_csv(SUMMARY / "f04_convergence.csv")
    conv = conv[(conv["case_label"] == "orbit") & conv["method"].str.startswith("GL")]
    orbit = along_orbit()
    long = drift()
    rows = []
    for tol, (_, env, its) in long.items():
        tenth = max(1, len(env) // 10)
        rows.append(
            {
                "iter_tol": tol,
                "iter_per_step": its,
                "dH_first_10pct": env[:tenth].max(),
                "dH_last_10pct": env[-tenth:].max(),
            }
        )
    save_table(pd.DataFrame(rows), "f15_iteration_tolerance")
    plot(conv, orbit, long)
    return conv, orbit, long


if __name__ == "__main__":
    main()
