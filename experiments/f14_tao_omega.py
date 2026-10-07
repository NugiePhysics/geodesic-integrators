"""F14: Tao's method against its coupling frequency omega, and the separation of its copies.

Left: error after 3 radial periods of the orbit (p, e) = (20, 0.5) and after the deflection of
a photon with b = 6, against omega h, at two step sizes. Right: separation |q_r - x_r| of the
two copies over 300 radial periods for several omega.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from _common import save_table

from geoint.experiments import RunSpec, run_many, solve
from geoint.experiments.analysis import envelope
from geoint.plotting import style
from geoint.testcases.schwarzschild import Deflection, Eccentric

ORBIT = Eccentric(20.0, 0.5, 3.0)
PHOTON = Deflection(6.0)
OMEGA_H = np.geomspace(1e-6, 3.0, 22)
LONG = Eccentric(20.0, 0.5, 300.0)
LONG_OMEGA = (1e-5, 1e-4, 1e-3, 1e-2)


def compute():
    T = ORBIT.reference()["T_tau"]
    settings = [(ORBIT, int(3 * n)) for n in (400, 1600)] + [(PHOTON, n) for n in (2000, 8000)]
    specs = []
    for case, n in settings:
        h = case.lam_end / n
        for method in ("Tao2", "Tao4"):
            specs += [RunSpec.make(case, "b", method, n_steps=n, omega=wh / h) for wh in OMEGA_H]
    df = run_many(specs)
    df["omega_h"] = df["opt_omega"] * df["h"]
    df["error"] = pd.to_numeric(df["error"], errors="coerce")
    separation = {}
    for omega in LONG_OMEGA:
        n = int(300 * 800)
        sol, _ = solve(RunSpec.make(LONG, "b", "Tao4", n_steps=n, omega=omega, save_every=40))
        z = sol.extended
        lam, sep = envelope(sol.lam, np.abs(z[:, 1] - z[:, 9]), T)
        separation[omega] = (lam, sep)
    return df, separation


def plot(df, separation):
    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.4), constrained_layout=True)
    panels = (
        (axes[0], ORBIT, "orbit (20, 0.5), 3 periods\n400 (coarse) and 1600 (fine) steps/period"),
        (axes[1], PHOTON, "deflection b = 6\n2000 (coarse) and 8000 (fine) steps"),
    )
    for ax, case, title in panels:
        sub = df[df["case"] == case.id]
        fine = sub["opt_n_steps"].max()
        for (method, n), g in sub.groupby(["method", "opt_n_steps"]):
            g = g.sort_values("omega_h")
            kind = "fine" if n == fine else "coarse"
            ax.plot(
                g["omega_h"],
                g["error"].replace(np.inf, np.nan),
                **style.method_style(method, linestyle="-" if n == fine else "--", markersize=3.5),
                label=f"{method}, {kind} step",
            )
        x = 1e-3 * case.lam_end / fine
        ax.axvline(x, color=style.MUTED, linewidth=0.8)
        ax.annotate(
            r"$\omega = 10^{-3}$",
            (x, 0.02),
            xycoords=("data", "axes fraction"),
            xytext=(3, 0),
            textcoords="offset points",
            color=style.INK_2,
            fontsize=7,
        )
        ax.set(xscale="log", yscale="log", xlabel=r"$\omega h$", title=title)
    axes[0].set_ylabel("error (rad)")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside lower center", ncol=4)
    ax = axes[2]
    shades = ["#86b6ef", "#3987e5", "#1c5cab", "#0d366b"]
    for (omega, (lam, sep)), color in zip(separation.items(), shades, strict=True):
        orbits = lam / ORBIT.reference()["T_tau"]
        ax.plot(orbits, sep, color=color, linewidth=1.2, label=rf"$\omega$ = {omega:g}")
    ax.set(
        xscale="log",
        yscale="log",
        xlabel="radial periods",
        ylabel=r"max $|q_r - x_r|$ per period",
        title="Tao4, 800 steps/period:\nseparation of the copies",
    )
    ax.legend(loc="upper left")
    style.save(fig, "F14_tao_omega")


def main():
    df, separation = compute()
    save_table(
        df,
        "f14_tao_omega",
        ["case", "method", "opt_n_steps", "h", "opt_omega", "omega_h", "error", "nfev", "status"],
    )
    plot(df, separation)
    return df, separation


if __name__ == "__main__":
    main()
