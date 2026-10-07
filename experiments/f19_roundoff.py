"""F19: round-off over 10^4 orbits, with and without compensated summation (pitfalls 5, 12).

(i) The circular orbit r_c = 10 in the equatorial plane: r and p_r are an exact fixed point of
every Runge-Kutta map, so the only error left is round-off, which accumulates in the
cyclic coordinate phi (phase error against the exact phi = (L / r_c^2) lambda).
(ii) The mass-shell error of GL3 on the orbit (20, 0.5): with 1000 steps per period its
truncation error in H is below the round-off level, so Brouwer's law (|dH| ~ N^1/2) or a
linear round-off drift becomes visible, depending on the summation.
(iii) The inclined photon-sphere orbit of F10 against the step size. With plain summation an
increment of r smaller than half an ulp of 3 (2.2e-16) is lost, r stays exactly at 3, and the
radial force, which is proportional to r - 3, stays zero: rounding freezes the instability.
The smaller the step, the smaller the increments and the longer the orbit appears to survive.
Compensated summation keeps the lost low-order bits and restores the physical behaviour.
"""

from __future__ import annotations

import math

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from _common import save_table

from geoint import analytic
from geoint.experiments import RunSpec, run_many
from geoint.experiments.analysis import fit_slope, log_binned_max
from geoint.experiments.longterm import series_many
from geoint.plotting import style
from geoint.testcases.schwarzschild import Circular, Eccentric, PhotonSphere

N = 10_000
CIRCULAR = Circular(10.0, float(N))
ECCENTRIC = Eccentric(20.0, 0.5, float(N))
STEPS = 1000
PHOTON = PhotonSphere(math.inf, 45.0)
PHOTON_ORBIT = 2 * math.pi * 9 / analytic.B_CRIT  # affine length of one orbit
PHOTON_STEPS = [250, 500, 1000, 2000, 4000, 8000, 16000]


def compute():
    T_c = CIRCULAR.lam_end / N
    circ_keys = [(m, comp) for m in ("RK4", "GL2", "GL3") for comp in (True, False)]
    circ = series_many(
        [
            RunSpec.make(
                CIRCULAR, "b", m, n_steps=N * STEPS, save_every=STEPS // 20, compensated=comp
            )
            for m, comp in circ_keys
        ],
        T_c,
    )
    ecc_keys = [(m, comp) for m in ("GL3", "RK4") for comp in (True, False)]
    T_e = ECCENTRIC.reference()["T_tau"]
    ecc = series_many(
        [
            RunSpec.make(
                ECCENTRIC, "b", m, n_steps=N * STEPS, save_every=STEPS // 20, compensated=comp
            )
            for m, comp in ecc_keys
        ],
        T_e,
    )
    L_over_r2 = CIRCULAR.reference()["L"] / CIRCULAR.r_c**2
    out_c = {}
    for k, s in zip(circ_keys, circ, strict=True):
        lam, y = s["traj_lam"], s["traj_y"]
        out_c[k] = (lam / T_c, np.abs(y[:, 3] - L_over_r2 * lam))
    out_e = {k: (s["lam_env"] / T_e, s["dH_env"]) for k, s in zip(ecc_keys, ecc, strict=True)}
    return out_c, out_e, photon_sphere()


def photon_sphere() -> pd.DataFrame:
    specs = []
    for m in ("RK4", "GL3"):
        for comp in (True, False):
            for spo in PHOTON_STEPS:
                n = int(round(PHOTON.lam_end / PHOTON_ORBIT * spo))
                opts = {"n_steps": n, "save_every": spo // 8, "compensated": comp}
                specs.append(RunSpec.make(PHOTON, "b", m, **opts))
    df = run_many(specs)
    df["steps_per_orbit"] = df["opt_n_steps"] / (PHOTON.lam_end / PHOTON_ORBIT)
    return df


def plot(circ, ecc, photon) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(13.2, 3.9), constrained_layout=True)
    for ax, data in zip(axes, (circ, ecc), strict=False):
        for (m, comp), (orbits, err) in data.items():
            x, y = log_binned_max(orbits, np.maximum(err, 1e-17))
            ax.plot(
                x,
                y,
                **style.method_style(m, "b" if comp else "a", markevery=4),
                label=f"{m}, {'compensated' if comp else 'plain'} sum",
            )
    for ax, (x0, y0, slopes) in zip(
        axes, ((30, 1e-15, (1, 2)), (10, 1e-15, (0.5, 1))), strict=False
    ):
        for p in slopes:
            length = 1.5 if p == 2 else 2.5
            style.reference_slope(ax, x0, y0, p, length=length, label=rf"$N^{{{p:g}}}$")
    axes[0].set(
        xscale="log",
        yscale="log",
        xlabel="orbits",
        ylabel=r"$|\phi - (L/r_c^2)\lambda|$",
        title="(i) circular r_c = 10: phase from round-off alone\n(RK4, GL2, GL3 coincide)",
    )
    axes[1].set(
        xscale="log",
        yscale="log",
        xlabel="radial periods",
        ylabel=r"max $|\delta H|$ per period",
        title="(ii) orbit (20, 0.5): mass shell",
    )
    ax = axes[2]
    for m in ("RK4", "GL3"):
        for comp in (True, False):
            g = photon[(photon["method"] == m) & (photon["opt_compensated"] == comp)]
            g = g.sort_values("steps_per_orbit")
            ax.plot(
                g["steps_per_orbit"],
                g["orbits"],
                **style.method_style(m, "b" if comp else "a"),
                label=f"{m}, {'compensated' if comp else 'plain'} sum",
            )
    cap = math.acosh(0.1 / np.finfo(float).eps) / (2 * math.pi)
    ax.axhline(cap, color=style.MUTED, linewidth=0.8)
    ax.annotate(
        rf"$\delta_0 = \epsilon_\mathrm{{mach}}$: {cap:.1f} orbits",
        (0.7, cap + 0.12),
        xycoords=("axes fraction", "data"),
        fontsize=7,
        color=style.INK_2,
    )
    ax.set(
        xscale="log",
        xlabel="steps per photon-sphere orbit",
        ylabel="orbits before |r - 3| > 0.1",
        title="(iii) photon sphere, inclined 45° (F10)",
    )
    axes[0].legend(loc="upper left")
    axes[1].legend(loc="upper left")
    axes[2].legend(loc="upper left")
    style.save(fig, "F19_roundoff")


def main():
    circ, ecc, photon = compute()
    rows = []
    for name, data in (("circular phase", circ), ("eccentric dH", ecc)):
        for (m, comp), (x, y) in data.items():
            late = x >= x.max() / 10
            p, _, _ = fit_slope(x[late], y[late], lo=1e-300, hi=np.inf)
            rows.append(
                {"quantity": name, "method": m, "compensated": comp, "final": y[-1], "exponent": p}
            )
    t = pd.DataFrame(rows)
    save_table(t, "f19_roundoff")
    save_table(
        photon,
        "f19_photon_sphere_summation",
        ["method", "opt_compensated", "steps_per_orbit", "orbits", "exit"],
    )
    plot(circ, ecc, photon)
    return t


if __name__ == "__main__":
    print(main().to_string())
