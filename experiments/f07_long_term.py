"""F7, F8, F9, F17 and T4: 10^4 radial periods of the orbit (p, e) = (20, 0.5) (TC5).

Every method in both formulations (Tao in (b) only): fixed-step methods with 1000 steps per
radial period, adaptive pairs at tol 1e-10. From each run: the maximum per period of the
relative errors in H (mass shell), E, L_z, and the phase error at every periapsis passage
against the exact Phi. T4 fits the growth exponents over the last decade of periods and
classifies them as bounded (0), random walk (1/2), linear (1) or quadratic (2). It also fits
the signed phase error as c1 N + c2 N^2: a constant frequency error plus a drifting one.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from _common import SUMMARY, markdown_table, save_table

from geoint.experiments import RunSpec
from geoint.experiments.analysis import (
    classify_growth,
    fit_phase_drift,
    fit_slope,
    log_binned_max,
)
from geoint.experiments.longterm import series_many
from geoint.experiments.sweeps import ALL, TAO_OMEGA, formulations_for
from geoint.integrators import get
from geoint.plotting import style
from geoint.testcases.schwarzschild import Eccentric

N_ORBITS = 10_000
CASE = Eccentric(20.0, 0.5, float(N_ORBITS))
STEPS_PER_ORBIT = 1000
TOL = 1e-10


def spec_for(method: str, form: str, **extra) -> RunSpec:
    if get(method).is_adaptive:
        options = {"rtol": TOL, "atol": TOL, "save_every": 1 if method == "DOP853" else 4}
    else:
        options = {"n_steps": N_ORBITS * STEPS_PER_ORBIT, "save_every": STEPS_PER_ORBIT // 20}
        if get(method).requires_hamiltonian:
            options["omega"] = TAO_OMEGA
    options.update(extra)
    return RunSpec.make(CASE, form, method, **options)


def runs():
    return [(m, f) for m in ALL for f in formulations_for(m)]


def compute():
    keys = runs()
    ref = CASE.reference()
    data = series_many([spec_for(m, f) for m, f in keys], ref["T_tau"])
    out = {}
    for (method, form), s in zip(keys, data, strict=True):
        k = np.arange(1, s["peri_phi"].size + 1)
        s = dict(s)
        s["orbits_env"] = s["lam_env"] / ref["T_tau"]
        s["peri_k"] = k
        s["phase_signed"] = s["peri_phi"] - k * ref["Phi"]
        s["phase_err"] = np.abs(s["phase_signed"])
        out[method, form] = s
    return out


def late_exponent(x, y) -> float:
    """Log-log slope over the last decade of ``x`` (all points, no error window)."""
    late = x >= x.max() / 10
    return fit_slope(x[late], y[late], lo=0, hi=np.inf)[0]


def growth_table(data) -> pd.DataFrame:
    rows = []
    for (method, form), s in data.items():
        p_dH = late_exponent(s["orbits_env"], s["dH_env"])
        p_ph = late_exponent(s["peri_k"], s["phase_err"])
        c1, c2, crossover = fit_phase_drift(s["peri_k"], s["phase_signed"])
        row = {
            "method": method,
            "formulation": form,
            "status": str(s["status"]),
            "nfev_per_orbit": float(s["nfev"]) / N_ORBITS,
            "dH_first_orbit": s["dH_env"][0],
            "dH_last_orbit": s["dH_env"][-1],
            "dH_exponent": p_dH,
            "dH_growth": classify_growth(p_dH),
            "phase_err_end": s["phase_err"][-1],
            "phase_exponent": p_ph,
            "phase_growth": classify_growth(p_ph),
            "phase_c1": c1,
            "phase_c2": c2,
            "phase_crossover": crossover,
            "dE_max": s["dE_env"].max(),
            "dLz_max": s["dLz_env"].max(),
        }
        for key in ("dE", "dLz"):
            values = s[f"{key}_env"]
            p = late_exponent(s["orbits_env"], values) if values.max() > 0 else np.nan
            row[f"{key}_exponent"] = p
        rows.append(row)
    return pd.DataFrame(rows)


def envelope_line(ax, x, y, method, form, label=None):
    """A per-period series, drawn as its maximum over logarithmic bins."""
    xb, yb = log_binned_max(x, np.maximum(y, 1e-18))
    ax.plot(xb, yb, **style.method_style(method, form, markevery=4), label=label)


def plot_dH(data):
    """F7: |dH| per radial period, formulation (b)."""
    fig, ax = plt.subplots(figsize=(6.6, 4.2), constrained_layout=True)
    for method in ALL:
        s = data[method, "b"]
        envelope_line(ax, s["orbits_env"], s["dH_env"], method, "b", label=method)
    for p, label in ((0, r"$N^0$"), (0.5, r"$N^{1/2}$"), (1, r"$N^1$")):
        style.reference_slope(ax, 10, 1e-14, p, length=2.5, label=label)
    ax.set(
        xscale="log",
        yscale="log",
        xlabel="radial periods N",
        ylabel=r"max $|\delta H|$ per period",
        title="Mass-shell error over 10⁴ periods of (20, 0.5), formulation (b)",
        ylim=(3e-17, 1e-3),
    )
    ax.legend(loc="upper left", ncol=4)
    style.save(fig, "F07_energy_long_term")


def plot_invariants(data):
    """F8: E, L_z and the norm in formulation (a); (b) keeps E and L_z bit for bit."""
    fig, axes = plt.subplots(1, 3, figsize=(10.8, 3.6), sharey=True, constrained_layout=True)
    panels = (
        ("dE_env", r"$|\Delta E| / E$"),
        ("dLz_env", r"$|\Delta L_z| / L_z$"),
        ("dH_env", r"$|g_{\mu\nu}u^\mu u^\nu + 1|$"),
    )
    for ax, (key, label) in zip(axes, panels, strict=True):
        for method in ALL:
            if (method, "a") not in data:
                continue
            s = data[method, "a"]
            envelope_line(ax, s["orbits_env"], s[key], method, "a", label=method)
        ax.set(xscale="log", yscale="log", xlabel="radial periods N", title=f"(a): {label}")
    axes[0].set_ylabel("max per period")
    exact = all(
        (s["dE_env"] == 0).all() and (s["dLz_env"] == 0).all()
        for (m, f), s in data.items()
        if f == "b"
    )
    note = (
        "In (b), E and L_z of every method stay exactly at their initial values "
        "(every sample, all 10⁴ periods)"
        if exact
        else "In (b), E and L_z drift"
    )
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside lower center", ncol=6)
    fig.suptitle(note, fontsize=8, color=style.INK_2)
    style.save(fig, "F08_invariants_formulation_a")
    return exact


def plot_phase(data):
    """F9: phase error at the N-th periapsis, both formulations."""
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 3.8), sharey=True, constrained_layout=True)
    for ax, form in zip(axes, ("a", "b"), strict=True):
        for method in ALL:
            if (method, form) not in data:
                continue
            s = data[method, form]
            envelope_line(ax, s["peri_k"], s["phase_err"], method, form, label=method)
        style.reference_slope(ax, 1, 1e-12, 1, length=2, label=r"$N^1$")
        style.reference_slope(ax, 1, 1e-12, 2, length=1.5, label=r"$N^2$")
        ax.set(
            xscale="log",
            yscale="log",
            xlabel="periapsis passage N",
            title=style.FORMULATIONS[form]["label"],
        )
    axes[0].set_ylabel(r"phase error $|\phi_N - N\Phi|$ (rad)")
    handles, labels = axes[1].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside lower center", ncol=8)
    style.save(fig, "F09_phase_error")


def plot_symmetric_vs_symplectic(data):
    """F17: GL2 in (a) (symmetric) and (b) (symplectic) against RK4."""
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 3.6), constrained_layout=True)
    for method in ("GL2", "RK4"):
        for form in ("a", "b"):
            s = data[method, form]
            label = f"{method} {style.FORMULATIONS[form]['label']}"
            envelope_line(axes[0], s["orbits_env"], s["dH_env"], method, form, label=label)
            envelope_line(axes[1], s["peri_k"], s["phase_err"], method, form, label=label)
    axes[0].set(
        xscale="log",
        yscale="log",
        xlabel="radial periods N",
        ylabel=r"max $|\delta H|$ per period",
        title="mass-shell error",
    )
    axes[1].set(
        xscale="log",
        yscale="log",
        xlabel="periapsis passage N",
        ylabel="phase error (rad)",
        title="phase error",
    )
    axes[0].legend(loc="upper left")
    a, b = data["GL2", "a"], data["GL2", "b"]
    c1 = [fit_phase_drift(s["peri_k"], s["phase_signed"])[0] for s in (a, b)]
    axes[0].annotate(
        f"GL2 (a) and (b) overlap: max |δH| = {a['dH_env'].max():.5g} and {b['dH_env'].max():.5g}",
        (0.03, 0.4),
        xycoords="axes fraction",
        fontsize=7,
        color=style.INK_2,
    )
    axes[1].annotate(
        f"GL2 phase drift: {c1[0]:+.2e} (a), {c1[1]:+.2e} (b) rad per orbit",
        (0.03, 0.92),
        xycoords="axes fraction",
        fontsize=7,
        color=style.INK_2,
    )
    style.save(fig, "F17_symmetric_vs_symplectic")


def main():
    data = compute()
    table = growth_table(data)
    save_table(table, "t4_growth")
    (SUMMARY / "t4_growth.md").write_text(markdown_table(table) + "\n")
    plot_dH(data)
    exact = plot_invariants(data)
    plot_phase(data)
    plot_symmetric_vs_symplectic(data)
    print("E and L_z exact in (b):", exact)
    return data, table


if __name__ == "__main__":
    _, table = main()
    with pd.option_context("display.width", 200):
        print(table.to_string())
