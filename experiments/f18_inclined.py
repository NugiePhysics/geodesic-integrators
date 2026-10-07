"""F18: inclined orbits (TC6) over 1000 radial periods: L^2 and the orbital plane.

The orbit (20, 0.5) in planes inclined by 30, 60 and 85 degrees. In formulation (b), L^2 =
p_theta^2 + p_phi^2 / sin^2 theta is a nonlinear invariant, the Schwarzschild analogue of
Carter's constant: no method conserves it exactly, and how its error grows is the preview of
Kerr. The tilt is the angle between the angular-momentum vector (L_x, L_y, L_z) at each
saved point (about four per period) and at the start.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from _common import save_table

from geoint.experiments import RunSpec
from geoint.experiments.analysis import classify_growth, fit_slope, log_binned_max
from geoint.experiments.longterm import series_many
from geoint.experiments.runner import formulation
from geoint.experiments.sweeps import TAO_OMEGA
from geoint.integrators import get
from geoint.plotting import style
from geoint.testcases.schwarzschild import Inclined, angular_momentum_vector

N_ORBITS = 1000
INCLINATIONS = (30.0, 60.0, 85.0)
METHODS = ("RK4", "DOP853", "GL2", "GL3", "Tao4")
STEPS_PER_ORBIT = 1000


def spec(method, form, inclination):
    case = Inclined(20.0, 0.5, inclination, float(N_ORBITS))
    if get(method).is_adaptive:
        opts = {"rtol": 1e-10, "atol": 1e-10, "save_every": 1}
    else:
        opts = {"n_steps": N_ORBITS * STEPS_PER_ORBIT, "save_every": STEPS_PER_ORBIT // 20}
        if get(method).requires_hamiltonian:
            opts["omega"] = TAO_OMEGA
    return RunSpec.make(case, form, method, **opts)


def keys():
    return [
        (m, f, i)
        for i in INCLINATIONS
        for m in METHODS
        for f in (("b",) if m.startswith("Tao") else ("a", "b"))
    ]


def compute():
    T = Inclined(20.0, 0.5, 60.0, 1.0).orbit.periods["T_tau"]
    ks = keys()
    data = series_many([spec(*k) for k in ks], T)
    out = {}
    for (m, f, i), s in zip(ks, data, strict=True):
        s = dict(s)
        form = formulation(f)
        L = np.array([angular_momentum_vector(*form.to_xp(y)) for y in s["traj_y"]])
        cross = np.linalg.norm(np.cross(L[0], L), axis=1)
        s["tilt"] = np.arctan2(cross, L @ L[0])
        # Longitude of the ascending node and inclination of the plane: with L_z and L^2
        # (nearly) conserved, the plane can only turn about the z axis.
        s["node"] = np.unwrap(np.arctan2(L[:, 0], -L[:, 1]))
        s["incl"] = np.arccos(L[:, 2] / np.linalg.norm(L, axis=1))
        s["tilt_orbits"] = s["traj_lam"] / T
        s["orbits_env"] = s["lam_env"] / T
        out[m, f, i] = s
    return out


def table(data) -> pd.DataFrame:
    rows = []
    for (m, f, i), s in data.items():
        late = s["orbits_env"] >= N_ORBITS / 10
        p, _, _ = fit_slope(s["orbits_env"][late], s["dL2_env"][late], lo=0, hi=np.inf)
        late = s["tilt_orbits"] >= N_ORBITS / 10
        p_tilt, _, _ = fit_slope(s["tilt_orbits"][late], s["tilt"][late], lo=0, hi=np.inf)
        orbits = s["tilt_orbits"][-1]
        rows.append(
            {
                "method": m,
                "formulation": f,
                "inclination": i,
                "status": str(s["status"]),
                "orbits": orbits,
                "dL2_max": s["dL2_env"].max(),
                "dL2_exponent": p,
                "dL2_growth": classify_growth(p),
                "tilt_end": s["tilt"][-1],
                "tilt_exponent": p_tilt,
                "node_drift_per_orbit": abs(s["node"][-1] - s["node"][0]) / orbits,
                "incl_change": abs(s["incl"][-1] - s["incl"][0]),
            }
        )
    return pd.DataFrame(rows)


def plot(data, t: pd.DataFrame) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(13.2, 3.9), constrained_layout=True)
    for m in METHODS:
        for f in ("b",) if m.startswith("Tao") else ("a", "b"):
            s = data[m, f, 60.0]
            st = style.method_style(m, f, markevery=4)
            x, y = log_binned_max(s["orbits_env"], np.maximum(s["dL2_env"], 1e-17))
            axes[0].plot(x, y, **st, label=f"{m} ({f})")
            x, y = log_binned_max(s["tilt_orbits"], np.maximum(s["tilt"], 1e-17))
            axes[1].plot(x, y, **st, label=f"{m} ({f})")
            g = t[(t["method"] == m) & (t["formulation"] == f)].sort_values("inclination")
            axes[2].plot(g["inclination"], g["node_drift_per_orbit"], **{**st, "markevery": None})
            failed = g[g["status"] != "completed"]
            axes[2].plot(
                failed["inclination"],
                failed["node_drift_per_orbit"],
                linestyle="none",
                marker="o",
                markersize=11,
                fillstyle="none",
                color=style.INK_2,
            )
    axes[0].set(
        xscale="log",
        yscale="log",
        xlabel="radial periods",
        ylabel=r"max $|\Delta L^2| / L^2$ per period",
        title="(i) inclination 60°: $L^2$",
    )
    axes[1].set(
        xscale="log",
        yscale="log",
        xlabel="radial periods",
        ylabel="tilt of the orbital plane (rad)",
        title="(ii) inclination 60°: orbital plane",
        xlim=(0.7, 1500),
    )
    axes[2].set(
        yscale="log",
        xlabel="inclination (degrees)",
        ylabel="spurious node precession (rad per period)",
        title="(iii) node drift against inclination",
        xticks=INCLINATIONS,
    )
    axes[2].annotate(
        "circled: Tao4's copies separated,\nthe orbit was captured",
        (0.5, 0.04),
        xycoords="axes fraction",
        fontsize=7,
        color=style.INK_2,
    )
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside lower center", ncol=9)
    style.save(fig, "F18_inclined")


def main():
    data = compute()
    t = table(data)
    save_table(t, "f18_inclined")
    plot(data, t)
    return t


if __name__ == "__main__":
    print(main().to_string())
