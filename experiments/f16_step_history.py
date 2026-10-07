"""F16: step size of the adaptive pairs along the orbit, against r.

DP5 and DOP853 at rtol = atol = 1e-10, formulation (b): the deflection of a photon with b = 6
(from r_far = 1000 to periapsis and back) and one radial period of the orbit (20, 0.5).
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from geoint.experiments import RunSpec, solve
from geoint.plotting import style
from geoint.testcases.schwarzschild import Deflection, Eccentric

CASES = {"deflection b = 6": Deflection(6.0), "orbit (20, 0.5)": Eccentric(20.0, 0.5, 1.0)}
TOL = 1e-10


def compute():
    out = {}
    for label, case in CASES.items():
        for method in ("DP5", "DOP853"):
            sol, _ = solve(RunSpec.make(case, "b", method, rtol=TOL, atol=TOL))
            h = np.diff(sol.lam)
            r = 0.5 * (sol.y[1:, 1] + sol.y[:-1, 1])
            incoming = np.diff(sol.y[:, 1]) < 0
            out[label, method] = (r, h, incoming, sol.n_steps)
    return out


def plot(data) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.4), constrained_layout=True)
    for ax, label in zip(axes, CASES, strict=True):
        for method in ("DP5", "DOP853"):
            r, h, incoming, n = data[label, method]
            st = style.method_style(method, "b", markersize=2.5, linewidth=0.8)
            ax.plot(r, h, color=st["color"], linewidth=0.8, alpha=0.9)
            ax.plot(
                r[:: max(1, len(r) // 60)],
                h[:: max(1, len(r) // 60)],
                linestyle="none",
                marker=st["marker"],
                color=st["color"],
                markersize=3,
                label=f"{method}, {n} steps",
            )
        ax.set(xscale="log", yscale="log", xlabel="r (M)", title=f"{label}, tol = {TOL:g}")
        ax.legend(loc="upper left")
    axes[0].set_ylabel("step size h")
    style.save(fig, "F16_step_history")


def main():
    data = compute()
    plot(data)
    return data


if __name__ == "__main__":
    main()
