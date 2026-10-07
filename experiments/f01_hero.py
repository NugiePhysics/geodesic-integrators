"""F1: the geodesics of the study in the equatorial plane (README and report title figure).

Left: a parallel bundle of light rays coming from the left, with impact parameters from 2 b_c
down to b_c (1 + 10^-6); the last one winds about twice around the photon sphere before it
escapes. Each ray starts at r = 40 and is rotated by the exact angle it would have swept
coming in from infinity, so that the bundle is parallel. Right: six radial periods of the
orbit (p, e) = (12, 0.5). Both are integrated with fine fixed steps (RK4) in formulation (b)
and drawn with the horizon (r = 2M), the photon sphere (3M) and the ISCO (6M).
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Circle

from geoint import analytic
from geoint.experiments import RunSpec, solve
from geoint.plotting import style
from geoint.testcases.schwarzschild import Deflection, Eccentric

style.apply()

R_FAR = 40.0
RAYS = [2.0, 1.4, 1.15, 1.05, 1.01]  # b / b_c
CRITICAL = 1.0 + 1e-6
ORBIT = Eccentric(12.0, 0.5, 6.0)
BLUES = ["#86b6ef", "#5a9be8", "#3987e5", "#1c5cab", "#0d366b"]
ACCENT = "#eb6834"


def trajectory(case, n_steps: int, rotate: float = 0.0):
    """Cartesian (x, y) along a fixed-step RK4 run fine enough to draw smooth curves."""
    sol, _ = solve(RunSpec.make(case, "b", "RK4", n_steps=n_steps))
    r, phi = sol.y[:, 1], sol.y[:, 3] + rotate
    return r * np.cos(phi), r * np.sin(phi)


def compute():
    rays = {}
    for f in [*RAYS, CRITICAL]:
        b = f * analytic.B_CRIT
        # Angle swept from infinity to r_far on the way in: the incoming asymptote lies at
        # phi = -swept. Rotating by pi + swept puts it on the negative x axis.
        swept = (analytic.delta_phi(b) - analytic.delta_phi(b, R_FAR)) / 2
        rays[f] = trajectory(Deflection(b, R_FAR), 40_000, rotate=np.pi + swept)
    orbit = trajectory(ORBIT, 60_000)
    return rays, orbit


def landmarks(ax) -> list:
    """Horizon, photon sphere and ISCO; returns legend handles."""
    ax.add_patch(Circle((0, 0), 2.0, color=style.INK, zorder=3))
    handles = [Line2D([], [], color=style.INK, marker="o", linestyle="none", label="horizon 2M")]
    for r, ls, name in ((3.0, "--", "photon sphere 3M"), (6.0, ":", "ISCO 6M")):
        ax.add_patch(Circle((0, 0), r, fill=False, color=style.MUTED, linestyle=ls, lw=0.9))
        handles.append(Line2D([], [], color=style.MUTED, linestyle=ls, lw=0.9, label=name))
    ax.set_aspect("equal")
    ax.grid(False)
    ax.set_xticks([])
    ax.set_yticks([])
    for side in ("left", "bottom"):
        ax.spines[side].set_visible(False)
    return handles


def plot(rays, orbit) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.6), constrained_layout=True)
    ax = axes[0]
    for (f, (x, y)), color in zip(rays.items(), [*BLUES, ACCENT], strict=True):
        label = r"$b = b_c(1 + 10^{-6})$" if f == CRITICAL else f"$b = {f:g}\\,b_c$"
        ax.plot(x, y, color=color, lw=1.3 if f == CRITICAL else 1.0, label=label)
    landmarks(ax)
    ax.set(xlim=(-16, 16), ylim=(-15, 13), title="Light rays near the photon sphere")
    ax.legend(loc="upper left", fontsize=7)

    ax = axes[1]
    ax.plot(*orbit, color=style.INK, lw=0.9)
    handles = landmarks(ax)
    ax.legend(handles=handles, loc="lower right", fontsize=7)
    ax.set(
        xlim=(-26, 26),
        ylim=(-26, 26),
        title=f"Orbit (p, e) = (12, 0.5), six periods\n"
        f"periapsis advance {ORBIT.reference()['precession']:.2f} rad per period",
    )
    style.save(fig, "F01_hero")


def main():
    plot(*compute())


if __name__ == "__main__":
    main()
