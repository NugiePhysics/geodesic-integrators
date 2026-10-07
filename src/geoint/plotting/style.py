"""One look for every figure: fixed colour and marker per method, recessive chrome.

Colours are the eight categorical slots of the validated reference palette (dataviz skill),
assigned to methods in a fixed order so a method keeps its colour in every figure. Every
series also has its own marker, because three slots are below 3:1 contrast on white and
colour must never be the only cue. SciPy runs reuse their method's colour with hollow markers.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[3]
FIGURES = ROOT / "figures"

INK = "#0b0b0b"
INK_2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"

METHODS = {
    "RK4": {"color": "#2a78d6", "marker": "o"},
    "DP5": {"color": "#eb6834", "marker": "s"},
    "DOP853": {"color": "#1baf7a", "marker": "D"},
    "GL1": {"color": "#eda100", "marker": "v"},
    "GL2": {"color": "#e87ba4", "marker": "^"},
    "GL3": {"color": "#008300", "marker": "<"},
    "Tao2": {"color": "#4a3aa7", "marker": "P"},
    "Tao4": {"color": "#e34948", "marker": "X"},
}
METHODS["scipy-RK45"] = {**METHODS["DP5"], "fillstyle": "none"}
METHODS["scipy-DOP853"] = {**METHODS["DOP853"], "fillstyle": "none"}

FORMULATIONS = {
    "a": {"linestyle": "--", "label": "(a) second order"},
    "b": {"linestyle": "-", "label": "(b) Hamiltonian"},
}


def method_style(method: str, formulation: str = "b", **overrides) -> dict:
    style = {
        "color": METHODS[method]["color"],
        "marker": METHODS[method]["marker"],
        "fillstyle": METHODS[method].get("fillstyle", "full"),
        "linestyle": FORMULATIONS[formulation]["linestyle"],
        "linewidth": 1.4,
        "markersize": 4.5,
    }
    style.update(overrides)
    return style


def apply() -> None:
    mpl.rcParams.update(
        {
            "figure.dpi": 110,
            "savefig.dpi": 200,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "font.family": "sans-serif",
            "font.size": 9,
            "axes.titlesize": 9.5,
            "axes.labelsize": 9,
            "legend.fontsize": 7.5,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "text.color": INK,
            "axes.labelcolor": INK_2,
            "axes.edgecolor": AXIS,
            "axes.linewidth": 0.8,
            "xtick.color": MUTED,
            "ytick.color": MUTED,
            "xtick.labelcolor": INK_2,
            "ytick.labelcolor": INK_2,
            "axes.grid": True,
            "grid.color": GRID,
            "grid.linewidth": 0.6,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "legend.frameon": False,
            "lines.linewidth": 1.4,
            "lines.markersize": 4.5,
            "axes.prop_cycle": mpl.cycler(color=[m["color"] for m in METHODS.values()][:8]),
            "mathtext.fontset": "dejavusans",
        }
    )


def save(fig, name: str) -> list[Path]:
    """Write ``figures/<name>.pdf`` (report) and ``figures/<name>.png`` (README)."""
    FIGURES.mkdir(exist_ok=True)
    paths = [FIGURES / f"{name}.pdf", FIGURES / f"{name}.png"]
    for path in paths:
        # A fixed creation date keeps the PDFs byte-identical between runs.
        metadata = {"CreationDate": None} if path.suffix == ".pdf" else None
        fig.savefig(path, bbox_inches="tight", metadata=metadata)
    plt.close(fig)
    return paths


def reference_slope(ax, x0, y0, slope, length=1.5, label=None, **kwargs):
    """Draw a short ``y ~ x^slope`` guide starting at ``(x0, y0)`` (log-log axes)."""
    x1 = x0 * 10**length
    y1 = y0 * (x1 / x0) ** slope
    ax.plot([x0, x1], [y0, y1], color=MUTED, linewidth=1.0, linestyle=":", **kwargs)
    if label:
        ax.annotate(
            label,
            (x1, y1),
            color=INK_2,
            fontsize=7,
            xytext=(3, 0),
            textcoords="offset points",
            va="center",
        )
