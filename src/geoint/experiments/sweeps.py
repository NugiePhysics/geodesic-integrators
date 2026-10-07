"""Building blocks for parameter sweeps: grids, per-method options, spec lists."""

from __future__ import annotations

import numpy as np

from ..integrators import get
from .runner import RunSpec

FIXED = ("RK4", "GL1", "GL2", "GL3", "Tao2", "Tao4")
ADAPTIVE = ("DP5", "DOP853")
ALL = (*FIXED[:1], *ADAPTIVE, *FIXED[1:])  # RK4, DP5, DOP853, GL1-3, Tao2, Tao4
#: Coupling frequency of Tao's method in all geodesic experiments (units 1/M), chosen in F14:
#: within a factor ~3 of the most accurate value at every step size tested, and large enough
#: to hold the two copies together over 10^4 orbits.
TAO_OMEGA = 1e-3


def step_grid(n_min: int, n_max: int, count: int = 10) -> list[int]:
    return sorted({int(round(n)) for n in np.geomspace(n_min, n_max, count)})


def tolerance_grid(tightest: float = 1e-13, loosest: float = 1e-3, count: int = 11) -> list[float]:
    return [float(t) for t in np.geomspace(loosest, tightest, count)]


def formulations_for(method: str, formulations=("a", "b")) -> tuple[str, ...]:
    """Tao needs the canonical form (b)."""
    return ("b",) if get(method).requires_hamiltonian else tuple(formulations)


def specs(case, methods, steps, tolerances, formulations=("a", "b"), omega=TAO_OMEGA, **extra):
    """Every (formulation, method, setting) combination for ``case``."""
    out = []
    for method in methods:
        adaptive = get(method).is_adaptive
        for form in formulations_for(method, formulations):
            for setting in tolerances if adaptive else steps:
                options = dict(extra)
                if adaptive:
                    options.update(rtol=setting, atol=setting)
                else:
                    options["n_steps"] = int(setting)
                    if get(method).requires_hamiltonian:
                        options["omega"] = omega
                out.append(RunSpec.make(case, form, method, **options))
    return out
