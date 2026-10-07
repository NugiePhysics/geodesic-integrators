"""Every integrator of the study under one name (table T1)."""

from __future__ import annotations

from .base import Integrator
from .embedded_rk import DOP853_INTEGRATOR, DP5_INTEGRATOR
from .explicit_rk import ExplicitRK
from .implicit_rk import GL1_INTEGRATOR, GL2_INTEGRATOR, GL3_INTEGRATOR
from .scipy_adapter import ScipyIntegrator
from .tao import TAO2, TAO4

INTEGRATORS: dict[str, Integrator] = {
    "RK4": ExplicitRK(),
    "DP5": DP5_INTEGRATOR,
    "DOP853": DOP853_INTEGRATOR,
    "GL1": GL1_INTEGRATOR,
    "GL2": GL2_INTEGRATOR,
    "GL3": GL3_INTEGRATOR,
    "Tao2": TAO2,
    "Tao4": TAO4,
    "scipy-RK45": ScipyIntegrator("RK45"),
    "scipy-DOP853": ScipyIntegrator("DOP853"),
}

FIXED_STEP = ("RK4", "GL1", "GL2", "GL3", "Tao2", "Tao4")
ADAPTIVE = ("DP5", "DOP853")


def get(name: str) -> Integrator:
    try:
        return INTEGRATORS[name]
    except KeyError:
        raise KeyError(f"unknown integrator {name!r}; known: {sorted(INTEGRATORS)}") from None
