"""Common machinery of the geodesic test cases."""

from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np

from ..formulations import Formulation, Hamiltonian
from ..integrators import Event, Problem, Solution

R_INDEX, THETA_INDEX, PHI_INDEX, T_INDEX = 1, 2, 3, 0
#: Radial velocity (a) or momentum (b); both vanish at the same points.
RADIAL_RATE_INDEX = 5
HORIZON_MARGIN = 1e-3


def periapsis_event() -> Event:
    """Radial turning point with ``r`` increasing afterwards (``u^r``, ``p_r`` cross 0 upward)."""
    return Event("periapsis", RADIAL_RATE_INDEX, 0.0, +1)


def capture_event(metric) -> Event:
    (surface,) = metric.stop_surfaces(HORIZON_MARGIN)
    return Event(surface.name, surface.index, surface.value, surface.direction, stop_after=1)


@dataclass(frozen=True)
class GeodesicCase:
    """A test case: initial data, integration span, events, exact reference, error function.

    Subclasses are frozen dataclasses whose fields are the case parameters, so that a case is
    fully described by ``(family, params)`` and can key the experiment cache.
    """

    family = "case"
    eps = 1

    @property
    def id(self) -> str:
        params = "-".join(f"{k}{v:g}" for k, v in asdict(self).items())
        return f"{self.family}-{params}"

    @property
    def params(self) -> dict:
        return asdict(self)

    def initial_xp(self, metric) -> tuple[np.ndarray, np.ndarray]:
        raise NotImplementedError

    @property
    def lam_end(self) -> float:
        raise NotImplementedError

    def events(self, metric) -> tuple[Event, ...]:
        return (capture_event(metric),)

    def problem(self, formulation: Formulation) -> Problem:
        x, p = self.initial_xp(formulation.metric)
        return Problem(
            formulation.rhs,
            formulation.from_xp(x, p),
            (0.0, self.lam_end),
            self.events(formulation.metric),
            hamiltonian=isinstance(formulation, Hamiltonian),
        )

    def reference(self) -> dict:
        return {}

    def errors(self, solution: Solution, formulation: Formulation) -> dict:
        """Errors against the exact reference; ``"error"`` is the headline number."""
        raise NotImplementedError

    def constraint_errors(self, solution: Solution, formulation: Formulation) -> dict:
        dH = formulation.relative_constraint(solution.y, self.eps)
        return {"dH_max": float(np.max(dH)), "dH_end": float(dH[-1])}
