"""Problem, event and solution types shared by every integrator.

Nothing here knows about metrics or geodesics (ADR 0001): a problem is an autonomous vector field
``rhs(y)``, an initial state, an interval of the independent variable ``lam`` and, optionally,
events defined on single state components. Integrators that need a canonical Hamiltonian
structure (Tao) additionally require ``y = (q, p)`` with ``rhs(y) = (dH/dp, -dH/dq)``.
"""

from __future__ import annotations

import time
from collections.abc import Callable
from dataclasses import dataclass, field

import numpy as np
from numpy.typing import NDArray

#: Driver status codes; a positive status ``k`` means "terminated by event ``k - 1``".
COMPLETED = 0
NONFINITE = -1
STEP_TOO_SMALL = -2
MAX_STEPS = -3
ITERATION_FAILED = -4

FAILURES = {
    NONFINITE: "failed: non-finite state",
    STEP_TOO_SMALL: "failed: step size too small",
    MAX_STEPS: "failed: maximum number of steps",
    ITERATION_FAILED: "failed: implicit iteration did not converge",
}


@dataclass(frozen=True)
class Event:
    """A crossing ``y[index] == value``, located to the order of the integrator.

    ``direction = +1`` counts only crossings where ``y[index]`` increases, ``-1`` only where it
    decreases, ``0`` both. With ``stop_after = k > 0`` the integration terminates at the k-th
    occurrence. A state that starts exactly on the surface does not count as a crossing.
    """

    name: str
    index: int
    value: float
    direction: int = 0
    stop_after: int = 0


@dataclass
class Problem:
    """``dy/dlam = rhs(y)`` on ``lam_span`` from ``y0``."""

    rhs: Callable[[NDArray[np.float64]], NDArray[np.float64]]
    y0: NDArray[np.float64]
    lam_span: tuple[float, float]
    events: tuple[Event, ...] = ()
    hamiltonian: bool = False

    def __post_init__(self):
        self.y0 = np.ascontiguousarray(self.y0, dtype=np.float64)
        lam0, lam1 = map(float, self.lam_span)
        if not lam1 > lam0:
            raise ValueError(f"only forward integration is supported, got lam_span {self.lam_span}")
        self.lam_span = (lam0, lam1)
        for event in self.events:
            if not 0 <= event.index < self.y0.size:
                raise ValueError(f"event {event.name!r} refers to component {event.index}")

    def event_arrays(self):
        """Events as the arrays the compiled drivers take."""
        ev = self.events
        return (
            np.array([e.index for e in ev], dtype=np.int64),
            np.array([e.value for e in ev], dtype=np.float64),
            np.array([e.direction for e in ev], dtype=np.int64),
            np.array([e.stop_after for e in ev], dtype=np.int64),
        )


@dataclass
class EventRecord:
    """All occurrences of one event: parameter values ``lam`` and states ``y``."""

    lam: NDArray[np.float64]
    y: NDArray[np.float64]


@dataclass
class Solution:
    """What an integrator returns. Costs are counted as defined in ADR 0003.

    ``lam``, ``y`` hold the saved points (the initial point, every ``save_every``-th step and
    the final point). For Tao's method ``y`` is the first copy ``(q, p)`` of the extended state,
    and ``extended`` holds the full state.
    """

    method: str
    lam: NDArray[np.float64]
    y: NDArray[np.float64]
    status: str
    nfev: int
    n_steps: int
    n_rejected: int = 0
    n_iter: int = 0
    njev: int = 0
    wall_time: float = 0.0
    events: dict[str, EventRecord] = field(default_factory=dict)
    extended: NDArray[np.float64] | None = None
    options: dict = field(default_factory=dict)
    #: Cumulative ``(nfev, n_iter)`` at each saved point (fixed-step methods), for per-step costs.
    counts: NDArray[np.int64] | None = None

    @property
    def success(self) -> bool:
        return not self.status.startswith("failed")

    @property
    def lam_end(self) -> float:
        return float(self.lam[-1])

    @property
    def y_end(self) -> NDArray[np.float64]:
        return self.y[-1]


def status_name(code: int, events: tuple[Event, ...]) -> str:
    if code == COMPLETED:
        return "completed"
    if code > 0:
        return events[code - 1].name
    return FAILURES[code]


def event_records(problem: Problem, ids, lams, ys) -> dict[str, EventRecord]:
    """Group the driver's flat event output by event name."""
    records = {}
    for k, event in enumerate(problem.events):
        mask = ids == k
        records[event.name] = EventRecord(lams[mask], ys[mask])
    return records


class Integrator:
    """Base class. Subclasses implement :meth:`solve` and set the class attributes."""

    name: str = "integrator"
    order: int = 0
    is_adaptive: bool = False
    is_symplectic: bool = False
    is_symmetric: bool = False
    is_implicit: bool = False
    requires_hamiltonian: bool = False
    #: Vector-field evaluations per step (per accepted attempt, or per iteration if implicit).
    evals_per_step: int = 0

    def solve(self, problem: Problem, **options) -> Solution:
        raise NotImplementedError

    def step(self, rhs, y, h: float, **options) -> NDArray[np.float64]:
        """One step of size ``h`` (either sign) from ``y``: the map whose geometric properties
        (symplecticity, symmetry) the tests check. Only for fixed-step methods."""
        raise NotImplementedError(f"{self.name} has no fixed one-step map")

    def properties(self) -> dict:
        """Row of table T1."""
        return {
            "method": self.name,
            "order": self.order,
            "implicit": self.is_implicit,
            "symplectic": self.is_symplectic,
            "symmetric": self.is_symmetric,
            "adaptive": self.is_adaptive,
            "evals_per_step": self.evals_per_step,
        }

    def __repr__(self) -> str:
        return f"{type(self).__name__}({self.name!r})"


def fixed_step_count(lam_span: tuple[float, float], h: float | None, n_steps: int | None):
    """Number of steps and the step size that divides ``lam_span`` exactly.

    Exactly one of ``h`` (the largest step allowed) and ``n_steps`` must be given.
    """
    span = lam_span[1] - lam_span[0]
    if (h is None) == (n_steps is None):
        raise ValueError("give exactly one of h and n_steps")
    if n_steps is None:
        if not h > 0:
            raise ValueError(f"h must be positive, got {h}")
        n_steps = max(1, int(np.ceil(span / h * (1 - 1e-12))))
    if n_steps < 1:
        raise ValueError(f"n_steps must be positive, got {n_steps}")
    return int(n_steps), span / n_steps


def run_fixed_step(kind, params, problem: Problem, h, n_steps, save_every, compensated):
    """Run the compiled fixed-step driver for method ``kind``.

    Returns the raw driver output, the step size actually used and the wall time.
    """
    from ._compiled import fixed

    n, h = fixed_step_count(problem.lam_span, h, n_steps)
    t0 = time.perf_counter()
    out = fixed(
        kind,
        problem.rhs,
        problem.y0,
        problem.lam_span[0],
        h,
        n,
        params,
        bool(compensated),
        *problem.event_arrays(),
        max(1, int(save_every)),
    )
    return out, h, time.perf_counter() - t0


def fixed_step_solution(name, problem, out, h, wall, options, extended_dim=None) -> Solution:
    """Wrap the raw output of a fixed-step driver. ``extended_dim`` splits off Tao's copy."""
    lam, y, status, nfev, niter, n_steps, ev_ids, ev_lams, ev_ys, counts = out
    extended = None
    if extended_dim is not None:
        extended, y, ev_ys = y, y[:, :extended_dim], ev_ys[:, :extended_dim]
    return Solution(
        method=name,
        lam=lam,
        y=y,
        status=status_name(status, problem.events),
        nfev=int(nfev),
        n_steps=int(n_steps),
        n_iter=int(niter),
        wall_time=wall,
        events=event_records(problem, ev_ids, ev_lams, ev_ys),
        extended=extended,
        options={"h": h, **options},
        counts=counts,
    )


def one_step(kind: int, params, rhs, y, h: float) -> NDArray[np.float64]:
    """Apply the compiled step of method ``kind`` once, without compensated summation."""
    from ._compiled import step

    y = np.ascontiguousarray(y, dtype=np.float64)
    return step(kind, rhs, y, np.zeros(y.size), float(h), params, False)[0]
