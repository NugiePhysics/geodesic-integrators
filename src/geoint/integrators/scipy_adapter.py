"""``scipy.integrate.solve_ivp`` behind the common interface: the correctness oracle.

Its wall time includes Python overhead and is not comparable with the compiled methods
(ADR 0003); its ``nfev`` is, and the compiled DP5 and DOP853 must reproduce it.
"""

from __future__ import annotations

import time

import numpy as np
from scipy.integrate import solve_ivp

from .base import EventRecord, Integrator, Problem, Solution


def _event_function(event):
    def g(_, y):
        return y[event.index] - event.value

    g.direction = event.direction
    g.terminal = event.stop_after if event.stop_after > 0 else False
    return g


class ScipyIntegrator(Integrator):
    """``solve_ivp`` with ``method`` ``"RK45"`` or ``"DOP853"``."""

    is_adaptive = True

    def __init__(self, method: str):
        self.method = method
        self.name = f"scipy-{method}"
        self.order = {"RK45": 5, "DOP853": 8}[method]
        self.evals_per_step = {"RK45": 6, "DOP853": 12}[method]

    def solve(self, problem: Problem, *, rtol: float = 1e-3, atol=1e-6) -> Solution:
        rhs = problem.rhs
        t0 = time.perf_counter()
        sol = solve_ivp(
            lambda _, y: rhs(y),
            problem.lam_span,
            problem.y0,
            method=self.method,
            rtol=rtol,
            atol=atol,
            events=[_event_function(e) for e in problem.events] or None,
        )
        wall = time.perf_counter() - t0
        if sol.status == -1:
            status = f"failed: {sol.message}"
        elif sol.status == 1:
            status = next(
                e.name
                for e, t in zip(problem.events, sol.t_events, strict=True)
                if e.stop_after and len(t) >= e.stop_after
            )
        else:
            status = "completed"
        dim = problem.y0.size
        t_events = sol.t_events or [[] for _ in problem.events]
        y_events = sol.y_events or [[] for _ in problem.events]
        events = {
            e.name: EventRecord(np.asarray(t), np.asarray(ys).reshape(-1, dim))
            for e, t, ys in zip(problem.events, t_events, y_events, strict=True)
        }
        return Solution(
            method=self.name,
            lam=sol.t,
            y=sol.y.T.copy(),
            status=status,
            nfev=int(sol.nfev),
            n_steps=int(sol.t.size - 1),
            wall_time=wall,
            events=events,
            options={"rtol": rtol, "atol": atol},
        )
