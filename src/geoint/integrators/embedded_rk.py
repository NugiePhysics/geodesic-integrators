"""Adaptive Dormand-Prince pairs DP5(4) and DOP853, compiled.

The step-size controller, the initial step selection and the error norms reproduce
``scipy.integrate.RK45`` and ``DOP853`` line for line (SciPy 1.18, ``_ivp/rk.py`` and
``_ivp/common.py``), so that the two implementations take the same steps and SciPy can serve
as the oracle. Only the summation order inside the dot products differs, which changes
results at the rounding level. Events are located by partial steps of the method itself.
"""

from __future__ import annotations

import time

import numpy as np

from ._compiled import adaptive, make_params
from .base import Integrator, Problem, Solution, event_records, status_name
from .tableaus import DOP853, DP5, EmbeddedTableau


class EmbeddedRK(Integrator):
    """Adaptive explicit pair with SciPy's controller; ``atol`` may be a per-component array."""

    is_adaptive = True

    def __init__(self, tableau: EmbeddedTableau):
        self.tableau = tableau
        self.name = tableau.name
        self.order = tableau.order
        self.evals_per_step = tableau.stages

    def solve(
        self,
        problem: Problem,
        *,
        rtol: float = 1e-3,
        atol=1e-6,
        first_step: float | None = None,
        max_step: float = np.inf,
        save_every: int = 1,
        max_steps: int = 10**8,
    ) -> Solution:
        t = self.tableau
        lam0, lam1 = problem.lam_span
        atol_arr = np.broadcast_to(np.asarray(atol, dtype=float), problem.y0.shape).copy()
        t0 = time.perf_counter()
        out = adaptive(
            problem.rhs, problem.y0, lam0, lam1, float(rtol), atol_arr,
            make_params(A=t.A, b=t.B, f_size=problem.y0.size), t.E, t.E3,
            t.name == "DOP853", t.error_estimator_order, first_step or 0.0, float(max_step),
            *problem.event_arrays(), max(1, int(save_every)), int(max_steps),
        )  # fmt: skip
        wall = time.perf_counter() - t0
        lam, y, status, nfev, n_steps, n_rej, ev_ids, ev_lams, ev_ys = out
        return Solution(
            method=self.name,
            lam=lam,
            y=y,
            status=status_name(status, problem.events),
            nfev=int(nfev),
            n_steps=int(n_steps),
            n_rejected=int(n_rej),
            wall_time=wall,
            events=event_records(problem, ev_ids, ev_lams, ev_ys),
            options={"rtol": rtol, "atol": atol},
        )


DP5_INTEGRATOR = EmbeddedRK(DP5)
DOP853_INTEGRATOR = EmbeddedRK(DOP853)
