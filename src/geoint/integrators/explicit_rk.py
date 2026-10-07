"""Explicit Runge-Kutta methods with a fixed step, from a Butcher tableau (classical RK4)."""

from __future__ import annotations

from ._compiled import EXPLICIT, make_params
from .base import Integrator, Problem, Solution, fixed_step_solution, one_step, run_fixed_step
from .tableaus import RK4, Tableau


class ExplicitRK(Integrator):
    """Explicit Runge-Kutta method with constant step size (``lam = lam0 + n h``)."""

    def __init__(self, tableau: Tableau = RK4, name: str = "RK4"):
        self.tableau = tableau
        self.name = name
        self.order = tableau.order
        self.evals_per_step = tableau.stages
        self._params = make_params(A=tableau.A, b=tableau.b)

    def solve(
        self,
        problem: Problem,
        *,
        h: float | None = None,
        n_steps: int | None = None,
        save_every: int = 1,
        compensated: bool = True,
    ) -> Solution:
        out, h, wall = run_fixed_step(
            EXPLICIT, self._params, problem, h, n_steps, save_every, compensated
        )
        return fixed_step_solution(self.name, problem, out, h, wall, {"compensated": compensated})

    def step(self, rhs, y, h: float):
        return one_step(EXPLICIT, self._params, rhs, y, h)
