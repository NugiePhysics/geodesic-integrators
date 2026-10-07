"""Gauss-Legendre collocation (implicit midpoint GL1, GL2, GL3): symplectic and symmetric.

Implementation after GNI VIII.5-VIII.6 and Hairer, McLachlan & Razakarivony (2008):

- the unknowns are the stage increments ``Z_i = Y_i - y``, solved by fixed-point iteration
  ``Z <- h A F(y + Z)``;
- by default the iteration runs to round-off *stagnation*: it stops when the increment
  ``||dZ||`` is zero, or no longer decreases while within 100 ulp of the stage increments
  (convergence of the fixed-point map need not be monotone); stopping earlier leaves an
  iteration error
  that is not symplectic and shows up as a linear drift of the energy (roadmap pitfall 10);
  ``iter_tol > 0`` instead stops at ``||dZ|| <= iter_tol``, for experiment F15;
- the starting guess extrapolates the collocation polynomial of the previous step;
- the update ``y + d^T Z`` with ``d = b^T A^{-1}`` needs no extra evaluation and is added with
  compensated summation.

The norm is ``max_k |dZ_k| / (1 + |y_k|)``.
"""

from __future__ import annotations

import numpy as np

from ._compiled import GAUSS, make_params
from .base import Integrator, Problem, Solution, fixed_step_solution, one_step, run_fixed_step
from .tableaus import GL1, GL2, GL3, GaussTableau


class GaussLegendre(Integrator):
    """``s``-stage Gauss-Legendre method of order ``2s`` with constant step size."""

    is_implicit = True
    is_symplectic = True
    is_symmetric = True

    def __init__(self, tableau: GaussTableau, name: str):
        self.tableau = tableau
        self.name = name
        self.order = tableau.order
        self.evals_per_step = tableau.stages

    def solve(
        self,
        problem: Problem,
        *,
        h: float | None = None,
        n_steps: int | None = None,
        save_every: int = 1,
        compensated: bool = True,
        iter_tol: float = 0.0,
        max_iter: int = 60,
    ) -> Solution:
        params = self._params(problem.y0.size, iter_tol, max_iter)
        out, h, wall = run_fixed_step(GAUSS, params, problem, h, n_steps, save_every, compensated)
        options = {"compensated": compensated, "iter_tol": iter_tol}
        return fixed_step_solution(self.name, problem, out, h, wall, options)

    def _params(self, dim: int, iter_tol: float = 0.0, max_iter: int = 60):
        t = self.tableau
        memory = np.zeros((t.stages + 1, dim))
        memory[t.stages, 0] = np.nan  # no previous step: no extrapolated starting guess
        return make_params(
            A=t.A, b=t.d, A_inv=t.A_inv, X=t.extrapolation, memory=memory,
            scalars=[iter_tol, max_iter, 0.0],
        )  # fmt: skip

    def step(self, rhs, y, h: float, iter_tol: float = 0.0):
        return one_step(GAUSS, self._params(np.size(y), iter_tol), rhs, y, h)


GL1_INTEGRATOR = GaussLegendre(GL1, "GL1")
GL2_INTEGRATOR = GaussLegendre(GL2, "GL2")
GL3_INTEGRATOR = GaussLegendre(GL3, "GL3")
