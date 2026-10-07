"""Tao's explicit symplectic method for non-separable Hamiltonians (Tao 2016, PRE 94, 043303).

The phase space is doubled to ``z = (q, p, x, y)`` with

    H~ = H(q, y) + H(x, p) + omega * (|q - x|^2 + |p - y|^2) / 2.

Each of the three parts has an exact, explicit flow:

- ``phi_A`` (``H(q, y)``): ``p -= delta dH/dq(q, y)``, ``x += delta dH/dp(q, y)``;
- ``phi_B`` (``H(x, p)``): ``q += delta dH/dp(x, p)``, ``y -= delta dH/dq(x, p)``;
- ``phi_C`` (coupling): a rotation by ``2 omega delta`` of ``(q - x, p - y)``, sums unchanged.

Both gradients of ``phi_A`` are taken at the same point, so one flow costs exactly one call of
the canonical ``rhs(q, p) = (dH/dp, -dH/dq)``, the same work as one evaluation of formulation
(b) (ADR 0003). The order-2 step is ``A(h/2) B(h/2) C(h) B(h/2) A(h/2)``; the order-4 step is
its triple jump (Yoshida 1990). Consecutive ``phi_A`` flows at an unchanged ``(q, y)`` reuse
the cached gradient, so a step costs 3 (order 2) or 9 (order 4) evaluations.

``coupling`` selects the coordinates whose copies are coupled by ``phi_C``; uncoupled
coordinates evolve under ``phi_A`` and ``phi_B`` only. The reported solution is the first copy
``(q, p)``; the full extended state is kept for the diagnostics ``|q - x|``, ``|p - y|``.
"""

from __future__ import annotations

import numpy as np

from ._compiled import GAMMA1, GAMMA2, TAO, make_params
from .base import Integrator, Problem, Solution, fixed_step_solution, one_step, run_fixed_step

STRANG = np.array([1.0])
TRIPLE_JUMP = np.array([GAMMA1, GAMMA2, GAMMA1])


def extend(y0: np.ndarray) -> np.ndarray:
    """Extended initial state ``(q, p, q, p)``."""
    n = y0.size // 2
    return np.concatenate([y0, y0[:n], y0[n:]])


class Tao(Integrator):
    """Tao's method of order 2 (Strang) or 4 (triple jump) with constant step size."""

    is_symplectic = True
    is_symmetric = True
    requires_hamiltonian = True

    def __init__(self, order: int):
        if order not in (2, 4):
            raise ValueError("Tao's method is implemented for orders 2 and 4")
        self.order = order
        self.name = f"Tao{order}"
        self.fractions = STRANG if order == 2 else TRIPLE_JUMP
        self.evals_per_step = 3 * self.fractions.size

    def solve(
        self,
        problem: Problem,
        *,
        omega: float,
        h: float | None = None,
        n_steps: int | None = None,
        coupling=None,
        save_every: int = 1,
        compensated: bool = True,
    ) -> Solution:
        if not problem.hamiltonian:
            raise ValueError("Tao's method needs a canonical Hamiltonian vector field")
        dim = problem.y0.size
        n = dim // 2
        mask = np.ones(n, dtype=np.bool_) if coupling is None else np.asarray(coupling, bool)
        if mask.shape != (n,):
            raise ValueError(f"coupling must have shape ({n},)")
        params = make_params(
            fractions=self.fractions, coupling=mask, scalars=[0.0, 0.0, omega], cache_size=dim
        )
        extended = Problem(problem.rhs, extend(problem.y0), problem.lam_span, problem.events, True)
        out, h, wall = run_fixed_step(TAO, params, extended, h, n_steps, save_every, compensated)
        options = {"omega": omega, "coupling": mask.tolist(), "compensated": compensated}
        return fixed_step_solution(self.name, extended, out, h, wall, options, extended_dim=dim)

    def step(self, rhs, z, h: float, omega: float = 20.0):
        """One step on the *extended* state ``z = (q, p, x, y)``."""
        n = np.size(z) // 4
        params = make_params(
            fractions=self.fractions, coupling=np.ones(n, np.bool_), scalars=[0.0, 0.0, omega],
            cache_size=2 * n,
        )  # fmt: skip
        return one_step(TAO, params, rhs, z, h)


TAO2 = Tao(2)
TAO4 = Tao(4)
