# 0003. Number of vector-field evaluations (nfev) as the primary cost metric

- Status: Accepted
- Date: 2026-10-06

## Context

RQ3 asks which method reaches a given accuracy most cheaply. Wall-clock time depends on language, implementation quality, interpreter overhead and hardware. Any comparison between SciPy (Python stepping loop) and compiled code in time would mostly measure the Python overhead ([pitfall 14](../pitfalls.md)). The numerical-ODE literature (Hairer–Nørsett–Wanner; SciMLBenchmarks work-precision diagrams) measures cost in evaluations of the vector field. That count is hardware-independent and reproducible exactly, and it is the dominant cost once the right-hand side is expensive, as it will be for Kerr and for generic $`\Gamma`$ contractions.

## Decision

**Primary cost = nfev**, with "one evaluation" defined as follows.

| Method family | One evaluation is | Counted |
|---|---|---|
| Explicit RK (RK4, DP5, DOP853) | one call of the formulation's `rhs(y)` | every stage of every attempted step, including rejected steps, the initial step-size selection, and the extra stages DOP853 needs for dense output |
| Implicit RK (Gauss–Legendre) | one call of `rhs(y)` | every stage in every fixed-point iteration, including the final iteration that only confirms convergence; Newton Jacobians are counted separately as `njev` |
| Tao (extended phase space) | one evaluation of the gradient pair $`(\partial H/\partial x, \partial H/\partial p)`$ at one phase-space point, the same work as one `rhs` of formulation (b) | each $`\phi_A`$ and $`\phi_B`$ sub-flow; consecutive $`\phi_A`$ flows merged across steps count once; the coupling rotation $`\phi_C`$ is free |
| SciPy `solve_ivp` | as reported by SciPy (`sol.nfev`) | our own DP5 and DOP853 must reproduce it within ±1 % (Phase 2 criterion) |

Recorded with every run, next to nfev: `n_steps`, `n_rejected`, `n_iter` (implicit methods), `njev`, and `wall_time`.

**Wall time** is compared only between our own Numba implementations. It is taken after JIT warm-up, as the median of at least five repeats, with `NUMBA_NUM_THREADS=1` and `OMP_NUM_THREADS=1`. CPU model, library versions and the git hash are stored with the result.

**Across formulations, nfev is not a fair currency.** An evaluation of (a), a $`\Gamma`$ contraction, does not cost the same as an evaluation of (b), a $`\partial g^{\mu\nu}`$ contraction. Cost comparisons between formulations therefore use wall time of our own implementations.

## Consequences

Positive:

- Work-precision diagrams (F6) and the cost table T5 are reproducible on any machine and directly comparable with the literature.
- SciPy runs take part in cost comparisons through nfev, even though their wall time is not meaningful.

Negative, and what we do about it:

- **nfev ignores other work**: the Newton linear solves of implicit methods, the stage-combination arithmetic, and memory traffic. For the cheap 8D Schwarzschild right-hand side these are not negligible. Table T5 therefore reports nfev **and** wall time, and the discussion addresses cases where the two rank methods differ.
- **Implicit-method nfev depends on the iteration stopping rule.** The rule (iterate to round-off stagnation) is fixed in Phase 3 and its effect is itself an experiment (F15).
- **Tao's extended state has 16 components**, but one gradient evaluation is charged as one `rhs` of (b). This is fair in floating-point work, since each sub-flow evaluates $`\nabla H`$ once at one 8D point. The doubled state shows up in wall time, which is reported alongside.

## Alternatives considered

- **Wall time only.** Dominated by implementation and interpreter overhead, not by the algorithm.
- **Number of steps.** Unfair between methods with different numbers of stages, and undefined for the iterations inside implicit methods.
- **Flop counting.** Precise in principle but impractical to instrument and to keep consistent across implementations.
