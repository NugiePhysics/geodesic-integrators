# Development Log

One entry per working session or phase: what was done, what was learned, and what went wrong and why. Newest entries go at the top.

---

## 2026-10-06 · Phase 0: setup and theory

### Done

- Created the repository with tooling adapted from `schwarzschild-geodesics`: hatchling packaging, ruff, pre-commit, GitHub Actions CI, dev container, MIT license, `CITATION.cff`. The package is `geoint`; the distribution is `geodesic-integrators`.
- Locked the environment with uv (`uv.lock`). Runtime dependencies are numpy, scipy, matplotlib, numba and mpmath; the `dev` dependency group adds pytest, hypothesis, sympy, pandas, ruff, pre-commit and nbmake; JAX is an optional extra. CI installs with `uv sync --locked`, so an out-of-date lock file fails the build.
- Supported Python is 3.12–3.14. The floor comes from NumPy 2.5, SciPy 1.18 and JAX, not from Numba (0.68 supports 3.10–3.15).
- Smoke test (`tests/test_smoke.py`): imports, plus the Numba pattern the design depends on, a jitted loop that receives a jitted RHS closure as an argument. It passes with and without `NUMBA_DISABLE_JIT`.
- Wrote [`theory/formulations.md`](theory/formulations.md): both formulations, the nine Christoffel symbols, $`\partial_\alpha g^{\mu\nu}`$, Hamilton's equations, the invariants $`E`$, $`L_z`$, $`L^2`$, $`H`$ in both sets of variables, and why $`p_t`$, $`p_\phi`$ are exact in formulation (b).
- Recorded ADRs [0001](decisions/0001-full-8d-phase-space.md) (8D state, physics-agnostic integrators), [0002](decisions/0002-numba-first.md) (Numba-first) and [0003](decisions/0003-nfev-as-cost-metric.md) (nfev as cost).

### Learned

- *The RQ2 hypothesis needs one qualification.* In formulation (b), $`p_t`$ and $`p_\phi`$ are conserved bit for bit by every Runge–Kutta method, explicit or implicit, adaptive or not, because their increments are exact floating-point zeros. **Tao's method is the exception.** Its coupling flow rotates $`(q - x, p - y)`$, and the cyclic positions of the two copies drift apart ($`\dot q^t = g^{tt}(x)p_t`$ vs $`\dot x^t = g^{tt}(q)y_t`$). Rounding aside, only $`p_t + y_t`$ is conserved. A throwaway NumPy Tao-2 run confirmed this, on the orbit $`r_\text{apo} = 20`$, $`L = 4.2`$ with $`h = 0.5`$, $`\omega = 20`$ and $`2\times10^4`$ steps. Coupling all eight components gave $`\max|\Delta E/E| = 1.7\times10^{-6}`$ for one copy and exactly 0 for the copy average. Coupling only $`(r, \theta)`$ gave exactly 0 for both, while the decoupled $`q^t - x^t`$ grew to $`2\times10^{-3}`$ with no effect on the dynamics. Phase 3 will turn this into a real test and decide which coupling to use. Details in theory §5.2.
- None of the invariants is quadratic in the state, because all have $`x`$-dependent coefficients. So the "Gauss methods preserve quadratic invariants" theorem gives nothing here. Symplectic methods can only offer bounded, not exact, conservation of $`H`$ and $`L^2`$.
- Captured closure values are compile-time constants in Numba. Per-run parameters (step size, tolerances, $`\omega`$) must be passed as arguments, or every sweep point recompiles (ADR 0002).

### Open

- CI has only run locally (Python 3.12 and 3.14, lint and tests). The exit criterion "CI green" needs the GitHub repository to be created and pushed.
- Core reading list ([`reading-list.md`](reading-list.md)): not started.
