# Development Log

One entry per working session or phase: what was done, what was learned, and what went wrong and why. Newest entries go at the top.

---

## 2026-10-07 · Phase 1: metric and formulations

### Done

- **Metric interface** (`geoint.metrics`). A metric is five Numba kernels returning full arrays: `g`, `g_inv`, `dg_inv` ($`\partial_i g^{jk}`$, derivative index first), `christoffel` ($`\Gamma^i_{jk}`$, upper index first) and `invariants(x, p)`. It also declares its horizon, cyclic coordinates and stopping surfaces. Recorded in [ADR 0004](decisions/0004-metric-interface.md).
- **`Schwarzschild`**, hand-written in terms of $`r - 2M`$ rather than $`f = 1 - 2M/r`$, so components stay accurate to full relative precision right down to the horizon. Kernels are compiled once per mass and shared between instances.
- **`tools/derive_metric.py`** derives $`g^{\mu\nu}`$, $`\partial g`$, $`\partial g^{-1}`$ and $`\Gamma`$ from $`g_{\mu\nu}`$ with SymPy. `uv run python tools/derive_metric.py schwarzschild` prints the 9 Christoffel symbols and 5 entries of $`\partial_\alpha g^{\mu\nu}`$ of the theory tables. The tests evaluate it at 50 digits.
- **Formulations** `SecondOrder` (a) and `Hamiltonian` (b). Both have a compiled `rhs(y)`, `from_xp`/`to_xp`, and `invariants(y)` (H, E, Lz, L2) for a single state or a whole trajectory, plus `constraint` and `relative_constraint` (theory §5.1). `Hamiltonian` also exposes `dH_dx` and `dH_dp` for Tao. Both contract over every index pair; neither assumes a diagonal metric.
- **Initial data** (`geoint.initial_conditions`). A generic solver for the mass shell as a quadratic in one momentum component (theory eq. 25), ready for off-diagonal metrics. The root is selected by the sign of the velocity component and evaluated in its cancellation-free form. Builders: `momentum_from_constants` (solve for $`p_r`$) and `turning_point_momentum` ($`p_r = 0`$, solve for $`E`$).
- 128 tests, which pass with and without `NUMBA_DISABLE_JIT` and on Python 3.12 and 3.14. Beyond the exit criteria they cover: exact symmetry of $`g`$, $`\partial g^{-1}`$ and $`\Gamma`$; exact zeros for cyclic coordinates and for $`\dot p_t`$, $`\dot p_\phi`$; reversibility $`\rho F(y) = -F(\rho y)`$ bit for bit (theory eq. 23); the (a)–(b) identity of theory eq. (17) at random off-shell points; agreement of the generic $`\Gamma`$ contraction with the hand-expanded eq. (8); and a hypothesis property test of the initial data.

### Exit criteria

| Criterion | Target | Measured |
|---|---|---|
| $`\Gamma`$ and $`\partial g^{\mu\nu}`$ vs SymPy at random points | $`\le 10^{-14}`$ | $`\le 6\times10^{-16}`$ relative per component, for $`M = 1`$ and $`2.5`$, with $`r - 2M`$ from $`10^{-8}M`$ to $`10^4 M`$ |
| (a) vs (b), $`r(\lambda)`$, $`\theta(\lambda)`$, $`\phi(\lambda)`$ with DOP853 at rtol $`10^{-13}`$ | $`\le 10^{-9}`$ | $`\le 2.8\times10^{-10}`$ (eccentric orbit $`r_\text{apo} = 20`$, $`L = 4.2`$, equatorial and inclined by 60°; photon with $`b = 6`$ from $`r = 1000`$) |
| Initial relative constraint | $`\le 10^{-15}`$ | $`\le 3.3\times10^{-16}`$ on all 39 orbits of the test-case matrix, in both formulations |

### Learned

- **The turning-point problem is about $`p_r`$, not about the constraint.** Near a turning point $`D = E^2 - V_\text{eff}`$ has an absolute error of a few ulp, so $`p_r`$ comes out of order $`\sqrt{\epsilon_\text{mach}}`$. For the circular orbit at $`r = 7`$, started from the closed-form $`E`$ and $`L`$, it was $`1.8\times10^{-8}`$ instead of 0. The constraint is still satisfied to rounding, because the error in $`p_r^2`$ is the rounding error of $`D`$. For a circular orbit, though, this is a spurious radial oscillation that TC3 would read as integrator error. The fix is to start at the turning point and solve for $`E`$, which involves no cancellation, rather than for $`p_r`$. The theory doc (§7) now says so. The roadmap's idea of factoring the radial polynomial is only needed for analytic references.
- **A reference must be more accurate than the thing it tests.** The first version of the eq. (17) test evaluated the SymPy $`\partial g_{\mu\nu}`$ in float64. Near the horizon its $`1 - 2M/r`$ form loses digits, and the test failed by a relative $`7\times10^{-12}`$ at $`r - 2M = 2\times10^{-5}`$. The kernels were right; the reference was wrong. All SymPy references are now evaluated in mpmath at 30–50 digits.
- **"Relative constraint $`\le 10^{-15}`$" only holds where the terms of $`H`$ are $`O(1)`$.** Close to the horizon they are of size $`E^2/f`$, and $`H`$ cannot be evaluated more accurately than a few ulp of that. The random-orbit property test therefore bounds $`|C|`$ by 16 ulp of $`\tfrac12\sum|g^{\mu\nu}p_\mu p_\nu|`$, not by an absolute number. It holds for $`r - 2M`$ down to $`10^{-6}`$.
- **Near-critical photons preview RQ4.** At $`b = 5.1962`$ ($`b - b_c \approx 5\times10^{-5}`$), (a) and (b) end up $`4.5\times10^{-7}`$ apart after about 1.75 windings, against $`10^{-10}`$ for $`b = 6`$. That is the $`e^{2\pi}`$-per-winding amplification of the photon sphere acting on $`10^{-11}`$-level differences, not a defect of either formulation. Cross-formulation agreement is therefore tested on well-conditioned orbits only (pitfall 19).
- **Check shapes at the Python boundary.** A post-push review found that `invariants` reshaped any input to `(-1, 8)`. A 16-component vector, which is exactly Tao's extended state in Phase 3, was therefore read silently as two states, and `(2, 3, 8)` batches came back flattened. `from_xp` and `to_xp` did not check lengths either. All three now reject wrong shapes, batches keep their leading shape, and regression tests cover both cases.
- **Cost baseline** (i5-7300HQ, single thread, called from a compiled loop): 135 ns per `rhs` of (a), 331 ns per `rhs` of (b). Formulation (b) calls two kernels (`g_inv` and `dg_inv`), so it pays twice for allocation and trigonometry. Before any cross-formulation wall-time comparison (ADR 0003), this gap must be removed or explained.

### Open

- The JAX autodiff cross-check of roadmap §5, Level 1, is deferred. SymPy at 50 digits plus fourth-order finite differences already cover the hand-written derivatives. JAX pays off once metrics are prototyped quickly (Kerr–(A)dS).
- The off-diagonal code paths are exercised only by the shell solver's tests on random Lorentzian forms. The first real metric with $`g^{t\phi} \neq 0`$ will be Kerr.
- Tag `v0.1` once CI is green on GitHub.

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
