# Development Log

One entry per working session or phase: what was done, what was learned, and what went wrong and why. Newest entries go at the top.

---

## 2026-10-07 · Phase 5: short integrations, convergence and work-precision

### Done

- Experiment scripts in `experiments/`, one per figure or table; `make figures` (`scripts/make_figures.py`) runs them in dependency order. Integrations are cached under `results/raw/` (git-ignored); small summary tables go to `results/summary/` (committed), figures to `figures/` as PDF (report) and PNG (README).
- Shared tools: `geoint.experiments.sweeps` (grids, per-method options, Tao's ω), `geoint.experiments.analysis` (slope fits over the asymptotic window, Pareto fronts, cost at a target error, per-period envelopes, growth classes), `geoint.plotting.style` (one colour and marker per method from the validated reference palette, with a marker on every series because three slots are below 3:1 contrast on white).
- Figures F2–F6 and F12–F16, tables T2, T3 and T5.

### What each figure shows

- **F14, Tao against ω.** The error of Tao4 grows in proportion to ω once $`\omega h \gtrsim 10^{-3}`$, and the method diverges for $`\omega h \gtrsim 1.4`$. Too small an ω decouples the copies: with $`\omega = 10^{-5}`$ they separate by up to 23 within 300 periods. **ω = 10⁻³ is used everywhere** ([ADR 0006](decisions/0006-tao-coupling.md)): it is within 1.2–5× of the best ω in all four settings tested, and it holds the copies within $`3\times10^{-6}`$ over 300 periods.
- **F4 and T3, convergence.** Every fixed-step method reaches its theoretical order in both formulations: on the orbit $`(7.5, 0.5)`$ the measured orders are within ±0.05; on the deflection $`b = 6`$ all are within ±0.08. GL1 and Tao2 are almost indistinguishable in (b). GL1 is about 35 times less accurate in (a) than in (b) at the same step.
- **F5, tolerance proportionality.** The global error falls like $`\text{tol}^{0.85\text{–}1.04}`$ for DP5 and DOP853, and our curves lie on SciPy's.
- **F6 and T5, work-precision.** On these short integrations **DOP853 is the cheapest method everywhere**. To reach $`10^{-10}`$ in (b) it needs $`1.4`$–$`3.2\times10^3`$ evaluations, against $`1.6`$–$`3.6\times10^4`$ for GL3 (the best fixed-step method) and $`3.5\times10^4`$–$`1.7\times10^5`$ for RK4. Tao4 costs 1.2–10 times as much as RK4 for the same error. GL1 and Tao2 never reach $`10^{-10}`$ in the sweep. Formulation (b) is cheaper than (a) for most method–case pairs, but not all: on the inclined orbit DOP853 and GL3 reach $`10^{-10}`$ more cheaply in (a). *Caveat on the time axis:* one evaluation costs 0.7–1.6 µs inside our drivers, 2–5 times the bare vector field (Phase 1: 135/331 ns), because the generic step code allocates small arrays. Times compare our implementations with each other, as ADR 0003 intends; they are not a measure of the algorithms' intrinsic cost.
- **F3, deflection at a fixed budget** ($`2\times10^4`$ evaluations). In (b) the adaptive pairs reach $`10^{-11}`$–$`10^{-12}`$ for every $`b`$. In (a) they lose up to three orders near $`b_c`$ ($`1.5\times10^{-8}`$ at $`b = 5.3`$), the drift of $`L = r^2u^\phi`$ found in Phase 4. Errors of all methods rise towards $`b_c`$ with the conditioning.
- **F2, deflection angle.** The integration matches the exact $`\Delta\phi(r_\text{far})`$ to $`10^{-14}`$ relative far from the hole and to $`10^{-11}`$ at $`b - b_c = 10^{-3}`$, where the exact problem amplifies errors. The fourth-order weak-field series is good to $`10^{-11}`$ at $`b \approx 2000`$ and off by more than 1 % below $`b \approx 10`$. Bozza's limit is accurate near $`b_c`$ and fails above $`b - b_c \approx 1`$.
- **F12, windings near $`b_c`$.** The windings follow $`-\ln(b/b_c - 1)/2\pi`$, as Bozza predicts. The integration error grows exactly like the condition number $`1/(b - b_c)`$: $`3.3\times10^{-13}/(b - b_c)`$ windings in (b), $`1.4\times10^{-10}/(b - b_c)`$ in (a), 440 times worse because of the drifting L.
- **F13, periapsis advance.** The integration agrees with $`\Phi - 2\pi`$ to $`\le 4\times10^{-13}`$ away from the separatrix, and to $`1.5\times10^{-8}`$ at $`p = 7.02`$, where the advance is 21.3 rad (3.4 extra turns per orbit). The weak-field $`6\pi M/p`$ is off by 87 % there and still by 2 % at $`p = 207`$.
- **F15, the implicit iteration.** Iterations per step grow from 2–5 at small $`h`$ to 10–23 at the largest; GL3 needs the fewest, and (a) needs more than (b). Along an orbit, steps near periapsis need about one iteration more. Stopping the iteration at a tolerance makes GL2 drift: over 1000 periods $`|\delta H|`$ grows tenfold with tolerances $`10^{-6}`$ and $`10^{-8}`$, threefold with $`10^{-10}`$, and not at all with $`10^{-12}`$ or with iteration to stagnation. A tolerance of $`10^{-12}`$ stays bounded at 3.5 iterations per step, against 6.7 for stagnation: the round-off rule is safe but not the cheapest safe rule.
- **F16, adaptive step sizes.** On the deflection the step grows roughly in proportion to $`r`$. DOP853 needs three to four times fewer steps than DP5 at the same tolerance. The tiny first steps come from the start-up of the controller at a turning point.
- **T2, validation.** Every check agrees with its analytic value to $`10^{-11}`$ or better, except the cases amplified by the photon sphere (about $`10^{-8}`$).

### Learned

- F2 first plotted $`\Delta\phi(r_\text{far}) - \pi`$ against the *asymptotic* angle, and the points fell below the curve at large $`b`$. The finite-$`r_\text{far}`$ curve now sits next to it (pitfall 18 again, this time in a figure).
- F15's first version reported a zero drift for iteration to stagnation. The last per-period window held a single sample. Envelopes now drop a short trailing window.
- The wall-time runs of F6 were first made while the interpreted test suite ran on another core. They were deleted from the cache and re-timed on an idle machine.

### Open

- `make figures` has not yet been run from an empty cache in one go; every script has been run, but at different times.
- Phase 6 has started (long-term runs, photon sphere, ISCO, inclined orbits, round-off), but is not committed: the scripts are untracked work in progress.

---

## 2026-10-07 · Phase 4: analytic references and test-case registry

### Done

- `geoint.analytic`: light deflection through three independent paths (the elliptic integral via the roots of $`2u^3 - u^2 + b^{-2}`$, Darwin's $`(r_0, Q)`$ form, direct mpmath quadrature), the weak-field series, Bozza's strong-deflection limit, and the angle swept from a photon turning point. Timelike orbits by $`(p, e)`$: conversions to $`(E, L)`$ and back, $`\Phi`$ in closed form, $`r(\phi)`$ through Jacobi functions, and $`\Phi`$, $`T_\tau`$, $`T_t`$ by an independent 30-digit quadrature in the angle $`\chi`$ of $`r = \tfrac12(r_a + r_p) - \tfrac12(r_a - r_p)\cos\chi`$, which removes both turning-point singularities.
- `geoint.testcases`: the TC1–TC6 matrix as frozen dataclasses (deflection, near-critical photons, photon sphere, circular, eccentric and inclined orbits). Each case builds the initial data for either formulation and computes errors on physical quantities (Δφ at $`r_\text{far}`$, phase at a known proper time, in-plane angle, $`\delta(L^2)`$, orbital-plane tilt).

### Exit criteria

| Criterion | Target | Measured |
|---|---|---|
| Darwin vs quadrature, $`b \in [5.2, 10^4]`$ | $`\le 10^{-12}`$ | identical in float64 (all three paths), for $`r_\text{far} = \infty, 10^3, 10^4`$ |
| Weak-field and Bozza limits | met | series error $`< 10^3 b^{-5}`$; Bozza error shrinks 100× per decade of $`b - b_c`$ |
| Precession from $`r(\phi)`$ vs the formula for $`\Phi`$ | match | $`\Phi`$ closed form vs quadrature $`1.2\times10^{-16}`$; $`r(0) = r(\Phi) = r_p`$, $`r(\Phi/2) = r_a`$; orbit equation satisfied |
| All TCs with tight DOP853 vs reference | $`\le 10^{-10}`$ | (b): $`\le 2\times10^{-11}`$ on TC1, TC3, TC4, TC6. (a): same on TC3, TC4, TC6; TC1 up to $`1.2\times10^{-8}`$ (explained below). TC1b and TC2 are ill-conditioned (pitfall 19) |

### Learned

- **I first transcribed Darwin's formula wrongly**, as $`-\pi + 4\sqrt{r_0/Q}\,F(\zeta, k)`$. The correct form has $`K(k) - F(\zeta, k)`$ (Iyer & Petters 2007, eq. 6). The two-path test caught it at once, with a disagreement of order 1.
- **The predecessor's "analytic" precession was not analytic.** Its README lists 2.127093985 rad for $`r_\text{apo} = 20`$, $`L = 4.2`$. Its test evaluated the elliptic formula at `r_peri = trajectory.r.min()`, the smallest *sampled* radius, which misses the true periapsis. The exact value is 2.12709400813 rad (closed form and a 40-digit quadrature agree). The old *numerical* value, 2.127094008, was right. This is the interpolation problem the roadmap set out to remove.
- **Formulation (a) loses accuracy on rays that go far out.** With DOP853 at $`10^{-13}`$, $`b = 5.3`$, $`r_\text{far} = 10^4`$: error $`1.2\times10^{-8}`$ in (a), $`3\times10^{-12}`$ in (b). In (a), $`L = r^2u^\phi`$ is rebuilt from a component that falls like $`1/r^2`$. An absolute tolerance cannot hold its relative accuracy, so L drifts by about $`10^{-10}`$ while E stays within $`4\times10^{-14}`$. The impact parameter $`b = L/E`$ drifts with it, and $`|d\Delta\phi/db|\,\delta b`$ reproduces the error within a factor of 2. In (b), L is a state variable that never changes. This is a first, quantitative answer to RQ2 for ray tracing.
- **The equatorial photon sphere is an exact fixed point of every Runge–Kutta map.** At $`r = 3`$ with $`p_r = 0`$, the $`(r, p_r)`$ components of the vector field vanish and $`t`$, $`\phi`$ are linear in λ, so no truncation error can perturb the orbit. With tight tolerances it survived 900 orbits; what eventually tips it is $`\cos(\pi/2) \approx 6\times10^{-17}`$ tilting the plane (pitfall 3). RQ4 therefore has to be measured on an *inclined* photon-sphere orbit, where the truncation error in $`L^2`$ kicks $`r`$. TC2 got an `inclination` parameter. A first look: 3.0 → 4.66 orbits survived as the tolerance goes from $`10^{-8}`$ to $`10^{-13}`$, about 0.32 orbits per decade against the predicted $`\ln 10/2\pi = 0.37`$.
- **Linear theory misses the photon-sphere exit by about 0.0035 orbits**, from the nonlinear terms at $`|r - 3| = 0.1`$. The exact exit angle is an elliptic integral from the turning point. Because $`u_2 = 1/r_0`$ is known exactly there, the other two roots follow without solving the cubic, which matters since $`b - b_c = O((r_0 - 3)^2)`$ falls below float64 resolution for $`r_0 - 3 < 10^{-8}`$. The measured growth rate is 0.1925, against $`\lambda_L = 1/(3\sqrt3) = 0.19245`$.
- `np.unwrap` lost whole turns of the in-plane angle when an adaptive step swept more than π. For prograde orbits the in-plane angle and $`\phi`$ never differ by more than $`\pi/2`$, so whole turns are now taken from $`\phi`$, which is continuous in the state.
- Circular orbits (TC3) are, like the photon sphere, fixed points of the radial dynamics. Their errors ($`\le 10^{-11}`$) measure round-off and the initial data, not truncation error, so they are not used for convergence or work-precision studies. Inclined orbits are.

---

## 2026-10-07 · Phase 3: symplectic integrators

### Done

- Gauss–Legendre GL1–GL3 (closed-form coefficients). Fixed-point iteration on the stage increments $`Z`$, run to round-off stagnation. Start value extrapolated from the previous step's collocation polynomial. Update $`y + d^\mathsf{T}Z`$ with $`d = b^\mathsf{T}A^{-1}`$, added with compensated summation. Optional `iter_tol` for experiment F15.
- Tao's method of order 2 (Strang) and 4 (triple jump). Every flow is one call of the canonical `rhs`, and a one-entry gradient cache merges consecutive $`\phi_A`$ flows, giving 3 and 9 evaluations per step. Optional coupling mask. Compensated summation inside the flows.
- One-step maps (`Integrator.step`) for the geometric tests, and cumulative per-step costs at the saved points (`Solution.counts`).

### Exit criteria (TC0)

| Criterion | Target | Measured |
|---|---|---|
| Order GL1/GL2/GL3 | 2/4/6 ± 0.1 | 2.00 / 4.00 / 6.00 (Kepler, $`e = 0.5`$) |
| Order Tao2/Tao4 | 2/4 ± 0.1 | 2.00 / 3.98 |
| $`\Psi'^\mathsf{T}J\Psi' = J`$ | $`\le 10^{-10}`$ (Tao in extended space); RK4 clearly fails | GL and Tao $`\le 2\times10^{-12}`$, limited by the finite-difference Jacobian; RK4 $`4.3\times10^{-3}`$ |
| $`\Psi_{-h}\circ\Psi_h = \mathrm{id}`$ | met | GL2, GL3, Tao $`\approx 2\times10^{-16}`$; GL1 $`1.5\times10^{-14}`$; RK4 $`2.7\times10^{-4}`$ |
| $`\lvert\Delta H\rvert`$ bounded, $`10^6`$ steps | met | Tao's non-separable example: GL1, GL2, Tao2, Tao4 flat; RK4 doubles from the first to the second half (linear drift) |

### Learned

- **Hairer's stagnation rule needs a guard.** "Stop when $`\lVert\Delta Z\rVert`$ no longer decreases" assumes monotone convergence. On a rotation-like problem the fixed-point map has complex eigenvalues, the norm of the increments wobbles, and the iteration stopped at $`10^{-11}`$ instead of round-off: GL1 was symmetric only to $`2.5\times10^{-11}`$. Stagnation now counts only within 100 ulp of the stage increments.
- **The symplectic structure of Tao's extended space pairs $`(q, p)`$ and $`(x, y)`$.** My first test used the standard $`J`$ on $`(q, p, x, y)`$ and reported a defect of 1.7. With the block-diagonal $`J`$ the defect is $`10^{-11}`$.
- **Tao on geodesics: couple only the non-cyclic coordinates** ([ADR 0006](decisions/0006-tao-coupling.md)). With all coordinates coupled, the time copies drift apart by up to $`10^4`$, the rotation feeds that into $`p_t`$, and the copies run away. With $`(r, \theta)`$ only, $`E`$ and $`L_z`$ stay exact in both copies.
- **Tao's accuracy depends strongly on ω.** On $`(p, e) = (20, 0.5)`$ the error grows linearly in $`\omega h`$ above about $`10^{-3}`$, and the method diverges for $`\omega h \gtrsim 1.4`$. Even at its best, Tao4 is about 60 times less accurate than RK4 at the same step, while costing 9 evaluations per step against 4. ω is chosen from F14.
- Iterating to stagnation costs about 7–11 iterations per step: GL2 takes 13.5 evaluations per step on Tao's example.

---

## 2026-10-07 · Phase 2: explicit integrators and experiment harness

### Done

- Integrator interface (`Problem`, `Event`, `Solution`, `Integrator`) without any physics import. Fixed-step RK4. Adaptive DP5(4) and DOP853 with SciPy's controller copied line by line. A `solve_ivp` adapter as the oracle. Registry `geoint.integrators.get(name)`.
- Events located by root finding on partial steps of the method itself, accurate to the method's order. Fixed steps use λ = λ₀ + n·h.
- TC0 toy problems: harmonic oscillator, Kepler, Tao's non-separable $`H`$.
- Experiment harness `geoint.experiments`: `RunSpec` → `run` → one flat row with errors, costs and provenance (git hash, versions, CPU), cached as JSON under `results/raw/` by a hash of the specification; `run_many` forks worker processes.

### Exit criteria

| Criterion | Target | Measured |
|---|---|---|
| RK4 order | 4 ± 0.1 | 4.05 (Kepler); 4.00 on the strong-field orbit $`(7.5, 0.5)`$ in both formulations |
| Own DP5/DOP853 vs SciPy | nfev ± 1 %, solution ≤ tol | nfev and step counts *identical* at tolerances $`10^{-4}`$–$`10^{-12}`$ (Kepler and a geodesic), compiled and interpreted; step sizes agree to $`\sim10^{-7}`$; final states to $`10^{-12}`$ |
| TC0 | pass | oscillator, Kepler, Tao's example |
| Predecessor's validation table | reproduced | $`d\phi/dt`$ at $`r_c = 10`$, precession (see Phase 4), capture threshold, deflection at $`b = 40`$, RK4 order |

### Learned

- **Compilation dominated everything.** The first design compiled one driver per (method kind, vector field), 10–20 s each, which is minutes per process. Numba's disk cache refused the drivers, because passing a jitted function as a value embeds a dispatcher pointer. Typed first-class vector fields, a single step function that dispatches on an integer kind, and one module for all compiled code fixed it: about 60 s once per machine, about 1 s per process afterwards ([ADR 0005](decisions/0005-compiled-integrator-core.md)). The cold-cache test suite runs in 84 s.
- RK4 on $`(p, e) = (20, 0.5)`$ is pre-asymptotic over the whole usable range (slopes 4.5 → 4.07 before round-off), so the order is measured on the strong-field orbit $`(7.5, 0.5)`$, where it is 4.00.
- **Started at a turning point, adaptive step sequences are rounding noise.** The no-JIT test run caught our DOP853 parting from SciPy after two steps (nfev 2.5 % apart). The orbit started at $`p_r = 0`$, where the first error estimates are at round-off level. The step-growth factor $`0.9\,\text{err}^{-1/8}`$ then amplifies rounding noise, so our compiled code, our interpreted code and SciPy all take different second steps, and only the totals land close. From a generic starting point all three agree on nfev exactly and on step sizes to $`\sim10^{-7}`$. The oracle test now starts there.
- SciPy's `nfev` with events includes 3 extra evaluations per located event (DOP853 dense output), and ours includes 20–50 (partial steps). The `nfev` comparison is therefore made without events.

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
