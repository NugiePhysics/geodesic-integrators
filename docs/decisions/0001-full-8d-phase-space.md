# 0001. Full 8D phase space and physics-agnostic integrators

- Status: Accepted
- Date: 2026-10-06

## Context

The predecessor project, `schwarzschild-geodesics`, integrates the reduced equatorial state $`(t, r, \phi, p_r)`$ and passes $`E`$ and $`L`$ in as fixed parameters. That was the right call for a single RK4 study, but it blocks this project in three ways:

1. **Conservation of $`E`$ and $`L`$ cannot be measured.** They are inputs, not outputs. Research questions RQ1 and RQ2 ask exactly how integrators treat these invariants.
2. **The metric interface only knows $`f(r)`$, $`f'(r)`$ and the horizon.** It cannot express Kerr, which has $`g^{t\phi} \neq 0`$ and genuine $`\theta`$-motion governed by the Carter constant. Adding Kerr without touching the integrators is an explicit acceptance test of the architecture (roadmap §4).
3. **Inclined orbits (TC6) do not exist** in an equatorial state.

Integrators must also be verifiable on toy problems (harmonic oscillator, Kepler, Tao's non-separable example) before any physics is involved. That rules out integrators that know about metrics.

## Decision

1. **State vector.** Formulation (a) integrates $`y = (x^\mu, u^\mu)`$ and formulation (b) integrates $`y = (x^\mu, p_\mu)`$. Both are contiguous `float64` arrays of shape `(8,)` ordered $`(t, r, \theta, \phi)`$ followed by the velocity or momentum in the same order. Tao's extended state is `(16,)`, ordered $`(q, p, x, y)`$.
2. **Physics-agnostic integrators.** An integrator sees only a vector field `rhs(y)`, or, for Hamiltonian splitting methods, the pair `dH_dx(x, p)`, `dH_dp(x, p)`, plus the dimension. Nothing in `geoint.integrators` imports from `geoint.metrics` or `geoint.formulations`. The dependency direction is `metrics → formulations → integrators → testcases → experiments → plotting`, with integrators importing no physics.
3. **General metric interface.** `Metric` provides $`g_{\mu\nu}`$, $`g^{\mu\nu}`$, $`\partial_\alpha g^{\mu\nu}`$ (shape 4×4×4), $`\Gamma^\mu_{\alpha\beta}`$ (shape 4×4×4), the horizon, stopping surfaces and the indices of cyclic coordinates. **No code may assume a diagonal metric.**
4. **No equatorial special case.** Equatorial orbits are 8D runs with $`\theta = \pi/2`$, $`p_\theta = 0`$.

## Consequences

Positive:

- $`E`$, $`L_z`$, $`L^2`$ and $`H`$ (or the norm $`g_{\mu\nu}u^\mu u^\nu`$) are measurable in both formulations ([theory §5](../theory/formulations.md#5-conserved-quantities)).
- Kerr becomes `metrics/kerr.py` plus analytic references and test cases, with zero lines changed in `integrators/` and `formulations/`.
- Every integrator is verified on TC0 independently of general relativity.

Negative, and what we do about it:

- **Twice the state size** of the reduced system, so roughly twice the arithmetic per step. Acceptable, because costs are only compared between methods on the same formulation ([ADR 0003](0003-nfev-as-cost-metric.md)).
- **Trivial invariants in (b).** $`p_t`$ and $`p_\phi`$ are conserved bit for bit by every Runge–Kutta method ([theory §5.2](../theory/formulations.md#5-conserved-quantities)). The conservation experiments must show $`E`$, $`L_z`$ in formulation (a) and $`H`$, $`L^2`$ in (b), and must not present the flat (b) curves as a result.
- **Unbounded cyclic coordinates.** $`t`$ and $`\phi`$ grow without bound and enter the error norm of adaptive methods, which silently loosens `rtol` (roadmap pitfall 5). This calls for a per-component `atol`, or for excluding $`t`$ from the norm.
- **Rounding-level tilt.** `cos(pi/2) ≈ 6e-17` gives equatorial orbits a $`p_\theta`$ of order $`10^{-16}`$ (pitfall 3). This is documented and tested as expected behaviour.
- **Polar singularity.** The singularity at $`\sin\theta = 0`$ becomes reachable, so inclinations are capped at 85° (pitfall 2).

## Alternatives considered

- **Reduced 4D state with $`E`$, $`L`$ as parameters** (the predecessor's design). Rejected for points 1–3 of the context.
- **6D state without $`t`$** (it never feeds back into the dynamics). Saves little, and loses coordinate time, which the photon-sphere Lyapunov exponent and $`\Omega = d\phi/dt`$ are measured in.
- **Cartesian Kerr–Schild coordinates.** No polar singularity and horizon-penetrating, but every analytic reference is in Schwarzschild/Boyer–Lindquist coordinates. Deferred to an extension (roadmap §8.3).
