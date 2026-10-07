# Extending to Kerr

The architecture has one acceptance test: **adding a new spacetime must not change a single line in `geoint/integrators/` or `geoint/formulations/`**. Integrators see only a vector field ([ADR 0001](decisions/0001-full-8d-phase-space.md)), and the formulations contract the full $`g^{\mu\nu}`$ and $`\partial_\alpha g^{\mu\nu}`$ without assuming a diagonal metric ([ADR 0004](decisions/0004-metric-interface.md)). This page lists what Kerr in Boyer–Lindquist coordinates needs, and where.

## 1. The metric (`geoint/metrics/kerr.py`)

Subclass `geoint.metrics.Metric` and provide five Numba kernels of a position `x = (t, r, theta, phi)`:

| Kernel | Shape | Entry |
|---|---|---|
| `g(x)` | (4, 4) | $`g_{ij}`$ |
| `g_inv(x)` | (4, 4) | $`g^{ij}`$, including $`g^{t\phi} \neq 0`$ |
| `dg_inv(x)` | (4, 4, 4) | $`\partial_i g^{jk}`$, derivative index first |
| `christoffel(x)` | (4, 4, 4) | $`\Gamma^i_{jk}`$, upper index first |
| `invariants(x, p)` | (n,) | constants of motion in canonical variables |

and the metadata:

- `invariant_names = ("E", "Lz", "Q")`, with Carter's constant $`Q = p_\theta^2 + \cos^2\theta\left(a^2(\mu^2 - E^2) + L_z^2/\sin^2\theta\right)`$ ($`\mu^2 = \epsilon`$). For $`a = 0`$ it reduces to $`L^2 - L_z^2`$, so the TC6 code path that monitors $`L^2`$ already exercises it.
- `cyclic_coords = (0, 3)`: the slices `dg_inv(x)[0]` and `dg_inv(x)[3]` must be exact zeros, so that $`p_t`$ and $`p_\phi`$ stay conserved bit for bit in formulation (b) (theory §5.2).
- `horizon`: the outer horizon $`r_+ = M + \sqrt{M^2 - a^2}`$.

`tools/derive_metric.py` derives $`g^{\mu\nu}`$, $`\partial g^{\mu\nu}`$ and $`\Gamma`$ from $`g_{\mu\nu}`$ with SymPy. Use it to check hand-written kernels at 50 digits, as `tests/unit/test_schwarzschild_metric.py` does for Schwarzschild, or as a starting point for code generation. Write the kernels in terms that stay accurate near the horizon ($`\Delta = r^2 - 2Mr + a^2`$ factored as $`(r - r_+)(r - r_-)`$), as Schwarzschild uses $`r - 2M`$.

## 2. Initial data

`geoint.initial_conditions.solve_quadratic_shell` already handles off-diagonal metrics: the mass shell is a quadratic in one momentum component with a linear term, $`A p_k^2 + 2B p_k + C = 0`$, and the root is chosen by the sign of the velocity component (theory §7). Only builders for Kerr's natural parameters (energy, $`L_z`$, $`Q`$) need adding.

## 3. Analytic references (`geoint/analytic/kerr*.py`)

- Prograde and retrograde ISCO and photon orbits (Bardeen, Press & Teukolsky 1972).
- Spherical photon orbits and the shadow boundary.
- Fundamental frequencies and orbits in closed form (Fujita & Hikida 2009, *Class. Quantum Grav.* 26, 135002) and null geodesics (Gralla & Lupsasca 2020, *Phys. Rev. D* 101, 044032).

As for Schwarzschild, every reference should have two independent evaluation paths ([pitfall 17](pitfalls.md)).

## 4. Test cases (`geoint/testcases/kerr.py`)

Subclass `geoint.testcases.base.GeodesicCase` as frozen dataclasses: `initial_xp`, `lam_end`, `events`, `reference` and `errors`. Natural first cases: equatorial circular orbits (prograde and retrograde), inclined orbits with $`\delta Q`$ as the error, and spherical photon orbits as the analogue of TC2.

## 5. The experiment harness

`geoint.experiments.runner.formulation` builds both formulations for Schwarzschild. It will take the metric from the run specification, which is the only change outside the new files. The cache key already contains the test case's family and parameters.

## Checklist

- [ ] Kernels agree with the SymPy derivation to about $`10^{-15}`$ at random points.
- [ ] `cyclic_coords` slices of `dg_inv` are exact zeros; $`p_t`$, $`p_\phi`$ conserved bit for bit in (b).
- [ ] Formulations (a) and (b) give the same trajectories to $`10^{-9}`$ with DOP853 at $`10^{-13}`$.
- [ ] Every new test case agrees with its reference to $`10^{-10}`$ with tight DOP853.
- [ ] `git diff --stat` shows no change under `src/geoint/integrators/` or `src/geoint/formulations/`.
