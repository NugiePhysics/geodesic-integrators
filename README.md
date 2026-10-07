# geodesic-integrators

[![CI](https://github.com/NugiePhysics/geodesic-integrators/actions/workflows/ci.yml/badge.svg)](https://github.com/NugiePhysics/geodesic-integrators/actions/workflows/ci.yml)
[![Docs](https://github.com/NugiePhysics/geodesic-integrators/actions/workflows/docs.yml/badge.svg)](https://nugiephysics.github.io/geodesic-integrators/)
[![Python](https://img.shields.io/badge/python-3.12%2B-blue)](pyproject.toml)
[![License: MIT](https://img.shields.io/badge/license-MIT-green)](LICENSE)

A reproducible study of **how numerical integrators behave on black-hole geodesics**: how their errors grow over long integrations, which physical invariants they conserve, and what each unit of accuracy costs. Classical Runge–Kutta (RK4 and the adaptive pairs DP5 and DOP853) is compared with geometric integrators (Gauss–Legendre collocation of order 2, 4 and 6, and Tao's explicit symplectic scheme of order 2 and 4) on timelike and null geodesics of the Schwarzschild metric. Each runs in two formulations, the second-order geodesic equation and Hamilton's equations, and every error is measured against a closed-form reference built from elliptic integrals and functions.

![Light rays near the photon sphere and a precessing orbit](figures/F01_hero.png)

**Read more:** [report (PDF)](report/main.pdf) · [one-page summary](report/summary.pdf) · [documentation](https://nugiephysics.github.io/geodesic-integrators/) · [development log](docs/devlog.md)

## Key findings

- **Linear against bounded (RQ1).** Over $`10^4`$ orbits the Runge–Kutta methods drift linearly in the mass shell (RK4: $`1.6\times10^{-8}`$) and quadratically in phase. The symplectic methods stay bounded (GL2: $`5.0\times10^{-12}`$) with linear phase error. Sixth-order Gauss–Legendre reaches round-off and then follows Brouwer's $`N^{1/2}`$ law.
- **Symmetric is enough (RQ2).** In Hamilton's form, energy and angular momentum are exact for every method. Gauss–Legendre applied to the second-order system, where it is symmetric but not symplectic, conserves as well as in canonical variables: the largest mass-shell error of GL2 is $`5.0138\times10^{-12}`$ in both.
- **The ranking turns over with length (RQ3).** DOP853 is the cheapest method for short integrations at every accuracy. Over 1000 orbits it is still cheapest for a phase error of $`10^{-4}`$, but only GL3 reaches $`10^{-8}`$ (about $`6\times10^3`$ evaluations per orbit).
- **Six orbits on the photon sphere (RQ4).** Each decade of accuracy buys $`\ln 10/2\pi = 0.37`$ orbits, up to about six orbits in double precision for every method. Plain floating-point summation fakes a longer survival. At the ISCO the time to plunge grows as $`\lvert\varepsilon\rvert^{-1/4}`$ in the size of the error, and the sign of the error decides whether the orbit plunges at all.

![Mass-shell error over 10^4 periods](figures/F07_energy_long_term.png)

## Quick start

```bash
git clone https://github.com/NugiePhysics/geodesic-integrators.git
cd geodesic-integrators
uv sync --all-extras    # the locked environment
make test               # fast test suite (about 30 s once compiled)
make all                # every figure and table, then the report
```

`make all` integrates everything that is not cached yet (about 25 minutes on four cores from an empty cache), writes the figures to [`figures/`](figures/) and the summary tables to [`results/summary/`](results/summary/), and compiles the report with [Tectonic](https://tectonic-typesetting.github.io/). A repeated run reproduces every table and figure byte for byte.

```python
from geoint.experiments import RunSpec, solve
from geoint.testcases.schwarzschild import Eccentric

case = Eccentric(20.0, 0.5, n_orbits=10)       # (p, e) = (20, 0.5), ten radial periods
solution, form = solve(RunSpec.make(case, "b", "GL3", n_steps=10_000))
print(case.errors(solution, form)["error"])      # phase error against the exact orbit
```

The [notebooks](notebooks/) walk through the formulations, convergence, long-term behaviour and unstable orbits.

## Methods

<!-- Table T1: results/summary/t1_methods.md -->
| Method | Order | Explicit | Symplectic | Symmetric | Adaptive | Evaluations per step |
|---|---|---|---|---|---|---|
| RK4 | 4 | yes | no | no | no | 4 |
| DP5 (Dormand–Prince 5(4)) | 5 | yes | no | no | yes | 6 |
| DOP853 | 8 | yes | no | no | yes | 12 |
| GL1 (implicit midpoint) | 2 | no | yes | yes | no | 1 per iteration (6.8 measured) |
| GL2 (Gauss–Legendre) | 4 | no | yes | yes | no | 2 per iteration (11.2 measured) |
| GL3 (Gauss–Legendre) | 6 | no | yes | yes | no | 3 per iteration (14.2 measured) |
| Tao2 (extended phase space) | 2 | yes | yes | yes | no | 3 |
| Tao4 (triple jump) | 4 | yes | yes | yes | no | 9 |

## Validation

Every test case agrees with its exact reference to $`10^{-11}`$ or better with DOP853 at tolerance $`10^{-13}`$, except where the photon sphere amplifies errors (about $`10^{-8}`$). Full table: [T2](results/summary/t2_validation.md).

| Check | Numerical | Exact | Difference |
|---|---|---|---|
| Deflection $`\Delta\phi(r_\text{far})`$, $`b = 6`$, $`r_\text{far} = 1000`$ | 4.84898089193 | 4.84898089193 | $`6\times10^{-13}`$ |
| Periapsis advance, $`(p, e) = (20, 0.5)`$ | 1.23386 | 1.23386 | $`5\times10^{-13}`$ |
| Circular orbit $`r_c = 6`$: $`d\phi/dt = r_c^{-3/2}`$ | 0.0680414 | 0.0680414 | $`3\times10^{-15}`$ |
| Photon sphere: orbits before $`\lvert r - 3\rvert > 0.1`$, $`r_0 = 3 + 10^{-6}`$ | 1.93916 | 1.93916 | $`1\times10^{-8}`$ |
| Growth rate of $`\lvert r - 3\rvert`$ (Lyapunov exponent $`1/3\sqrt3`$) | 0.192515 | 0.19245 | $`7\times10^{-5}`$ |

## Project layout

```
src/geoint/
  metrics/          Metric interface and hand-written Schwarzschild kernels (Numba)
  formulations/     second-order (a) and Hamiltonian (b) vector fields, invariants
  integrators/      RK4, DP5, DOP853, Gauss-Legendre, Tao; one compiled core; SciPy oracle
  analytic/         exact references: deflection, bound orbits, photon sphere
  testcases/        TC0-TC6: initial data, events, references, error functions
  experiments/      run specifications, cache, sweeps, growth-law fits
  plotting/         one colour and marker per method
experiments/        one script per figure or table (make figures)
results/summary/    small CSV and Markdown tables (committed)
figures/            PDF (report) and PNG (README, docs)
report/             LaTeX report and one-page summary (make report)
docs/               theory, decisions (ADRs), development log, documentation site
tools/              SymPy derivation of Christoffel symbols and metric derivatives
```

The integrators never import physics. Adding a metric (Kerr is next) must not change a line in `integrators/` or `formulations/`; [Extending to Kerr](docs/extending.md) lists the steps.

## Development

```bash
uv run pre-commit install   # ruff lint and format on every commit
make lint                   # the same checks as CI
make test-nojit             # the tests with the Numba kernels run as plain Python
make test-slow              # the long integrations
make docs-serve             # the documentation site, live
```

## Citation

If you use this software, please cite it as described in [`CITATION.cff`](CITATION.cff).

## Related

This project grew out of [schwarzschild-geodesics](https://github.com/NugiePhysics/schwarzschild-geodesics), a Hamiltonian RK4 simulation of equatorial Schwarzschild geodesics whose validation table is reused here as a regression test.

## License

[MIT](LICENSE) © 2026 Nugie Saputra
