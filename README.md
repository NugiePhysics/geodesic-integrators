# geodesic-integrators

[![CI](https://github.com/NugiePhysics/geodesic-integrators/actions/workflows/ci.yml/badge.svg)](https://github.com/NugiePhysics/geodesic-integrators/actions/workflows/ci.yml)
[![Python](https://img.shields.io/badge/python-3.12%2B-blue)](pyproject.toml)
[![License: MIT](https://img.shields.io/badge/license-MIT-green)](LICENSE)

A reproducible study of **how numerical integrators behave on black-hole geodesics**: error growth over long integrations, conservation of physical invariants, and cost per unit of accuracy. Classical Runge–Kutta (RK4, RK45, DOP853) is compared with geometric integrators (Gauss–Legendre collocation and Tao's explicit symplectic scheme) on timelike and null geodesics of the Schwarzschild metric, in two formulations:

- **(a) second order**, the geodesic equation $`\ddot x^\mu = -\Gamma^\mu_{\alpha\beta}\dot x^\alpha \dot x^\beta`$ in $`(x^\mu, u^\mu)`$;
- **(b) Hamiltonian**, $`H = \tfrac12 g^{\mu\nu}p_\mu p_\nu`$ in $`(x^\mu, p_\mu)`$.

Every result is measured against a closed-form reference (elliptic integrals and functions), and the code is laid out so that Kerr can be added without touching the integrators.

> **Status: Phase 5 of 7 (short-integration experiments).** All integrators (RK4, DP5, DOP853, Gauss–Legendre GL1–GL3, Tao of order 2 and 4) are implemented and verified, every test case has an exact reference, and the convergence and work-precision figures are in [`figures/`](figures/) (`make figures` regenerates them). Long-term experiments come next. Progress, results and what went wrong are logged in [`docs/devlog.md`](docs/devlog.md); design decisions in [`docs/decisions/`](docs/decisions/).

## Development

The environment is locked with [uv](https://docs.astral.sh/uv/):

```bash
git clone https://github.com/NugiePhysics/geodesic-integrators.git
cd geodesic-integrators
uv sync --all-extras        # runtime + dev group + optional JAX
uv run pre-commit install   # ruff lint + format on every commit
make test                   # fast test suite
make lint                   # the same checks as CI
NUMBA_DISABLE_JIT=1 uv run pytest   # run the Numba kernels as plain Python, for debugging
```

## Related

This project grew out of [schwarzschild-geodesics](https://github.com/NugiePhysics/schwarzschild-geodesics), a Hamiltonian RK4 simulation of equatorial Schwarzschild geodesics whose validation table is reused here as a regression test.

## License

[MIT](LICENSE) © 2026 Nugie Saputra
