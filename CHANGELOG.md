# Changelog

Versions follow the phases of the project. Details, including what went wrong, are in the [development log](docs/devlog.md).

## 1.0.0 (2026-10-07)

Report, documentation and release.

- Report (`report/main.tex`, RevTeX, 27 pages) and one-page research summary, compiled with Tectonic by `make report`. Every table and headline number in them is generated from `results/summary/` by `scripts/make_tables.py`.
- Documentation site (MkDocs Material): theory, integrators, test cases, results by research question, pitfalls, reproducibility, a guide to adding Kerr, and the API reference. It is deployed to GitHub Pages by `docs.yml`.
- Four walkthrough notebooks, executed in CI with nbmake.
- `make all` regenerates every figure and table and the report. The weekly `reproduce.yml` workflow runs it from a clean checkout.
- Figure F1 (the geodesics of the study) and tables T1 (integrators), T7 (recommendations) and T8 (computing environment).

## 0.9.0 (2026-10-07)

Exact references, short and long experiments (Phases 4–6).

- Exact references for deflection, bound orbits and the photon sphere; the test-case registry TC1–TC6; the experiment harness with cached runs and provenance.
- Short integrations: convergence, tolerance proportionality, work-precision, deflection, windings, precession, Tao's ω, implicit iteration, step sizes (F2–F6, F12–F16, T2, T3, T5).
- Long integrations and unstable orbits: growth laws, invariants, symmetric against symplectic, work-precision over 1000 periods, photon sphere, ISCO, inclined orbits, round-off (F7–F11, F17–F20, T4, T5b, T6).
- Tao's method fails on the ISCO and on inclined orbits with ω = 10⁻³ (ADR 0006).

## 0.3.0 (2026-10-07)

Gauss–Legendre GL1–GL3 and Tao's method of order 2 and 4, verified on toy problems.

## 0.2.0 (2026-10-07)

The compiled integrator core: RK4, DP5(4) and DOP853 (step for step identical to SciPy), events by partial steps, the SciPy oracle.

## 0.1.0 (2026-10-07)

The metric interface, hand-written Schwarzschild kernels checked against SymPy, both formulations, and mass-shell initial data.
