# Reading List

The core reading, in the suggested order, with what each piece is needed for. A box is ticked when the reading is done and its key points are in the devlog. Bibliographic details are checked against NASA ADS before they enter `report/refs.bib`.

## Core (Phase 0 exit criterion)

| Done | Reading | Take-away for this project | Needed by |
|---|---|---|---|
| [ ] | Hairer, Nørsett & Wanner, *Solving ODEs I*, II.1–II.4 | Order conditions, local vs global error, embedded pairs and step-size control | Phase 2 |
| [ ] | Hairer, Lubich & Wanner, *Geometric Numerical Integration* (GNI), Ch. I | The motivating examples: why symplectic methods win on long runs | Phase 3 |
| [ ] | GNI VI.1–VI.4 | Symplectic maps, symplectic RK (Gauss), symmetry | Phase 3 |
| [ ] | Tao (2016), PRE 94, 043303 | Extended phase space, the role of $`\omega`$, composition to order 4 | Phase 3 |
| [ ] | GNI IX.1–IX.3, IX.8 | Backward error analysis, modified Hamiltonian, why $`\lvert\Delta H\rvert`$ stays bounded | Phases 3, 6 |
| [ ] | GNI X.3 | Linear vs quadratic growth of the phase error on integrable systems (RQ1) | Phase 6 |
| [ ] | GNI VIII.5–VIII.6 and Hairer, McLachlan & Razakarivony (2008), BIT 48, 231 | Implementing implicit RK: fixed-point iteration to stagnation, compensated summation, Brouwer's law | Phase 3 |

## During Phase 4 (analytic references)

| Done | Reading | Take-away | Needed by |
|---|---|---|---|
| [ ] | Chandrasekhar, *The Mathematical Theory of Black Holes*, Ch. 3 | Exact Schwarzschild orbits in elliptic functions | Phase 4 |
| [ ] | Darwin (1959), Proc. R. Soc. A 249, 180; Iyer & Petters (2007), GRG 39, 1563 | Exact light-deflection angle | Phase 4 |
| [ ] | Bozza (2002), PRD 66, 103001 | Strong-deflection limit near $`b_c`$ | Phase 4 |
| [ ] | Cutler, Kennefick & Poisson (1994), PRD 50, 3816 | $`(p, e)`$ parametrization and the separatrix | Phase 4 |
| [ ] | DLMF Ch. 19 and 22 | Elliptic integrals and Jacobi functions; the $`k`$ vs $`m = k^2`$ convention (pitfall 17) | Phase 4 |

## Background, as needed

- Carroll, *Spacetime and Geometry*, Ch. 3 and 5; MTW §25.2; Poisson, *A Relativist's Toolkit*, Ch. 1. Already condensed in [`theory/formulations.md`](theory/formulations.md).
- Yoshida (1990), Phys. Lett. A 150, 262; McLachlan & Quispel (2002), Acta Numerica 11, 341: composition methods.
- Christian & Chan (2021), ApJ 909, 67 (FANTASY); Wang, Sun, Liu & Wu (2021), ApJ 907, 66; Seyrich & Lukes-Gerakopoulos (2012), PRD 86, 124013: prior work on geodesic integrators.
- Cardoso et al. (2009), PRD 79, 064016: Lyapunov exponents of circular orbits (TC2).
- Leimkuhler & Reich, *Simulating Hamiltonian Dynamics* (2004): a gentler introduction than GNI.
