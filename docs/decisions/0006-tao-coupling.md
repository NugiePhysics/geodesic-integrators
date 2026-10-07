# 0006. Tao's method: couple only the non-cyclic coordinates, report the first copy

- Status: Accepted
- Date: 2026-10-07

## Context

Tao's method integrates two copies `(q, p)` and `(x, y)` of the state. A rotation by 2ωδ in the planes `(q_k − x_k, p_k − y_k)` couples them. Two choices are left open by the method itself: which coordinates to couple, and which combination of the copies to report.

Phase 0 found with a throwaway script that coupling all eight components lets E of one copy fluctuate. Phase 3 measured the effect on the orbit (p, e) = (20, 0.5) over three radial periods with Tao4:

| Coupling | ω | Steps per orbit | Error in φ at the end |
|---|---|---|---|
| all coordinates | 0.1 | 400, 1600 | integration fails (copies separate) |
| all coordinates | 1, 20 | 1600 | 1.7e-4, 9.8e-6 |
| all coordinates | 5, 100 | 1600 | integration fails |
| r and θ only | 0.1 | 400, 1600 | 1.0e-2, 3.0e-5 |
| r and θ only | 1, 20, 100 | 1600 | 2.2e-4, 2.2e-4, 2.3e-4 |

With full coupling the time copies drift apart (up to 10⁴ in t), because their rates are evaluated at different copies. The rotation then converts this difference into a change of p_t, so E is no longer conserved and the copies can run away. With only (r, θ) coupled, p_t and p_φ stay identical in both copies to the last bit. The cyclic positions t, φ of the two copies drift apart without any effect on the dynamics.

## Decision

1. **Couple only the coordinates the metric depends on.** The integrator stays physics-agnostic: it takes a boolean `coupling` mask. The experiment harness builds the mask from `Metric.cyclic_coords`, which gives `(False, True, True, False)` for Schwarzschild and Kerr.
2. **Report the first copy `(q, p)`**, and keep the full extended state in `Solution.extended` for the diagnostics ‖q − x‖, ‖p − y‖. Events act on the first copy.
3. **Choose ω per step size, from data.** Phase 5 (figure F14) measures the error against ω at fixed h, and the copy separation over long runs. The value used in the comparisons is taken from there and recorded in the devlog.

## Consequences

- E and L_z are conserved exactly by Tao in formulation (b), as for every Runge–Kutta method. The RQ2 statement "conserved by all methods in (b)" then holds without exception.
- The coupling mask is the only place where Tao uses knowledge of the metric. It enters through the harness, not through the integrator.
- Reporting one copy rather than the average keeps the reported state a point on a trajectory of the extended Hamiltonian flow. The average is not.

## Alternatives considered

- **Couple all coordinates** (Tao 2016; FANTASY). Rejected on the evidence above.
- **Report the copy average.** It removes the O(h²) asymmetry between the copies, but it is no longer the image of a symplectic map, and the reported E would be an average.
- **Fix ω once** (for example ω = 20, as in Tao's examples). The good range depends on h, because the rotation angle 2ωγh enters through trigonometric functions: with large ωh the coupling can become nearly the identity, and the copies decouple.
