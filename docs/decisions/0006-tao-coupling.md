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

## Update (Phase 6): ω = 10⁻³ fails on the marginally stable ISCO

Phase 5 fixed ω = 10⁻³ for all geodesic experiments (figure F14): within a factor 1.2–5 of the most accurate ω on the orbits tested, and the copies stay within 3×10⁻⁶ over 300 periods. Phase 6 found a case where this choice fails. On the inclined ISCO orbit (figure F11, 45°), the radial dynamics is marginal: a difference between the copies along r is not rotated back by the orbital motion, as it is on an eccentric orbit, but drifts. The coupling rotation of 2ωh per step is too weak to stop it.

| ω | Steps per orbit | Orbits to plunge | max ‖r_q − r_x‖ | max \|δH\| |
|---|---|---|---|---|
| 10⁻³ | 100, 1000 | 3.2, 5.5 | 1.7, 1.7 | 2.5×10⁻², 6.5×10⁻² |
| 10⁻² | 100, 1000 | 12.5, 52.9 | 3×10⁻⁴, 0.17 | 2.2×10⁻⁵, 3.4×10⁻³ |
| 0.1 | 100, 1000 | 9.1, 106 | 1.5×10⁻³, 6.9×10⁻⁸ | 6.4×10⁻⁵, 4.4×10⁻⁹ |
| 1 | 100, 1000 | 6.0, 63.2 | 8×10⁻⁵, 1.7×10⁻⁸ | 4.4×10⁻⁴, 4.7×10⁻⁸ |

The comparisons keep the single global ω = 10⁻³, so that every figure shows the same method. F11 adds Tao4 with ω = 0.1 as a separate curve. The lesson for the report: ω is a problem-dependent parameter. A value tuned on a stable orbit does not transfer to a marginally stable one. Choosing ω is a real cost of Tao's method that the Gauss–Legendre methods do not have.

The inclined orbits of TC6 (figure F18, (p, e) = (20, 0.5), 1000 steps per period, 1000 periods) show the same weakness. At 30° Tao4 behaves: L² is bounded at 7×10⁻⁸ and the copies stay within 2×10⁻⁶. At 60° the copies separate after about 10 periods and the orbit is captured after 41. At 85° it is captured after 3 periods. Coupling φ as well makes every inclination worse. A larger ω helps only partly. At 60°, ω = 10⁻² holds the copies within 2×10⁻⁴ over 300 periods, but L² still drifts by 6×10⁻⁵, about 3000 times more than GL2. At 85° no ω from 10⁻³ to 1 keeps the run stable. Near the pole all methods are under-resolved at this step (RK4: δL² = 8×10⁻²), but only Tao fails outright.
