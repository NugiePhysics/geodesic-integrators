# Test cases and references

Each test case (`geoint.testcases.schwarzschild`) fixes the initial data in both formulations, the span of the integration, the events that stop it, an exact reference and an error. Errors are always measured on **physical quantities** (an angle at an event, the phase at a known proper time), never on the raw state, because the states of formulations (a) and (b) differ. Units: $`G = c = M = 1`$.

| ID | Case | Parameters | Reference | Headline error |
|---|---|---|---|---|
| TC0 | Toy problems | harmonic oscillator; Kepler; $`H = \tfrac12(q^2 + 1)(p^2 + 1)`$ | exact solution or invariant | order, symplecticity, $`\lvert\Delta H\rvert`$ |
| TC1 | Light deflection | $`b \in [5.3, 10^3]`$, $`r_\text{far} = 10^3`$–$`10^4`$, $`E = 1`$, $`L = b`$ | elliptic integral to the same $`r_\text{far}`$ | $`\lvert\Delta\phi - \Delta\phi_\text{exact}(r_\text{far})\rvert`$ |
| TC1b | Near-critical photons | $`b = b_c(1 + 10^{-k})`$, $`k = 1, \dots, 10`$ | windings $`\approx -\ln(b/b_c - 1)/2\pi`$ (Bozza 2002) | windings |
| TC2 | Photon sphere | $`r_0 = 3 + 10^{-k}`$ or $`r_0 = 3`$ in an inclined plane | elliptic integral from the turning point; $`\lambda_L = 1/(3\sqrt3)`$ | orbits until $`\lvert r - 3\rvert > 0.1`$ |
| TC3 | Circular orbits | $`r_c \in \{6, 7, 10, 20, 100\}`$ | $`E`$, $`L`$, $`\Omega = r_c^{-3/2}`$ | phase $`\phi - (L/r_c^2)\lambda`$ |
| TC3m | Inclined ISCO | $`r = 6`$, 45° | marginal stability | orbits until $`\lvert r - 6\rvert > 1`$ |
| TC4 | Eccentric orbits | $`(p, e) \in \{(100, 0.5), (20, 0.5), (7.5, 0.5)\}`$ | $`\Phi`$ and $`r(\phi)`$ in closed form | phase at $`\lambda = NT_\tau`$; shape $`\lvert r - r(\phi)\rvert`$ |
| TC5 | Long term | TC4 $`(20, 0.5)`$ and TC3 $`r_c = 10`$ for $`10^3`$–$`10^4`$ periods | $`\phi`$ at the $`N`$-th periapsis $`= N\Phi`$ | growth laws of $`\lvert\delta H\rvert`$ and the phase |
| TC6 | Inclined orbits | $`(20, 0.5)`$ at 30°, 60°, 85° | $`L^2`$, $`L_z`$, the orbital plane | $`\delta(L^2)`$, tilt of the plane |

## References

**Light deflection.** With $`u = 1/r`$ and the roots $`u_1 < 0 < u_2 < u_3`$ of $`2u^3 - u^2 + b^{-2}`$, a photon from $`r_\text{far}`$ back to $`r_\text{far}`$ sweeps

```math
\Delta\phi(r_\text{far}) = \frac{4F(\psi\,|\,m)}{\sqrt{2(u_3 - u_1)}}, \qquad
m = \frac{u_2 - u_1}{u_3 - u_1}, \qquad
\sin^2\psi = \frac{(u_3 - u_1)(u_2 - u_\text{far})}{(u_2 - u_1)(u_3 - u_\text{far})}.
```

For $`r_\text{far} \to \infty`$ this is Darwin's (1959) angle (Iyer & Petters 2007, eq. 6). Three independent evaluations agree to double precision for $`b \in [5.2, 10^4]`$: this form, Darwin's form in terms of the periapsis, and a direct 30-digit quadrature. The comparison always uses the same $`r_\text{far}`$ as the integration ([pitfall 18](pitfalls.md)).

**Bound orbits.** For the orbit $`(p, e)`$ (Cutler, Kennefick & Poisson 1994), $`r_p = p/(1 + e)`$, $`r_a = p/(1 - e)`$, and $`(du/d\phi)^2 = 2(u - u_1)(u - u_2)(u - u_3)`$ with $`u_1 = (1 - e)/p`$, $`u_2 = (1 + e)/p`$ and $`u_3 = \tfrac12 - 2/p`$. Then

```math
\Phi = \frac{4K(m)}{\sqrt{2(u_3 - u_1)}}, \qquad
u(\phi) = u_1 + (u_2 - u_1)\,\mathrm{cd}^2\!\left(\sqrt{\tfrac{u_3 - u_1}{2}}\,\phi \,\middle|\, m\right),
\qquad m = \frac{u_2 - u_1}{u_3 - u_1}.
```

The radial periods in proper and coordinate time, and $`\Phi`$ a second time, come from a 30-digit quadrature in the angle $`\chi`$ of $`r = \tfrac12(r_a + r_p) - \tfrac12(r_a - r_p)\cos\chi`$, which removes both turning-point singularities. Conversions: $`E^2 = [(p - 2)^2 - 4e^2]/[p(p - 3 - e^2)]`$, $`L^2 = p^2/(p - 3 - e^2)`$; the separatrix is $`p = 6 + 2e`$.

**Photon sphere.** An offset $`\delta_0`$ grows like $`\cosh(\lambda_L t)`$ with $`\lambda_L = 1/(3\sqrt3)`$, a factor $`e^{2\pi}`$ per orbit, so linear theory predicts $`\operatorname{arccosh}(0.1/\delta_0)/2\pi`$ orbits before $`\lvert r - 3\rvert > 0.1`$. The exact count is an elliptic integral from the turning point. With $`r_0`$ known exactly the other two roots follow without solving a cubic, which matters because $`b - b_c = O(\delta_0^2)`$ falls below double precision for $`\delta_0 < 10^{-8}`$.

**Fixed points.** The equatorial photon sphere and every equatorial circular orbit are exact fixed points of every Runge–Kutta map: the radial components of the vector field vanish, and $`t`$, $`\phi`$ are linear in $`\lambda`$. Their errors measure round-off, not truncation. To test truncation the experiments tilt the orbital plane, so that errors in $`L^2`$ reach $`r`$.

## Error measures

- Mass shell: $`\delta H = \lvert H + \tfrac12\rvert/\tfrac12`$ (timelike), $`\lvert H\rvert/E^2`$ (null). The latter is invariant under the free rescaling of null $`\lambda`$. In (a) the same quantity is the norm $`\lvert g_{\mu\nu}u^\mu u^\nu + \epsilon\rvert`$.
- $`\delta E = \lvert E - E_0\rvert/E_0`$, and likewise $`\delta L_z`$ and $`\delta(L^2)`$.
- Order of convergence: least-squares slope of $`\log(\text{error})`$ against $`\log h`$ over the asymptotic range, above the round-off floor.
- Growth laws: the log-log slope of the per-period maximum over the last decade of the run, classified as 0 (bounded), ½ (random walk), 1 (linear) or 2 (quadratic). For the phase also the signed fit $`\phi_N - N\Phi = c_1N + c_2N^2`$.
- Cost: vector-field evaluations, including implicit iterations ([ADR 0003](decisions/0003-nfev-as-cost-metric.md)).

## Validation

Every reference is checked by two independent paths, and every test case agrees with a tight DOP853 integration (tolerance $`10^{-13}`$, formulation (b)) to $`10^{-11}`$ or better. The exceptions are the cases amplified by the photon sphere, which agree to about $`10^{-8}`$.

<!-- include: results/summary/t2_validation.md -->

*Table T2 ([data](../results/summary/t2_validation.csv)).*
