# Integrators

Every integrator receives a vector field and an initial state and knows nothing about the metric ([ADR 0001](decisions/0001-full-8d-phase-space.md)). All of them run in one compiled core ([ADR 0005](decisions/0005-compiled-integrator-core.md)) and are available by name from `geoint.integrators.get`.

<!-- include: results/summary/t1_methods.md -->

*Table T1 ([data](../results/summary/t1_methods.csv)). For the implicit methods the measured cost per step comes from $`10^4`$ periods of the orbit $`(20, 0.5)`$ at 1000 steps per period.*

## Explicit Runge–Kutta

- **RK4**, the classical four-stage method with a fixed step. After $`n`$ steps the parameter is $`\lambda_0 + nh`$, never a running sum.
- **DP5(4)** (Dormand & Prince 1980) and **DOP853** (Hairer, Nørsett & Wanner), the embedded pairs behind SciPy's `RK45` and `DOP853`. They are reimplemented in compiled code with SciPy's step-size controller copied line by line. From a generic starting point they take the same steps as SciPy and use the same number of evaluations at tolerances $`10^{-4}`$–$`10^{-12}`$; step sizes agree to $`10^{-7}`$. SciPy itself is kept as an oracle (`scipy-RK45`, `scipy-DOP853`).

Events (periapsis passages, crossings of a radius) are located by root finding on partial steps of the method itself, so they are as accurate as the method. Linear interpolation would cap the measured order at two (pitfall 7).

## Gauss–Legendre collocation

The $`s`$-stage Gauss methods have order $`2s`$, and they are symplectic and symmetric. GL1 is the implicit midpoint rule. One step solves for the stage increments $`Z_i = h\sum_j a_{ij}F(y_n + Z_j)`$:

- by **fixed-point iteration**, started from the collocation polynomial of the previous step;
- **until the correction stops decreasing**, Hairer's rule for round-off-limited integrations, with a guard: stagnation counts only once the correction is within 100 ulp of the increments, because on rotation-like problems the correction can wobble long before round-off;
- the update $`y_{n+1} = y_n + d^\mathsf{T}Z`$ with $`d = b^\mathsf{T}A^{-1}`$ avoids a final evaluation and is added with **compensated (Kahan) summation**.

An optional `iter_tol` stops the iteration at an absolute tolerance instead. Experiment F15 shows that tolerances above $`10^{-12}`$ make GL2 drift (pitfall 10).

| | $`s = 2`$ (GL2) | $`s = 3`$ (GL3) |
|---|---|---|
| nodes $`c`$ | $`\tfrac12 \mp \tfrac{\sqrt3}{6}`$ | $`\tfrac12 - \tfrac{\sqrt{15}}{10},\ \tfrac12,\ \tfrac12 + \tfrac{\sqrt{15}}{10}`$ |
| weights $`b`$ | $`\tfrac12,\ \tfrac12`$ | $`\tfrac{5}{18},\ \tfrac49,\ \tfrac{5}{18}`$ |

The coefficients $`a_{ij}`$ are computed in closed form in `geoint.integrators.tableaus`.

## Tao's method

For a nonseparable $`H(q, p)`$, Tao (2016) integrates two copies with

```math
\tilde H(q, p, x, y) = H(q, y) + H(x, p) + \frac{\omega}{2}\left(\lvert q - x\rvert^2 + \lvert p - y\rvert^2\right).
```

The flows of the first two terms are explicit shifts, and the third rotates $`(q - x, p - y)`$ by the angle $`2\omega\delta`$. The second-order method is the Strang composition $`\phi_A^{h/2}\phi_B^{h/2}\phi_C^{h}\phi_B^{h/2}\phi_A^{h/2}`$, the fourth-order one its triple jump (Yoshida 1990). Consecutive $`\phi_A`$ flows share a cached gradient, which gives 3 and 9 evaluations per step.

Two choices are not fixed by the method ([ADR 0006](decisions/0006-tao-coupling.md)):

- **Only the non-cyclic coordinates $`(r, \theta)`$ are coupled.** Coupling $`t`$ and $`\phi`$ lets the time copies drift apart and feeds that drift into $`p_t`$, so $`E`$ is no longer conserved.
- **$`\omega = 10^{-3}`$** (figure F14) in every experiment. It fails on the marginally stable ISCO and on orbits inclined by 60° or more, where the copies separate.

The reported state is the first copy $`(q, p)`$; the full extended state is kept for diagnostics.

## Cost

The cost of a run is its number of vector-field evaluations, including every fixed-point iteration and every gradient of Tao's method ([ADR 0003](decisions/0003-nfev-as-cost-metric.md)). Wall times are compared only between our compiled implementations, with one thread, after compilation, as a median of repeats.
