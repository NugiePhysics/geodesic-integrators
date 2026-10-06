# The Two Formulations of the Geodesic Problem

This document derives the two systems of ordinary differential equations that every integrator in this project is tested on:

- **(a) the second-order formulation**, the geodesic equation in the state $`y = (x^\mu, u^\mu) \in \mathbb{R}^8`$;
- **(b) the Hamiltonian formulation**, Hamilton's equations for $`H = \tfrac12 g^{\mu\nu}p_\mu p_\nu`$ in the state $`y = (x^\mu, p_\mu) \in \mathbb{R}^8`$.

It then lists the conserved quantities ($`E`$, $`L_z`$, $`L^2`$, $`H`$) in both sets of variables and explains which of them an integrator can or cannot preserve exactly. That last point decides which conservation plots are meaningful (research question RQ2 of the [roadmap](../../ROADMAP.md)).

**Conventions.** Signature $`(-,+,+,+)`$, geometrized units $`G = c = 1`$, Schwarzschild coordinates $`x^\mu = (t, r, \theta, \phi)`$. A dot is the derivative with respect to the affine parameter $`\lambda`$. $`\epsilon = 1`$ for timelike geodesics (then $`\lambda = \tau`$, the proper time) and $`\epsilon = 0`$ for null geodesics. The code sets $`M = 1`$; it is kept explicit here so that dimensions can be checked.

The derivation of §§1–4 holds for any static, spherically symmetric metric function $`f(r)`$. Schwarzschild is the special case $`f = 1 - 2M/r`$, $`f' = 2M/r^2`$. The code does **not** inherit this restriction: its `Metric` interface takes the full $`g^{\mu\nu}(x)`$ and $`\partial_\alpha g^{\mu\nu}(x)`$ so that Kerr ($`g^{t\phi} \neq 0`$) fits without changes ([ADR 0001](../decisions/0001-full-8d-phase-space.md)).

---

## 1. The Metric

```math
ds^2 = -f(r)\,dt^2 + \frac{dr^2}{f(r)} + r^2\left(d\theta^2 + \sin^2\theta\,d\phi^2\right),
\qquad f(r) = 1 - \frac{2M}{r}. \tag{1}
```

The metric is diagonal, so its inverse is the component-wise reciprocal:

```math
g_{\mu\nu} = \operatorname{diag}\!\left(-f,\ \frac1f,\ r^2,\ r^2\sin^2\theta\right),
\qquad
g^{\mu\nu} = \operatorname{diag}\!\left(-\frac1f,\ f,\ \frac1{r^2},\ \frac1{r^2\sin^2\theta}\right). \tag{2}
```

Coordinate singularities: $`f = 0`$ at the horizon $`r = 2M`$ ($`g_{rr}, g^{tt} \to \infty`$) and $`\sin\theta = 0`$ on the polar axis ($`g^{\phi\phi} \to \infty`$). Both are properties of the coordinates, not of the spacetime, but both are real numerical hazards (roadmap §7, pitfalls 1–2).

---

## 2. Geodesics from a Variational Principle

Affinely parametrized geodesics are the stationary curves of the action $`S = \int \mathcal{L}\,d\lambda`$ with

```math
\mathcal{L}(x, \dot x) = \tfrac12\, g_{\mu\nu}(x)\, \dot x^\mu \dot x^\nu . \tag{3}
```

The Euler–Lagrange equations are

```math
\frac{d}{d\lambda}\left(g_{\mu\nu}\dot x^\nu\right) = \tfrac12\, \partial_\mu g_{\alpha\beta}\, \dot x^\alpha \dot x^\beta . \tag{4}
```

Expanding the left-hand side, $`g_{\mu\nu}\ddot x^\nu + \partial_\alpha g_{\mu\beta}\,\dot x^\alpha \dot x^\beta`$, symmetrizing the last term in $`(\alpha, \beta)`$ and multiplying by $`g^{\sigma\mu}`$ gives the geodesic equation

```math
\ddot x^\sigma + \Gamma^\sigma_{\alpha\beta}\,\dot x^\alpha \dot x^\beta = 0,
\qquad
\Gamma^\sigma_{\alpha\beta} = \tfrac12\, g^{\sigma\mu}\left(\partial_\alpha g_{\mu\beta} + \partial_\beta g_{\mu\alpha} - \partial_\mu g_{\alpha\beta}\right). \tag{5}
```

$`\mathcal{L}`$ does not depend explicitly on $`\lambda`$, so $`\mathcal{L}`$ itself is conserved. Its value fixes the type of the geodesic and the normalization of $`\lambda`$:

```math
g_{\mu\nu}\,\dot x^\mu \dot x^\nu = -\epsilon . \tag{6}
```

For null geodesics (6) leaves the scale of $`\lambda`$ free. This project fixes it by setting $`E = 1`$, so that the angular momentum equals the impact parameter, $`L = b`$.

---

## 3. Formulation (a): The Second-Order System

Writing $`u^\mu = \dot x^\mu`$ turns (5) into eight first-order equations,

```math
\dot x^\mu = u^\mu, \qquad \dot u^\mu = -\Gamma^\mu_{\alpha\beta}(x)\, u^\alpha u^\beta, \tag{7}
```

with the constraint $`\mathcal{N}(x, u) \equiv g_{\mu\nu}(x)\,u^\mu u^\nu = -\epsilon`$.

### 3.1 Christoffel symbols

For the metric (1) exactly nine Christoffel symbols are independent and non-zero (thirteen counting the symmetry $`\Gamma^\sigma_{\alpha\beta} = \Gamma^\sigma_{\beta\alpha}`$). Each follows from (5) with a single non-zero metric derivative, for example $`\Gamma^r_{tt} = \tfrac12 g^{rr}(-\partial_r g_{tt}) = \tfrac12 f f'`$.

| Symbol | General $`f(r)`$ | Schwarzschild |
|---|---|---|
| $`\Gamma^t_{tr} = \Gamma^t_{rt}`$ | $`\dfrac{f'}{2f}`$ | $`\dfrac{M}{r(r-2M)}`$ |
| $`\Gamma^r_{tt}`$ | $`\tfrac12 f f'`$ | $`\dfrac{M(r-2M)}{r^3}`$ |
| $`\Gamma^r_{rr}`$ | $`-\dfrac{f'}{2f}`$ | $`-\dfrac{M}{r(r-2M)}`$ |
| $`\Gamma^r_{\theta\theta}`$ | $`-rf`$ | $`-(r-2M)`$ |
| $`\Gamma^r_{\phi\phi}`$ | $`-rf\sin^2\theta`$ | $`-(r-2M)\sin^2\theta`$ |
| $`\Gamma^\theta_{r\theta} = \Gamma^\theta_{\theta r}`$ | $`\dfrac1r`$ | $`\dfrac1r`$ |
| $`\Gamma^\theta_{\phi\phi}`$ | $`-\sin\theta\cos\theta`$ | $`-\sin\theta\cos\theta`$ |
| $`\Gamma^\phi_{r\phi} = \Gamma^\phi_{\phi r}`$ | $`\dfrac1r`$ | $`\dfrac1r`$ |
| $`\Gamma^\phi_{\theta\phi} = \Gamma^\phi_{\phi\theta}`$ | $`\cot\theta`$ | $`\cot\theta`$ |

### 3.2 Explicit equations

Inserting the table into (7) (the factor 2 comes from the two orderings of the mixed lower indices):

```math
\begin{align}
\dot u^t &= -\frac{f'}{f}\, u^t u^r, \tag{8a}\\
\dot u^r &= -\tfrac12 f f'\,(u^t)^2 + \frac{f'}{2f}\,(u^r)^2 + r f\left[(u^\theta)^2 + \sin^2\theta\,(u^\phi)^2\right], \tag{8b}\\
\dot u^\theta &= -\frac2r\, u^r u^\theta + \sin\theta\cos\theta\,(u^\phi)^2, \tag{8c}\\
\dot u^\phi &= -\frac2r\, u^r u^\phi - 2\cot\theta\, u^\theta u^\phi . \tag{8d}
\end{align}
```

For Schwarzschild, (8a) and (8b) read

```math
\dot u^t = -\frac{2M}{r(r-2M)}\,u^t u^r, \qquad
\dot u^r = -\frac{M(r-2M)}{r^3}(u^t)^2 + \frac{M}{r(r-2M)}(u^r)^2 + (r-2M)\left[(u^\theta)^2 + \sin^2\theta\,(u^\phi)^2\right]. \tag{9}
```

The right-hand side is a quadratic form in $`u`$ with $`x`$-dependent coefficients. The cost of one evaluation is dominated by the contraction $`\Gamma^\mu_{\alpha\beta}u^\alpha u^\beta`$; a generic implementation (the one Kerr will use) does the full $`4\times4\times4`$ contraction, whereas (8) exploits sparsity.

---

## 4. Formulation (b): Hamilton's Equations

### 4.1 Legendre transform

The canonical momentum conjugate to $`x^\mu`$ is the covariant velocity,

```math
p_\mu = \frac{\partial \mathcal{L}}{\partial \dot x^\mu} = g_{\mu\nu}\,\dot x^\nu , \tag{10}
```

and the Legendre transform of (3) is the super-Hamiltonian (MTW §25.2)

```math
H(x, p) = p_\mu \dot x^\mu - \mathcal{L} = \tfrac12\, g^{\mu\nu}(x)\, p_\mu p_\nu . \tag{11}
```

Hamilton's equations are

```math
\dot x^\mu = \frac{\partial H}{\partial p_\mu} = g^{\mu\nu} p_\nu, \qquad
\dot p_\mu = -\frac{\partial H}{\partial x^\mu} = -\tfrac12\, \partial_\mu g^{\alpha\beta}\, p_\alpha p_\beta . \tag{12}
```

Because $`H = \mathcal{L}`$ numerically along a solution, the constraint (6) becomes

```math
H(x, p) = -\frac{\epsilon}{2}. \tag{13}
```

Only $`g^{\mu\nu}`$ and its first derivatives enter (12); no Christoffel symbols are needed. This is why formulation (b) is the natural one for Tao's method, which needs $`\partial H/\partial x`$ and $`\partial H/\partial p`$ as separate functions.

### 4.2 Schwarzschild

```math
H = \frac12\left[-\frac{p_t^2}{f} + f\,p_r^2 + \frac{p_\theta^2}{r^2} + \frac{p_\phi^2}{r^2\sin^2\theta}\right]. \tag{14}
```

Of the $`4 \times 4 \times 4`$ array $`\partial_\alpha g^{\mu\nu}`$ only five entries are non-zero:

```math
\partial_r g^{tt} = \frac{f'}{f^2}, \quad
\partial_r g^{rr} = f', \quad
\partial_r g^{\theta\theta} = -\frac{2}{r^3}, \quad
\partial_r g^{\phi\phi} = -\frac{2}{r^3\sin^2\theta}, \quad
\partial_\theta g^{\phi\phi} = -\frac{2\cos\theta}{r^2\sin^3\theta}. \tag{15}
```

Hamilton's equations (12) then read

```math
\begin{align}
\dot t &= -\frac{p_t}{f}, & \dot p_t &= 0, \tag{16a}\\
\dot r &= f\,p_r, & \dot p_r &= -\frac{f'}{2}\left(\frac{p_t^2}{f^2} + p_r^2\right) + \frac{1}{r^3}\left(p_\theta^2 + \frac{p_\phi^2}{\sin^2\theta}\right), \tag{16b}\\
\dot\theta &= \frac{p_\theta}{r^2}, & \dot p_\theta &= \frac{\cos\theta}{r^2\sin^3\theta}\,p_\phi^2, \tag{16c}\\
\dot\phi &= \frac{p_\phi}{r^2\sin^2\theta}, & \dot p_\phi &= 0. \tag{16d}
\end{align}
```

For Schwarzschild the radial force is $`\dot p_r = -M p_t^2/(r-2M)^2 - M p_r^2/r^2 + L^2/r^3`$, with $`L^2`$ from (20) below. In the equatorial plane this is eq. (24) of the derivation in `schwarzschild-geodesics`.

$`H`$ is **not separable**: the "kinetic" term $`\tfrac12 g^{\mu\nu}(x)p_\mu p_\nu`$ depends on position, so $`H \neq T(p) + V(x)`$ and the Störmer–Verlet/leapfrog family cannot be applied in its explicit form. This is why the project uses implicit Gauss–Legendre methods and Tao's extended phase space (roadmap §2.4).

### 4.3 Equivalence of (a) and (b)

Differentiating $`g^{\alpha\beta}g_{\beta\gamma} = \delta^\alpha_\gamma`$ gives $`\partial_\mu g^{\alpha\beta} = -g^{\alpha\rho}g^{\beta\sigma}\,\partial_\mu g_{\rho\sigma}`$. With $`p = g\,u`$ the momentum equation in (12) becomes

```math
\dot p_\mu = \tfrac12\, \partial_\mu g_{\rho\sigma}\, u^\rho u^\sigma . \tag{17a}
```

Differentiating $`p_\mu = g_{\mu\nu}u^\nu`$ along the curve gives instead

```math
\dot p_\mu = \partial_\rho g_{\mu\sigma}\, u^\rho u^\sigma + g_{\mu\nu}\,\dot u^\nu . \tag{17b}
```

Equating the two and symmetrizing in $`(\rho, \sigma)`$ yields $`g_{\mu\nu}\dot u^\nu = -\tfrac12(\partial_\rho g_{\mu\sigma} + \partial_\sigma g_{\mu\rho} - \partial_\mu g_{\rho\sigma})u^\rho u^\sigma`$, which is (7). The two formulations therefore describe the same curves; they differ only in the variables the integrator sees.

Equation (17b) is also the form of the **unit test** that ties the two right-hand sides together in Phase 1: for random $`(x, u)`$ and $`p = g(x)u`$, the position part of (b) must equal $`u`$, and its momentum part must equal $`\partial_\rho g_{\mu\nu}u^\rho u^\nu + g_{\mu\nu}\dot u^\nu`$ with $`\dot u`$ taken from (a).

---

## 5. Conserved Quantities

### 5.1 Definitions

If $`\xi^\mu`$ is a Killing vector, $`\xi_\mu u^\mu = \xi^\mu p_\mu`$ is constant along every geodesic. In formulation (b) this is immediate when $`\xi = \partial_\mu`$ for a coordinate $`x^\mu`$ that does not appear in $`H`$ (a *cyclic* coordinate): then $`\dot p_\mu = -\partial H/\partial x^\mu = 0`$. Schwarzschild has four Killing vectors: $`\partial_t`$ (staticity) and the three generators of rotations.

| Quantity | Formulation (a) | Formulation (b) | Origin |
|---|---|---|---|
| Energy $`E`$ | $`f\,u^t`$ | $`-p_t`$ | Killing vector $`\partial_t`$ |
| Axial angular momentum $`L_z`$ | $`r^2\sin^2\theta\,u^\phi`$ | $`p_\phi`$ | Killing vector $`\partial_\phi`$ |
| Total angular momentum squared $`L^2`$ | $`r^4\left[(u^\theta)^2 + \sin^2\theta\,(u^\phi)^2\right]`$ | $`p_\theta^2 + \dfrac{p_\phi^2}{\sin^2\theta}`$ | all three rotations |
| Mass shell | $`\tfrac12 g_{\mu\nu}u^\mu u^\nu`$ | $`H = \tfrac12 g^{\mu\nu}p_\mu p_\nu`$ | $`\lambda`$-independence of $`\mathcal{L}`$, $`H`$ |

with the constant values $`E_0`$, $`L_{z,0}`$, $`L^2_0`$ and $`-\epsilon/2`$ respectively.

**$`L^2`$ from the rotation generators.** The Killing vectors of $`SO(3)`$ are

```math
\xi_x = -\sin\phi\,\partial_\theta - \cot\theta\cos\phi\,\partial_\phi, \qquad
\xi_y = \cos\phi\,\partial_\theta - \cot\theta\sin\phi\,\partial_\phi, \qquad
\xi_z = \partial_\phi . \tag{18}
```

The cross terms in $`L_x^2 + L_y^2`$ cancel, leaving

```math
L^2 = L_x^2 + L_y^2 + L_z^2 = p_\theta^2 + \cot^2\theta\,p_\phi^2 + p_\phi^2 = p_\theta^2 + \frac{p_\phi^2}{\sin^2\theta}. \tag{19}
```

A direct check from (16): $`\tfrac{d}{d\lambda}L^2 = 2p_\theta\dot p_\theta - 2p_\phi^2\cos\theta\,\dot\theta/\sin^3\theta = 0`$.

**Link to Kerr.** Carter's constant reduces for $`a = 0`$ to $`Q = p_\theta^2 + \cot^2\theta\,p_\phi^2 = L^2 - L_z^2`$. Testing $`\delta(L^2)`$ on inclined Schwarzschild orbits (test case TC6) therefore exercises exactly the code path that will monitor $`Q`$ in Kerr.

**Radial equation.** Substituting $`E`$, $`L^2`$ into (13) and using $`u^r = f p_r`$:

```math
(u^r)^2 = E^2 - V_\text{eff}(r), \qquad V_\text{eff}(r) = f(r)\left(\epsilon + \frac{L^2}{r^2}\right). \tag{20}
```

This is the form used to set up initial data (§7) and to find turning points analytically.

**Relative errors used in the experiments** (roadmap §1.4): timelike $`\delta H = |H + \tfrac12|/\tfrac12`$; null $`\delta H = |H|/E^2`$, which is invariant under the rescaling $`p \to \alpha p`$ that the free normalization of null $`\lambda`$ allows; $`\delta E = |E - E_0|/E_0`$, and likewise for $`L_z`$ and $`L^2`$.

### 5.2 Why $`p_t`$ and $`p_\phi`$ are conserved *exactly* in formulation (b)

Every Runge–Kutta method, explicit or implicit, advances the state as

```math
y_{n+1} = y_n + h\sum_{i=1}^s b_i\, F(Y_i), \qquad Y_i = y_n + h\sum_{j=1}^s a_{ij}\, F(Y_j). \tag{21}
```

In formulation (b) the $`p_t`$- and $`p_\phi`$-components of $`F`$ are $`-\tfrac12 \partial_t g^{\alpha\beta}p_\alpha p_\beta`$ and $`-\tfrac12 \partial_\phi g^{\alpha\beta}p_\alpha p_\beta`$. Their coefficients $`\partial_t g^{\alpha\beta}`$ and $`\partial_\phi g^{\alpha\beta}`$ are identically zero in the code, so these components evaluate to the floating-point number $`0`$, not merely to something small. The increment $`h\sum_i b_i \cdot 0`$ is then exactly zero and

```math
p_t^{(n+1)} = p_t^{(n)}, \qquad p_\phi^{(n+1)} = p_\phi^{(n)} \qquad \text{bit for bit.} \tag{22}
```

The argument does not depend on the order, on the step size, on adaptivity, on how many fixed-point iterations an implicit method performs (the stage values $`Y_i`$ also have unchanged $`p_t`$, $`p_\phi`$), or on compensated summation. It is the floating-point version of the general theorem that all Runge–Kutta methods preserve linear invariants (GNI Thm. IV.1.5).

**Consequence for the experiments.** Plots of $`\delta E`$ and $`\delta L_z`$ in formulation (b) are flat at zero for every Runge–Kutta-type method. That is a property of the variables, not a sign of quality, and must not be presented as a result. The informative invariants in (b) are $`H`$ and, for inclined orbits, $`L^2`$. The informative $`E`$ and $`L_z`$ tests live in formulation (a), where $`E = f(r)\,u^t`$ and $`L_z = r^2\sin^2\theta\,u^\phi`$ are *nonlinear* functions of the state. There, every method shows a non-zero error, and how that error grows is the actual comparison (RQ2).

**Exception: Tao's extended phase space.** Tao's method integrates two copies $`(q, p)`$ and $`(x, y)`$ with $`\tilde H = H(q, y) + H(x, p) + \omega\,\tfrac12\left(|q - x|^2 + |p - y|^2\right)`$. The flows of $`H(q, y)`$ and $`H(x, p)`$ leave $`p_t`$ and $`y_t`$ untouched for the reason above. The coupling flow, however, rotates the pair $`(q - x,\ p - y)`$ by the angle $`2\omega\delta`$. The time copies drift apart, $`q^t \neq x^t`$, because $`\dot q^t = g^{tt}(x)\,p_t`$ and $`\dot x^t = g^{tt}(q)\,y_t`$ are evaluated at different copies. So the rotation feeds $`q^t - x^t`$ into $`p_t - y_t`$. Only the sums $`p_t + y_t`$ and $`p_\phi + y_\phi`$ are invariant under all three sub-flows (up to rounding in the rotation). For Tao, therefore:

- $`E`$ and $`L_z`$ of each individual copy fluctuate, at the level of the copy separation;
- the copy average $`\tfrac12(p_t + y_t)`$ is conserved to rounding error;
- every reported $`E`$ must state which of the two it is.

An implementation choice that removes the effect is to apply the coupling rotation only to the non-cyclic components $`(r, \theta)`$. The cyclic momenta of the two copies then stay identical, and the decoupled cyclic positions never feed back into the dynamics. Whether this changes Tao's accuracy is to be decided in Phase 3. It is recorded here so that the RQ2 claim "conserved by all methods" is not overstated.

### 5.3 What each class of method can preserve

| Invariant | Structure in (a) | Structure in (b) | Preserved exactly by |
|---|---|---|---|
| $`E`$, $`L_z`$ | nonlinear in $`(r, \theta, u)`$ | linear, and $`\dot p_t = \dot p_\phi \equiv 0`$ | (b): every RK method, bit for bit (§5.2). (a): no method |
| $`L^2`$ | nonlinear | nonlinear ($`\theta`$-dependent) | no method |
| $`H`$ / mass shell | nonlinear ($`g`$ depends on $`x`$) | nonlinear ($`g^{-1}`$ depends on $`x`$) | no method in this study |

Gauss–Legendre methods preserve every *quadratic* invariant $`y^\mathsf{T} C\, y`$ (GNI Thm. IV.2.1), but none of the invariants above is quadratic in the state, because all of them carry $`x`$-dependent coefficients. What symplectic methods offer for $`H`$ and $`L^2`$ is not exact conservation but a **bounded** error, $`O(h^p)`$ with no secular drift over exponentially long times, from backward error analysis (GNI IX.8). Methods that preserve $`H`$ exactly (discrete-gradient or averaged-vector-field methods) exist but are outside the scope of this study.

---

## 6. Symplecticity and Reversibility

**Formulation (b)** is a canonical Hamiltonian system: its flow preserves $`\omega = dx^\mu \wedge dp_\mu`$. Symplectic Runge–Kutta methods (Gauss–Legendre) applied in $`(x, p)`$ produce symplectic maps, and backward error analysis guarantees near-conservation of $`H`$.

**Formulation (a)** is the same flow written in the non-canonical variables $`(x, u)`$. The change of variables $`(x, u) \mapsto (x, g(x)\,u)`$ pulls $`\omega`$ back to

```math
dx^\mu \wedge d\left(g_{\mu\nu}u^\nu\right) = g_{\mu\nu}\, dx^\mu \wedge du^\nu + \partial_\rho g_{\mu\nu}\, u^\nu\, dx^\mu \wedge dx^\rho ,
```

a symplectic form with $`x`$-dependent coefficients. A symplectic Runge–Kutta method is guaranteed to preserve a *constant* symplectic form, and only for vector fields that are Hamiltonian with respect to it. System (a) is Hamiltonian only with respect to the form above, which a Runge–Kutta method has no reason to preserve. In (a) Gauss–Legendre is therefore **symmetric but not symplectic** in any useful sense.

**Reversibility.** Both systems are reversible with respect to the linear involution $`\rho`$ that flips the velocity or momentum, $`\rho(x, u) = (x, -u)`$ and $`\rho(x, p) = (x, -p)`$. Writing either system as $`\dot y = F(y)`$:

```math
\rho\, F(y) = -F(\rho\, y), \tag{23}
```

because $`\dot x`$ is odd and $`\dot u`$ (or $`\dot p`$) is even in the velocity. Symmetric methods ($`\Phi_{-h} = \Phi_h^{-1}`$: implicit midpoint, Gauss–Legendre, Strang and triple-jump compositions; *not* RK4, RK45 or DOP853) are then $`\rho`$-reversible. For integrable reversible systems with a non-degenerate torus structure, reversible methods show the same long-time behaviour as symplectic ones: bounded errors in the action variables and linear growth of the phase error (GNI XI). Schwarzschild geodesics are integrable ($`E`$, $`L_z`$, $`L^2`$, $`H`$ in involution). The comparison "Gauss–Legendre in (a) versus in (b)" (figure F17) thus tests whether reversibility alone is enough here.

---

## 7. Initial Data and Stopping

The constraint (13) is quadratic in the momenta. Given a position $`x^\mu`$, the constants $`E`$, $`L_z`$ and either $`p_\theta`$ or $`L^2`$, it is solved for the remaining component, by default $`p_r`$:

```math
p_r = \pm\frac{\sqrt{E^2 - V_\text{eff}(r)}}{f(r)} . \tag{24}
```

Near a turning point $`E^2 - V_\text{eff}`$ is a difference of nearly equal numbers. Phase 1 will evaluate it in a cancellation-free form, for example by factoring the radial polynomial $`R(r) = r^3\left(E^2 - V_\text{eff}\right)`$ with its known roots, or by starting exactly at a turning point with $`p_r = 0`$. The state for formulation (a) then follows from $`u^\mu = g^{\mu\nu}p_\nu`$, so both formulations start from the same physical point.

Integrations stop at $`r = 2M(1 + \delta)`$ ("captured"), at $`r = r_\text{far}`$ ("escaped"), or at a requested $`\lambda`$ or event count. Events (periapsis passages $`p_r = 0`$ with $`\dot p_r > 0`$, crossings of $`r_\text{far}`$) are located by root finding on dense output, never by linear interpolation (roadmap §7, pitfall 7).

---

## 8. Relation to the Reduced System of `schwarzschild-geodesics`

Setting $`\theta = \pi/2`$, $`p_\theta = 0`$ in (16) and treating $`E = -p_t`$, $`L = p_\phi`$ as fixed parameters gives the four-equation system $`(t, r, \phi, p_r)`$ of the earlier project. That reduction is exact for the true flow but removes $`p_t`$ and $`p_\phi`$ from the state, so their conservation cannot be measured. It also hard-codes a diagonal, static, spherically symmetric metric. This project keeps all eight components, so the same code handles inclined orbits and, later, Kerr ([ADR 0001](../decisions/0001-full-8d-phase-space.md)).

In floating point $`\cos(\pi/2) \approx 6\times10^{-17} \neq 0`$. An equatorial orbit integrated in 8D therefore acquires $`p_\theta = O(10^{-16})`$ through (16c). This is a tilt of the orbital plane at the rounding level, not a bug (pitfall 3).

---

## 9. Reference Values

Closed-form numbers used by the tests ($`M = 1`$):

| Quantity | Value |
|---|---|
| Horizon | $`r_h = 2`$ |
| Photon sphere, critical impact parameter | $`r = 3`$, $`b_c = 3\sqrt3 \approx 5.196152`$ |
| Lyapunov exponent of the photon sphere (coordinate time) | $`\lambda_L = 1/(3\sqrt3)`$, equal to its orbital frequency |
| Circular timelike orbit at $`r_c`$ | $`E = \dfrac{1 - 2/r_c}{\sqrt{1 - 3/r_c}}`$, $`L = \dfrac{\sqrt{r_c}}{\sqrt{1 - 3/r_c}}`$, $`\Omega = \dfrac{d\phi}{dt} = r_c^{-3/2}`$ |
| ISCO | $`r = 6`$, $`E = \sqrt{8/9}`$, $`L = 2\sqrt3`$ |

---

## References

- S. Carroll, *Spacetime and Geometry* (2004), Ch. 3 (geodesics, Christoffel symbols) and Ch. 5 (Schwarzschild, Killing vectors).
- C. Misner, K. Thorne, J. Wheeler, *Gravitation* (1973), §25.2 (super-Hamiltonian).
- E. Poisson, *A Relativist's Toolkit* (2004), Ch. 1.
- S. Chandrasekhar, *The Mathematical Theory of Black Holes* (1983), Ch. 3.
- B. Carter, *Global structure of the Kerr family of gravitational fields*, Phys. Rev. 174, 1559 (1968).
- E. Hairer, C. Lubich, G. Wanner, *Geometric Numerical Integration* (2nd ed., 2006): IV.1–IV.2 (linear and quadratic invariants), VI.4 (symplectic RK), IX.8 (long-time energy conservation), XI (reversible methods).
- M. Tao, *Explicit symplectic approximation of nonseparable Hamiltonians*, Phys. Rev. E 94, 043303 (2016).

Bibliographic details to be checked against NASA ADS before they enter `report/refs.bib`.
