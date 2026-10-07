# 0005. One compiled integrator core: typed vector fields, kind dispatch, events by partial steps

- Status: Accepted
- Date: 2026-10-07

## Context

[ADR 0002](0002-numba-first.md) puts the integrator loops in Numba and composes them with vector fields passed as arguments. The first implementation did exactly that: a generic fixed-step driver taking a step function and an `rhs`. Every pair of driver and vector field then compiled separately, in 10–20 s each. Tests and experiments use about five vector fields (three toy problems and two formulations) and four kinds of methods, which would cost minutes of compilation in every process.

The obvious remedy, Numba's on-disk cache, failed for two reasons found while trying it:

1. A cached function cannot receive another jitted function *as a value*. Numba embeds a pointer to the dispatcher object and refuses to cache ("dynamic globals").
2. The cache is invalidated only when the file that defines a cached function changes, not when a callee in another file does. Stale machine code would survive edits to helper modules.

Events need their own decision. The predecessor located them by linear interpolation, which caps the measured order at 2 (roadmap pitfall 7). SciPy uses its dense output, but RK4, Gauss–Legendre and Tao have no interpolant of matching order.

## Decision

1. **Vector fields are first-class function values of one type**, `float64[::1] -> float64[::1]`. The compiled entry points `fixed(...)` and `adaptive(...)` declare that type in explicit signatures. They are compiled once and cached on disk, and they accept any vector field without recompiling. A lazily compiled `rhs` is converted automatically.
2. **One universal step function dispatches on an integer method kind** (`EXPLICIT`, `GAUSS`, `TAO`, `EMBEDDED`). All methods read their parameters from one tuple `P` of fixed type, and each ignores the fields it does not use. No jitted function is ever passed as a value inside the core.
3. **All compiled integrator code lives in one module**, `integrators/_compiled.py`, so that file-based cache invalidation is correct. The other integrator modules are thin Python classes.
4. **Events are data**, `Event(name, index, value, direction, stop_after)`, like the metric's `StopSurface` ([ADR 0004](0004-metric-interface.md)). A sign change of `y[index] - value` within a step is located by Illinois regula falsi on θ ∈ [0, 1], where every trial value is a genuine step of size θh *of the same method*. Event times are therefore accurate to the method's order for every integrator. A state that starts on the surface is not a crossing.
5. **Fixed steps use λ = λ₀ + n·h.** The step is shrunk so that n steps divide the interval exactly. The state update uses compensated (Kahan) summation by default; this can be switched off, for the Phase 6 experiment.
6. **The adaptive controller copies SciPy line by line** (initial step, SAFETY 0.9, factors 0.2 and 10, DOP853's combined error estimate), with the coefficients imported from SciPy itself.

## Consequences

Positive:

- A cold start compiles the core once, in about 60 s, and every later process starts in about 1 s. Solving a new problem costs only the compilation of its vector field.
- Our DP5 and DOP853 take the same steps and the same number of evaluations as SciPy, at every tolerance tested.
- Event locations converge at the order of the method, so convergence studies on event quantities (deflection angle, periapsis phase) are not limited by the locator.

Negative, and what we do about it:

- **Calls through a function pointer cannot be inlined.** They cost about 30 ns per call, measured on the oscillator. All our methods pay this equally, so wall-time comparisons between them stay fair (ADR 0003).
- **Locating an event costs evaluations**, typically 20–50 per event (about 6–10 trial partial steps), against SciPy's 3 for its DOP853 dense output. This is a few per cent of the total in event-heavy runs. It is included in `nfev`. The SciPy comparisons that check `nfev` are run without events.
- **The parameter tuple is untyped by method.** A wrong field fails only at run time. `make_params` builds it, and every method is covered by tests.
- **The disk cache can go stale only within `_compiled.py`** if Numba misses a change. `find src -name '*.nb[ic]' -delete` resets it.

## Alternatives considered

- **A driver per vector field without caching.** 10–20 s per pair, minutes per process.
- **Closures that capture the step function.** Closures cannot be cached either.
- **Dense-output event location for every method.** It needs an interpolant of the method's order. For Gauss–Legendre the collocation polynomial has only stage order s, which is too low for GL2 and GL3.
- **Hermite interpolation between steps.** Order 3 locally. That is enough for RK4 but not for GL3, DOP853 or long-term phase measurements.
