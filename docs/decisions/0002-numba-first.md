# 0002. Numba-first numerical core, SciPy as oracle, JAX optional

- Status: Accepted
- Date: 2026-10-06

## Context

The integrators advance an 8-component state through $`10^5`$–$`10^8`$ small steps. One right-hand-side evaluation is on the order of a hundred floating-point operations. In pure Python with NumPy, the per-call overhead of microseconds dwarfs that work, so wall-clock comparisons would measure the interpreter rather than the algorithms. SciPy's `solve_ivp` has the same issue: its stepping loop is Python.

The project needs:

- its own implementations of DP5(4), DOP853, Gauss–Legendre and Tao, to time them fairly against each other;
- an independent, trusted reference implementation, to check correctness;
- a way to step through the numerics in a debugger.

Keeping a pure-Python reference next to a compiled copy doubles maintenance, and the two copies drift apart.

## Decision

1. **One code path in Numba.** Metric functions, right-hand sides and integrator loops are `@njit` functions. For debugging, `NUMBA_DISABLE_JIT=1` runs exactly the same source as plain Python (`make test-nojit`).
2. **Composition by closures, not classes.** A formulation is a factory that returns a jitted `rhs`, closing over jitted metric functions. Integrator loops are `@njit` functions that take that `rhs` as a first-class argument. `jitclass` is avoided. Python classes are thin wrappers that hold jitted callables plus metadata (order, symplectic, ...). `tests/test_smoke.py` checks this pattern against the installed Numba.
3. **Closures capture only structural parameters.** Values fixed for a whole study (for example $`M`$, later the Kerr spin $`a`$) are compiled in as constants. Per-run values ($`h`$, tolerances, Tao's $`\omega`$, initial data) are passed as arguments, so that sweeps do not trigger recompilation.
4. **No `fastmath`.** It would allow reassociation, which breaks compensated summation and makes results depend on the compiler.
5. **SciPy is the correctness oracle**, not a timing competitor. `solve_ivp` (RK45, DOP853) checks our own embedded methods (same steps, same nfev) and provides event-location cross-checks.
6. **JAX is an optional extra** (`pip install geodesic-integrators[jax]`). It is used for autodiff verification of the hand-written $`\partial g^{\mu\nu}`$ and $`\Gamma`$, for fast prototyping of new metrics (Kerr–(A)dS), and later for `vmap` ray tracing. Any module using it must call `jax.config.update("jax_enable_x64", True)` first.

## Consequences

Positive:

- One source of truth for the numerics, and wall-clock timings that compare algorithms, not interpreters.
- Bugs can be chased in interpreted mode without a second implementation.

Negative, and what we do about it:

- **Compilation latency.** Compiling each new type signature takes seconds. Benchmarks warm up before timing, heavy kernels use `cache=True`, and the test suite accepts a fixed compile overhead.
- **Numba's subset of Python and NumPy.** No SciPy, limited exceptions and strings inside kernels. Kernels report failure through integer status codes, which the Python wrapper turns into exceptions or `Solution.status`.
- **Interpreted and compiled runs need not agree bit for bit**, because libm and LLVM intrinsics and summation order can differ. The Level-3 cross-check (roadmap §5) therefore compares the two modes to a few ulp per step rather than demanding identical output.
- **Numba lags new CPython releases.** It decides the supported Python range together with NumPy, SciPy and JAX. Today the floor is 3.12 (NumPy ≥ 2.5, SciPy ≥ 1.18 and JAX require it). Numba 0.68 covers 3.10–3.15.

## Alternatives considered

- **Pure NumPy.** Overhead-dominated for 8D states, as described above.
- **Cython, or C++ through pybind11.** Fast, but adds a build system and a second language for a single-author research code.
- **JAX as the main path.** Adaptive step control and implicit iterations need `lax.while_loop` with awkward control flow, compile times are long, and float64 is opt-in. It is excellent for derivatives and batching, which is why it stays as an optional extra.
- **Julia (DifferentialEquations.jl).** Superb integrator library, but the point here is to implement and dissect the methods, in the Python stack used across astrophysics.
