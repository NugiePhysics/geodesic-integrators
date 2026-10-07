"""All compiled integrator code, in one module.

Numba's on-disk cache is invalidated only when the file that defines a cached function changes,
not when a callee in another file does. Keeping every compiled function of the integrators in
this one file makes ``cache=True`` safe (ADR 0005).

*Step functions* have the signature ``step(rhs, y, err, h, P, comp) -> (y_new, err_new, nfev,
niter, ok)``; ``err`` is the compensation vector of compensated (Kahan) summation, used only if
``comp`` is true (GNI VIII.5). The universal :func:`step` dispatches on an integer method
``kind``, and every method reads its parameters from one tuple ``P`` of fixed type (see
:data:`PARAMS`). Passing step *functions* as values would make Numba embed dispatcher pointers,
which cannot be cached.

The *entry points* at the end have explicit signatures in which ``rhs`` is a first-class
function type. They are therefore compiled once, for every vector field with the signature
``float64[::1] -> float64[::1]``, and cached on disk.
"""

from __future__ import annotations

import numpy as np
from numba import njit, types

VEC = types.float64[::1]
MAT = types.float64[:, ::1]
IVEC = types.int64[::1]
BVEC = types.boolean[::1]
#: Type of every vector field the integrators accept.
RHS = types.FunctionType(VEC(VEC))
f8, i8, b1 = types.float64, types.int64, types.boolean
EVENTS = (IVEC, VEC, IVEC, IVEC)

#: Method kinds understood by :func:`step`.
EXPLICIT, GAUSS, TAO, EMBEDDED = 0, 1, 2, 3

#: Parameter tuple ``P``; each method reads only its own fields, the rest are empty arrays.
#: ``0 A`` (MAT) tableau or Gauss matrix · ``1 b`` (VEC) weights, or Gauss ``d = b A^-1`` ·
#: ``2 A_inv`` (MAT) · ``3 X`` (MAT) Gauss extrapolation · ``4 memory`` (MAT) previous Gauss
#: stages and step · ``5 fractions`` (VEC) Tao composition · ``6 coupling`` (BVEC) Tao ·
#: ``7 scalars`` (VEC) ``[iter_tol, max_iter, omega]`` · ``8, 9, 10`` (VEC, VEC, BVEC) Tao
#: gradient cache · ``11 f`` (VEC) derivative at the step start (embedded pairs).
PARAMS = types.Tuple((MAT, VEC, MAT, MAT, MAT, VEC, BVEC, VEC, VEC, VEC, BVEC, VEC))

# --- state update, events, fixed-step driver --------------------------------------------------


EVENT_THETA_TOL = 1e-14


@njit
def update(y, err, delta, comp):
    """``y + delta``, with compensated summation if ``comp``."""
    if comp:
        d = delta + err
        y_new = y + d
        return y_new, (y - y_new) + d
    return y + delta, err


@njit
def all_finite(y):
    for v in y:  # noqa: SIM110 (Numba cannot compile all() over a generator)
        if not np.isfinite(v):
            return False
    return True


@njit
def crossing(g0, g1, direction):
    """Whether ``g`` changes sign from ``g0`` to ``g1`` in the requested direction.

    A start exactly on the surface (``g0 == 0``) is not a crossing, so an orbit started at a
    turning point does not report a spurious event at ``lam = 0`` (roadmap pitfall 7).
    """
    up = g0 < 0.0 and g1 >= 0.0
    down = g0 > 0.0 and g1 <= 0.0
    if direction > 0:
        return up
    if direction < 0:
        return down
    return up or down


@njit
def locate(kind, rhs, y, h, P, index, value, g0, g1, y1):
    """Fraction ``theta`` of the step at which ``y[index] == value``, and the state there.

    Every trial value of ``theta`` is evaluated with a genuine step of size ``theta h`` of the
    same method, so the event is located to the order of the method, not of an interpolant.
    The root is bracketed in ``[0, 1]`` and found with the Illinois variant of regula falsi.
    """
    zero = np.zeros(y.size)
    a, b = 0.0, 1.0
    ga, gb = g0, g1
    theta, yt = 1.0, y1
    if g1 == 0.0:
        return theta, yt, 0
    nfev = 0
    side = 0
    for _ in range(100):
        theta = (a * gb - b * ga) / (gb - ga)
        if not (a < theta < b):
            theta = 0.5 * (a + b)
        yt, _, fe, _, _ = step(kind, rhs, y, zero, theta * h, P, False)
        nfev += fe
        gt = yt[index] - value
        if gt == 0.0:
            break
        if (gt > 0.0) == (gb > 0.0):
            b, gb = theta, gt
            if side == -1:
                ga *= 0.5
            side = -1
        else:
            a, ga = theta, gt
            if side == 1:
                gb *= 0.5
            side = 1
        if b - a <= EVENT_THETA_TOL:
            break
    return theta, yt, nfev


@njit
def scan_events(kind, rhs, y, y_new, h, P, ev_index, ev_value, ev_dir):
    """All event crossings within one step, sorted by ``theta``."""
    n_ev = ev_index.size
    ids = np.empty(n_ev, np.int64)
    thetas = np.empty(n_ev)
    states = np.empty((n_ev, y.size))
    n = 0
    nfev = 0
    for e in range(n_ev):
        g0 = y[ev_index[e]] - ev_value[e]
        g1 = y_new[ev_index[e]] - ev_value[e]
        if crossing(g0, g1, ev_dir[e]):
            theta, ye, fe = locate(kind, rhs, y, h, P, ev_index[e], ev_value[e], g0, g1, y_new)
            nfev += fe
            ids[n] = e
            thetas[n] = theta
            states[n] = ye
            n += 1
    order = np.argsort(thetas[:n])
    return ids[:n][order], thetas[:n][order], states[:n][order], nfev


@njit
def grow(lams, ys, ids):
    """Event buffers with twice the capacity."""
    n = lams.size
    new_lams = np.empty(2 * n)
    new_ys = np.empty((2 * n, ys.shape[1]))
    new_ids = np.empty(2 * n, np.int64)
    new_lams[:n] = lams
    new_ys[:n] = ys
    new_ids[:n] = ids
    return new_lams, new_ys, new_ids


@njit
def process_events(
    kind, rhs, y, y_new, lam, h, P, ev_index, ev_value, ev_dir, ev_stop, counts, buf
):
    """Record the events of one step. Returns the index of a terminating event or -1.

    ``buf`` is ``(lams, ys, ids, n_found)``; it is returned updated because buffers may grow.
    """
    lams, ys, ids, n_found = buf
    hit_ids, thetas, states, nfev = scan_events(
        kind, rhs, y, y_new, h, P, ev_index, ev_value, ev_dir
    )
    stop, stop_lam, stop_y = -1, lam, y_new
    for k in range(hit_ids.size):
        e = hit_ids[k]
        if n_found == lams.size:
            lams, ys, ids = grow(lams, ys, ids)
        lams[n_found] = lam + thetas[k] * h
        ys[n_found] = states[k]
        ids[n_found] = e
        n_found += 1
        counts[e] += 1
        if ev_stop[e] > 0 and counts[e] >= ev_stop[e]:
            stop, stop_lam, stop_y = e, lam + thetas[k] * h, states[k]
            break
    return stop, stop_lam, stop_y, nfev, (lams, ys, ids, n_found)


@njit
def fixed_step_driver(
    kind, rhs, y0, lam0, h, n_steps, P, comp, ev_index, ev_value, ev_dir, ev_stop, save_every
):
    """March ``n_steps`` steps of size ``h``; ``lam = lam0 + n h``, never ``lam += h``."""
    dim = y0.size
    cap = n_steps // save_every + 3
    lam_out = np.empty(cap)
    y_out = np.empty((cap, dim))
    counts_out = np.zeros((cap, 2), np.int64)  # cumulative (nfev, niter) at each saved point
    lam_out[0] = lam0
    y_out[0] = y0
    n_out = 1
    buf = (np.empty(16), np.empty((16, dim)), np.empty(16, np.int64), 0)
    counts = np.zeros(ev_index.size, np.int64)
    y = y0.copy()
    err = np.zeros(dim)
    lam = lam0
    nfev = 0
    niter = 0
    status = 0
    n = 0
    while n < n_steps:
        y_new, err_new, fe, it, ok = step(kind, rhs, y, err, h, P, comp)
        nfev += fe
        niter += it
        if not ok:
            status = -4
            break
        if not all_finite(y_new):
            status = -1
            break
        if ev_index.size > 0:
            stop, stop_lam, stop_y, fe, buf = process_events(
                kind, rhs, y, y_new, lam, h, P, ev_index, ev_value, ev_dir, ev_stop, counts, buf
            )
            nfev += fe
            if stop >= 0:
                n += 1
                lam = stop_lam
                y = stop_y
                status = 1 + stop
                break
        n += 1
        lam = lam0 + n * h
        y = y_new
        err = err_new
        if n % save_every == 0 and n < n_steps:
            lam_out[n_out] = lam
            y_out[n_out] = y
            counts_out[n_out, 0] = nfev
            counts_out[n_out, 1] = niter
            n_out += 1
    lam_out[n_out] = lam
    y_out[n_out] = y
    counts_out[n_out, 0] = nfev
    counts_out[n_out, 1] = niter
    n_out += 1
    lams, ys, ids, n_found = buf
    return (
        lam_out[:n_out],
        y_out[:n_out],
        status,
        nfev,
        niter,
        n,
        ids[:n_found],
        lams[:n_found],
        ys[:n_found],
        counts_out[:n_out],
    )


# --- explicit Runge-Kutta --------------------------------------------------------------------


@njit
def explicit_rk_step(rhs, y, err, h, P, comp):
    """One step of the explicit method with tableau ``(A, b) = (P[0], P[1])``."""
    A, b = P[0], P[1]
    s, n = b.size, y.size
    K = np.empty((s, n))
    K[0] = rhs(y)
    for i in range(1, s):
        dy = np.zeros(n)
        for j in range(i):
            if A[i, j] != 0.0:
                dy += A[i, j] * K[j]
        K[i] = rhs(y + h * dy)
    delta = np.zeros(n)
    for i in range(s):
        if b[i] != 0.0:
            delta += b[i] * K[i]
    y_new, err_new = update(y, err, h * delta, comp)
    return y_new, err_new, s, 0, True


# --- Gauss-Legendre --------------------------------------------------------------------------


#: An iteration that stops decreasing has converged only if its increment is within this many
#: units of rounding of the stage increments; above it the non-monotone convergence of the
#: fixed-point map (complex eigenvalues of ``h A J``) would otherwise stop it early.
STAGNATION_ULPS = 100.0
EPS = np.finfo(np.float64).eps


@njit
def gauss_step(rhs, y, err, h, P, comp):
    A, d, A_inv, X, memory = P[0], P[1], P[2], P[3], P[4]
    iter_tol, max_iter = P[7][0], int(P[7][1])
    s, n = d.size, y.size
    Z = np.zeros((s, n))
    if memory[s, 0] == h:
        for i in range(s):
            for j in range(s):
                Z[i] += X[i, j] * memory[j]
    F = np.empty((s, n))
    nfev = 0
    it = 0
    previous = np.inf
    first = -1.0
    dnorm = np.inf
    ok = False
    while it < max_iter:
        it += 1
        for i in range(s):
            F[i] = rhs(y + Z[i])
        nfev += s
        dnorm = 0.0
        znorm = 0.0
        for i in range(s):
            z_new = np.zeros(n)
            for j in range(s):
                z_new += A[i, j] * F[j]
            z_new *= h
            for k in range(n):
                scale = 1.0 + abs(y[k])
                dnorm = max(dnorm, abs(z_new[k] - Z[i, k]) / scale)
                znorm = max(znorm, abs(z_new[k]) / scale)
            Z[i] = z_new
        if first < 0.0:
            first = dnorm
        if dnorm == 0.0:
            ok = True
            break
        if iter_tol > 0.0:
            if dnorm <= iter_tol:
                ok = True
                break
        elif dnorm >= previous and dnorm <= STAGNATION_ULPS * EPS * (1.0 + znorm):
            ok = True
            break
        if not np.isfinite(dnorm) or dnorm > 1e3 * max(first, 1e-300):
            break
        previous = dnorm
    if it == max_iter and np.isfinite(dnorm) and dnorm <= 1e-10:
        ok = True  # converged as far as rounding allows, but never stagnated cleanly
    delta = np.zeros(n)
    for i in range(s):
        delta += d[i] * Z[i]
    for i in range(s):
        memory[i] = 0.0
        for j in range(s):
            memory[i] += A_inv[i, j] * Z[j]
    memory[s, 0] = h
    y_new, err_new = update(y, err, delta, comp)
    return y_new, err_new, nfev, it, ok


# --- Tao -------------------------------------------------------------------------------------


GAMMA1 = 1.0 / (2.0 - 2.0 ** (1.0 / 3.0))
GAMMA2 = 1.0 - 2.0 * GAMMA1


@njit
def _add(z, e, k, inc, comp):
    """``z[k] += inc``, compensated if ``comp``."""
    if comp:
        d = inc + e[k]
        new = z[k] + d
        e[k] = (z[k] - new) + d
        z[k] = new
    else:
        z[k] += inc


@njit
def _gradient(rhs, v, cache):
    """``rhs(v)`` with a one-entry cache ``cache = (v_last, rhs_last, valid)``; returns nfev."""
    v_last, f_last, valid = cache
    if valid[0] and np.array_equal(v, v_last):
        return f_last, 0
    f = rhs(v)
    v_last[:] = v
    f_last[:] = f
    valid[0] = True
    return f, 1


@njit
def _flow_A(rhs, z, e, delta, n, comp, cache):
    v = np.concatenate((z[0:n], z[3 * n : 4 * n]))  # (q, y)
    f, fe = _gradient(rhs, v, cache)
    for k in range(n):
        _add(z, e, n + k, delta * f[n + k], comp)  # p -= delta dH/dq
        _add(z, e, 2 * n + k, delta * f[k], comp)  # x += delta dH/dp
    return fe


@njit
def _flow_B(rhs, z, e, delta, n, comp, cache):
    v = np.concatenate((z[2 * n : 3 * n], z[n : 2 * n]))  # (x, p)
    f, fe = _gradient(rhs, v, cache)
    for k in range(n):
        _add(z, e, k, delta * f[k], comp)  # q += delta dH/dp
        _add(z, e, 3 * n + k, delta * f[n + k], comp)  # y -= delta dH/dq
    return fe


@njit
def _flow_C(z, e, delta, n, omega, coupling, comp):
    angle = omega * delta
    sin2 = np.sin(2.0 * angle)
    cos2_m1 = -2.0 * np.sin(angle) ** 2  # cos(2 angle) - 1 without cancellation
    for k in range(n):
        if not coupling[k]:
            continue
        s = z[k] - z[2 * n + k]  # q - x
        w = z[n + k] - z[3 * n + k]  # p - y
        dq = 0.5 * (cos2_m1 * s + sin2 * w)
        dp = 0.5 * (cos2_m1 * w - sin2 * s)
        _add(z, e, k, dq, comp)
        _add(z, e, 2 * n + k, -dq, comp)
        _add(z, e, n + k, dp, comp)
        _add(z, e, 3 * n + k, -dp, comp)


@njit
def tao_step(rhs, z, err, h, P, comp):
    fractions, coupling, omega = P[5], P[6], P[7][2]
    cache = (P[8], P[9], P[10])
    n = z.size // 4
    z_new = z.copy()
    e = err.copy() if comp else np.zeros(z.size)
    nfev = 0
    for g in fractions:
        delta = g * h
        nfev += _flow_A(rhs, z_new, e, 0.5 * delta, n, comp, cache)
        nfev += _flow_B(rhs, z_new, e, 0.5 * delta, n, comp, cache)
        _flow_C(z_new, e, delta, n, omega, coupling, comp)
        nfev += _flow_B(rhs, z_new, e, 0.5 * delta, n, comp, cache)
        nfev += _flow_A(rhs, z_new, e, 0.5 * delta, n, comp, cache)
    return z_new, (e if comp else err), nfev, 0, True


# --- adaptive Dormand-Prince -----------------------------------------------------------------


SAFETY = 0.9
MIN_FACTOR = 0.2
MAX_FACTOR = 10.0


@njit
def _rms(x):
    return np.sqrt(np.sum(x * x) / x.size)


@njit
def _stages(rhs, y, f, h, A, n_stages, K):
    """Stages ``K[0..n_stages-1]`` of one step (``K[0] = f`` is given); returns nfev."""
    n = y.size
    K[0] = f
    for s in range(1, n_stages):
        dy = np.zeros(n)
        for j in range(s):
            if A[s, j] != 0.0:
                dy += A[s, j] * K[j]
        K[s] = rhs(y + dy * h)
    return n_stages - 1


@njit
def _combine(K, B, n_stages):
    out = np.zeros(K.shape[1])
    for j in range(n_stages):
        out += B[j] * K[j]
    return out


@njit
def embedded_partial_step(rhs, y, err, h, P, comp):
    """Step of size ``h`` without error control, for event location: ``A, B, f = P[0, 1, 11]``."""
    A, B, f = P[0], P[1], P[11]
    n_stages = B.size
    K = np.empty((n_stages, y.size))
    nfev = _stages(rhs, y, f, h, A, n_stages, K)
    return y + h * _combine(K, B, n_stages), err, nfev, 0, True


@njit
def select_initial_step(rhs, y0, f0, interval, max_step, order, rtol, atol):
    """``scipy.integrate._ivp.common.select_initial_step`` (HNW II.4); one evaluation."""
    scale = atol + np.abs(y0) * rtol
    d0 = _rms(y0 / scale)
    d1 = _rms(f0 / scale)
    h0 = 1e-6 if d0 < 1e-5 or d1 < 1e-5 else 0.01 * d0 / d1
    h0 = min(h0, interval)
    f1 = rhs(y0 + h0 * f0)
    d2 = _rms((f1 - f0) / scale) / h0
    if d1 <= 1e-15 and d2 <= 1e-15:
        h1 = max(1e-6, h0 * 1e-3)
    else:
        h1 = (0.01 / max(d1, d2)) ** (1.0 / (order + 1))
    return min(100 * h0, h1, interval, max_step)


@njit
def _error_norm(K, h, scale, E, E3, dop853):
    n = scale.size
    err5 = np.zeros(n)
    for j in range(E.size):
        if E[j] != 0.0:
            err5 += E[j] * K[j]
    if not dop853:
        return _rms(err5 * h / scale)
    err3 = np.zeros(n)
    for j in range(E3.size):
        if E3[j] != 0.0:
            err3 += E3[j] * K[j]
    err5 /= scale
    err3 /= scale
    e5 = np.sum(err5 * err5)
    e3 = np.sum(err3 * err3)
    if e5 == 0.0 and e3 == 0.0:
        return 0.0
    return abs(h) * e5 / np.sqrt((e5 + 0.01 * e3) * n)


@njit
def adaptive_driver(
    rhs, y0, lam0, lam1, rtol, atol, P, E, E3, dop853, err_order, first_step, max_step,
    ev_index, ev_value, ev_dir, ev_stop, save_every, max_steps,
):  # fmt: skip
    A, B = P[0], P[1]
    dim = y0.size
    n_stages = B.size
    exponent = -1.0 / (err_order + 1)
    cap = 64
    lam_out = np.empty(cap)
    y_out = np.empty((cap, dim))
    lam_out[0] = lam0
    y_out[0] = y0
    n_out = 1
    buf = (np.empty(16), np.empty((16, dim)), np.empty(16, np.int64), 0)
    counts = np.zeros(ev_index.size, np.int64)
    K = np.empty((n_stages + 1, dim))

    y = y0.copy()
    lam = lam0
    f = rhs(y)
    nfev = 1
    if first_step > 0:
        h_abs = first_step
    else:
        h_abs = select_initial_step(rhs, y, f, lam1 - lam0, max_step, err_order, rtol, atol)
        nfev += 1
    n_steps = 0
    n_rejected = 0
    status = 0
    while lam < lam1:
        if n_steps >= max_steps:
            status = -3
            break
        min_step = 10.0 * abs(np.nextafter(lam, np.inf) - lam)
        if h_abs > max_step:
            h_abs = max_step
        elif h_abs < min_step:
            h_abs = min_step
        accepted = False
        rejected = False
        while not accepted:
            if h_abs < min_step:
                status = -2
                break
            lam_new = lam + h_abs
            if lam_new > lam1:
                lam_new = lam1
            h = lam_new - lam
            h_abs = abs(h)
            nfev += _stages(rhs, y, f, h, A, n_stages, K)
            y_new = y + h * _combine(K, B, n_stages)
            f_new = rhs(y_new)
            nfev += 1
            K[n_stages] = f_new
            scale = atol + np.maximum(np.abs(y), np.abs(y_new)) * rtol
            error_norm = _error_norm(K, h, scale, E, E3, dop853)
            if error_norm < 1.0:
                if error_norm == 0.0:
                    factor = MAX_FACTOR
                else:
                    factor = min(MAX_FACTOR, SAFETY * error_norm**exponent)
                if rejected:
                    factor = min(1.0, factor)
                h_abs *= factor
                accepted = True
            else:
                h_abs *= max(MIN_FACTOR, SAFETY * error_norm**exponent)
                rejected = True
                n_rejected += 1
        if status != 0:
            break
        if not all_finite(y_new):
            status = -1
            break
        n_steps += 1
        if ev_index.size > 0:
            P[11][:] = f
            stop, stop_lam, stop_y, fe, buf = process_events(
                EMBEDDED, rhs, y, y_new, lam, h, P, ev_index, ev_value, ev_dir, ev_stop, counts, buf
            )
            nfev += fe
            if stop >= 0:
                lam = stop_lam
                y = stop_y
                status = 1 + stop
                break
        lam = lam_new
        y = y_new
        f = f_new
        if n_steps % save_every == 0 and lam < lam1:
            if n_out == lam_out.size:
                new_lam = np.empty(2 * n_out)
                new_y = np.empty((2 * n_out, dim))
                new_lam[:n_out] = lam_out
                new_y[:n_out] = y_out
                lam_out, y_out = new_lam, new_y
            lam_out[n_out] = lam
            y_out[n_out] = y
            n_out += 1
    if n_out == lam_out.size:
        new_lam = np.empty(n_out + 1)
        new_y = np.empty((n_out + 1, dim))
        new_lam[:n_out] = lam_out
        new_y[:n_out] = y_out
        lam_out, y_out = new_lam, new_y
    lam_out[n_out] = lam
    y_out[n_out] = y
    n_out += 1
    lams, ys, ids, n_found = buf
    return (
        lam_out[:n_out], y_out[:n_out], status, nfev, n_steps, n_rejected,
        ids[:n_found], lams[:n_found], ys[:n_found],
    )  # fmt: skip


# --- dispatch and typed entry points ------------------------------------------------------


@njit
def step(kind, rhs, y, err, h, P, comp):
    """One step of the method ``kind`` (see the constants at the top of the module)."""
    if kind == EXPLICIT:
        return explicit_rk_step(rhs, y, err, h, P, comp)
    if kind == GAUSS:
        return gauss_step(rhs, y, err, h, P, comp)
    if kind == TAO:
        return tao_step(rhs, y, err, h, P, comp)
    return embedded_partial_step(rhs, y, err, h, P, comp)


@njit((i8, RHS, VEC, f8, f8, i8, PARAMS, b1, *EVENTS, i8), cache=True)
def fixed(kind, rhs, y0, lam0, h, n, P, comp, ev_i, ev_v, ev_d, ev_s, save_every):
    """Fixed-step integration with method ``kind``; ``lam = lam0 + n h``."""
    return fixed_step_driver(kind, rhs, y0, lam0, h, n, P, comp, ev_i, ev_v, ev_d, ev_s, save_every)


@njit((RHS, VEC, f8, f8, f8, VEC, PARAMS, VEC, VEC, b1, i8, f8, f8, *EVENTS, i8, i8), cache=True)
def adaptive(
    rhs, y0, lam0, lam1, rtol, atol, P, E, E3, dop853, err_order, first_step, max_step,
    ev_i, ev_v, ev_d, ev_s, save_every, max_steps,
):  # fmt: skip
    """Adaptive integration with SciPy's controller; ``P`` holds ``A``, ``B`` and room for ``f``."""
    return adaptive_driver(
        rhs, y0, lam0, lam1, rtol, atol, P, E, E3, dop853, err_order, first_step, max_step,
        ev_i, ev_v, ev_d, ev_s, save_every, max_steps,
    )  # fmt: skip


def make_params(
    A=None, b=None, A_inv=None, X=None, memory=None, fractions=None, coupling=None,
    scalars=None, cache_size=0, f_size=0,
):  # fmt: skip
    """Build the parameter tuple ``P``; omitted fields become empty arrays of the right type."""

    def mat(a):
        return np.zeros((0, 0)) if a is None else np.ascontiguousarray(a, dtype=np.float64)

    def vec(a):
        return np.zeros(0) if a is None else np.ascontiguousarray(a, dtype=np.float64)

    return (
        mat(A), vec(b), mat(A_inv), mat(X), mat(memory), vec(fractions),
        np.ones(0, np.bool_) if coupling is None else np.ascontiguousarray(coupling, np.bool_),
        vec(scalars), np.zeros(cache_size), np.zeros(cache_size), np.zeros(1, np.bool_),
        np.zeros(f_size),
    )  # fmt: skip
