"""Long integrations reduced to compact series, cached as ``.npz``.

A 10^4-orbit run saves 10^5 states or more. Instead of keeping them, :func:`series` reduces a
run to per-period maxima of the invariant errors and the list of periapsis passages, which is
all the growth-law analysis (F7-F9, F17, T4) needs, and caches that under
``results/raw/series/``. :func:`series_many` spreads runs over worker processes.
"""

from __future__ import annotations

import multiprocessing
import os

import numpy as np

from .analysis import envelope
from .runner import RAW, RunSpec, formulation, solve


def _path(spec: RunSpec):
    return RAW / "series" / f"{spec.key()}.npz"


def reduce(spec: RunSpec, period: float) -> dict[str, np.ndarray]:
    """Integrate ``spec`` and reduce the solution; ``period`` sets the envelope windows."""
    sol, form = solve(spec)
    inv = form.invariants(sol.y)
    eps = spec.case.eps
    out = {
        "status": np.array(sol.status),
        "nfev": np.array(sol.nfev),
        "n_steps": np.array(sol.n_steps),
        "n_iter": np.array(sol.n_iter),
        "wall_time": np.array(sol.wall_time),
        "lam_end": np.array(sol.lam_end),
    }
    errors = {
        "dH": form.relative_constraint(sol.y, eps),
        "dE": np.abs(inv["E"] - inv["E"][0]) / abs(inv["E"][0]),
        "dLz": np.abs(inv["Lz"] - inv["Lz"][0]) / max(abs(inv["Lz"][0]), 1e-300),
        "dL2": np.abs(inv["L2"] - inv["L2"][0]) / inv["L2"][0],
    }
    for name, values in errors.items():
        lam, env = envelope(sol.lam, values, period)
        out["lam_env"] = lam
        out[f"{name}_env"] = env
    if sol.extended is not None:
        n = sol.y.shape[1] // 2
        z = sol.extended
        lam, sep = envelope(sol.lam, np.abs(z[:, 1] - z[:, 2 * n + 1]), period)
        out["copy_sep_env"] = sep
    peri = sol.events.get("periapsis")
    if peri is not None:
        out["peri_lam"] = peri.lam
        out["peri_y"] = peri.y
        out["peri_t"] = peri.y[:, 0]
        out["peri_phi"] = peri.y[:, 3]
        out["peri_r"] = peri.y[:, 1]
    stride = max(1, sol.lam.size // 4000)
    out["traj_lam"] = sol.lam[::stride]
    out["traj_y"] = sol.y[::stride]
    out["end_y"] = sol.y_end
    return out


def series(spec: RunSpec, period: float, *, refresh: bool = False) -> dict[str, np.ndarray]:
    path = _path(spec)
    if path.exists() and not refresh:
        with np.load(path, allow_pickle=False) as data:
            return {k: data[k] for k in data.files}
    out = reduce(spec, period)
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(path, **out)
    return out


def _star(args):
    spec, period = args
    series(spec, period)
    return spec.key()


def series_many(specs, period: float, *, processes: int | None = None):
    """Series for every spec (cached), computing missing ones on up to ``processes`` workers."""
    specs = list(specs)
    todo = [s for s in specs if not _path(s).exists()]
    processes = processes or min(4, os.cpu_count() or 1)
    if todo:
        formulation("a"), formulation("b")  # compile before forking
        if processes > 1 and len(todo) > 1:
            with multiprocessing.get_context("fork").Pool(processes) as pool:
                pool.map(_star, [(s, period) for s in todo], chunksize=1)
        else:
            for s in todo:
                series(s, period)
    return [series(s, period) for s in specs]
