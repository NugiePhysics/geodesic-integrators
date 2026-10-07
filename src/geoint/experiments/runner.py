"""Experiment harness: run specification -> integration -> one tidy table row, cached.

A :class:`RunSpec` names a test case, a formulation (``"a"`` second order, ``"b"``
Hamiltonian), an integrator and its options. :func:`run` integrates it and returns a flat dict:
the specification, the case's errors against its exact reference, the costs of ADR 0003
(``nfev``, steps, rejections, implicit iterations, wall time) and the provenance (git hash,
library versions, CPU). Rows are cached as JSON under ``results/raw/`` keyed by a hash of the
specification, so re-running a figure script only integrates what is new. :func:`run_many`
spreads specifications over processes.
"""

from __future__ import annotations

import dataclasses
import functools
import hashlib
import json
import math
import multiprocessing
import os
import platform
import subprocess
import time
from collections.abc import Iterable
from pathlib import Path

import numpy as np

from .. import __version__
from ..formulations import Formulation, Hamiltonian, SecondOrder
from ..integrators import Solution, get
from ..metrics import Schwarzschild
from ..testcases.base import GeodesicCase
from ..testcases.schwarzschild import make_case

#: Bump to invalidate every cached row (e.g. after a change in the error definitions).
CACHE_VERSION = 1
ROOT = Path(__file__).resolve().parents[3]
RAW = ROOT / "results" / "raw"


@functools.cache
def formulation(name: str) -> Formulation:
    metric = Schwarzschild()
    return {"a": SecondOrder, "b": Hamiltonian}[name](metric)


@dataclasses.dataclass(frozen=True)
class RunSpec:
    case: GeodesicCase
    formulation: str
    method: str
    options: tuple = ()

    @classmethod
    def make(cls, case, formulation, method, **options):
        return cls(case, formulation, method, tuple(sorted(options.items())))

    def config(self) -> dict:
        return {
            "case": self.case.id,
            "family": self.case.family,
            "case_params": self.case.params,
            "formulation": self.formulation,
            "method": self.method,
            "options": dict(self.options),
        }

    def key(self) -> str:
        text = json.dumps({**self.config(), "cache": CACHE_VERSION}, sort_keys=True, default=str)
        return hashlib.sha1(text.encode()).hexdigest()[:16]


def solve(spec: RunSpec) -> tuple[Solution, Formulation]:
    """Integrate ``spec`` (uncached) and return the solution and the formulation used."""
    form = formulation(spec.formulation)
    options = dict(spec.options)
    method = get(spec.method)
    if method.requires_hamiltonian and "coupling" not in options:
        # Couple only the non-cyclic coordinates (ADR 0006); integrators know no physics.
        cyclic = form.metric.cyclic_coords
        options["coupling"] = [mu not in cyclic for mu in range(4)]
    return method.solve(spec.case.problem(form), **options), form


@functools.cache
def provenance() -> dict:
    def git(*args):
        try:
            out = subprocess.run(["git", *args], cwd=ROOT, capture_output=True, text=True)
            return out.stdout.strip()
        except OSError:
            return ""

    import numba
    import scipy

    cpu = platform.processor()
    try:
        with open("/proc/cpuinfo") as fh:
            cpu = next(line.split(":", 1)[1].strip() for line in fh if "model name" in line)
    except (OSError, StopIteration):
        pass
    return {
        "git_hash": git("rev-parse", "--short", "HEAD"),
        "git_dirty": bool(git("status", "--porcelain", "--untracked-files=no")),
        "geoint": __version__,
        "python": platform.python_version(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "numba": numba.__version__,
        "cpu": cpu,
    }


def _clean(value):
    if isinstance(value, float | np.floating):
        value = float(value)
        return value if math.isfinite(value) else str(value)
    if isinstance(value, np.integer | np.bool_):
        return value.item()
    return value


def run(spec: RunSpec, *, refresh: bool = False, repeats: int = 1, cache: bool = True) -> dict:
    """One table row for ``spec``; ``repeats > 1`` times it as the median of the repeats."""
    path = RAW / f"{spec.key()}.json"
    if cache and not refresh and path.exists():
        return json.loads(path.read_text())
    walls = []
    for _ in range(max(1, repeats)):
        solution, form = solve(spec)
        walls.append(solution.wall_time)
    if repeats > 1:
        walls = walls[1:]  # the first repetition may include compiling the vector field
    errors = spec.case.errors(solution, form)
    flat_options = {f"opt_{k}": v for k, v in spec.options}
    row = {
        **{k: v for k, v in spec.config().items() if k not in ("options", "case_params")},
        **{f"case_{k}": v for k, v in spec.case.params.items()},
        **flat_options,
        "status": solution.status,
        "nfev": solution.nfev,
        "n_steps": solution.n_steps,
        "n_rejected": solution.n_rejected,
        "n_iter": solution.n_iter,
        "wall_time": float(np.median(walls)),
        "h": solution.options.get("h", math.nan),
        **{k: v for k, v in errors.items() if np.isscalar(v)},
        **provenance(),
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S"),
    }
    row = {k: _clean(v) for k, v in row.items()}
    if cache:
        RAW.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(row, indent=1, default=str))
    return row


def _run_star(args):
    spec, kwargs = args
    return run(spec, **kwargs)


def run_many(specs: Iterable[RunSpec], *, processes: int | None = None, **kwargs):
    """Rows for all ``specs`` as a :class:`pandas.DataFrame`, using up to ``processes`` workers.

    Timing runs (``repeats > 1``) are executed serially, so that they do not compete for cores.
    """
    import pandas as pd

    specs = list(specs)
    if processes is None:
        processes = 1 if kwargs.get("repeats", 1) > 1 else min(4, os.cpu_count() or 1)
    todo = [s for s in specs if kwargs.get("refresh") or not (RAW / f"{s.key()}.json").exists()]
    if processes > 1 and len(todo) > 1:
        # Fork after the compiled code is loaded, so workers do not recompile.
        formulation("a"), formulation("b")
        with multiprocessing.get_context("fork").Pool(processes) as pool:
            pool.map(_run_star, [(s, kwargs) for s in todo], chunksize=1)
    return pd.DataFrame([run(s, **kwargs) for s in specs])


def case_from_row(row) -> GeodesicCase:
    params = {k[len("case_") :]: v for k, v in row.items() if k.startswith("case_")}
    return make_case(row["family"], **params)
