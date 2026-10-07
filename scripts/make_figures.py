"""Regenerate every figure and summary table from scratch (``make figures``).

Runs the scripts in ``experiments/`` in dependency order. Integrations are cached under
``results/raw/`` (git-ignored), so a second run only redraws; delete that directory, or pass
``--only`` with script names, to recompute selectively.

    uv run python scripts/make_figures.py                 # everything
    uv run python scripts/make_figures.py --only f04 f15  # selected scripts
"""

from __future__ import annotations

import argparse
import importlib
import sys
import time
from pathlib import Path

EXPERIMENTS = Path(__file__).resolve().parents[1] / "experiments"

ORDER = [
    # Phase 5: short integrations
    "f14_tao_omega",  # fixes Tao's omega for everything below
    "f04_convergence",
    "f05_tolerance",
    "f06_work_precision",
    "f02_deflection",
    "f03_deflection_cost",
    "f12_windings",
    "f13_precession",
    "f15_implicit",  # reads the F4 table
    "f16_step_history",
    "t2_validation",
]


def main(argv=None) -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--only", nargs="*", default=None, help="script name prefixes")
    args = parser.parse_args(argv)
    sys.path.insert(0, str(EXPERIMENTS))
    for name in ORDER:
        if args.only and not any(name.startswith(prefix) for prefix in args.only):
            continue
        t0 = time.perf_counter()
        importlib.import_module(name).main()
        print(f"{name}: {time.perf_counter() - t0:.0f} s", flush=True)


if __name__ == "__main__":
    main()
