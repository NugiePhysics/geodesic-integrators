"""T8: the computing environment of the results (reproducibility).

Versions and CPU come from the same provenance record that every cached run carries. Timing
conventions follow ADR 0003: one thread per process, compiled code warmed up, median of
repeats, timing runs executed one at a time.
"""

from __future__ import annotations

import importlib.metadata
import os
import platform

import pandas as pd
from _common import SUMMARY, markdown_table, save_table

from geoint.experiments import provenance

PACKAGES = ("numpy", "scipy", "numba", "llvmlite", "mpmath", "matplotlib", "pandas")


def table() -> pd.DataFrame:
    prov = provenance()
    rows = [
        ("CPU", prov["cpu"]),
        ("logical cores", str(os.cpu_count())),
        ("operating system", platform.platform(terse=True)),
        ("Python", prov["python"]),
        *((name, importlib.metadata.version(name)) for name in PACKAGES),
        ("geoint", prov["geoint"]),
        ("threads per process", "1 (NUMBA_NUM_THREADS = OMP_NUM_THREADS = 1)"),
        ("parallel sweeps", "up to 4 forked worker processes; timing runs one at a time"),
        ("timing", "median of 2 repeats after a warm-up run (ADR 0003)"),
        ("floating point", "IEEE 754 double; compensated summation in fixed-step updates"),
    ]
    return pd.DataFrame(rows, columns=["item", "value"])


def main():
    t = table()
    save_table(t, "t8_environment")
    (SUMMARY / "t8_environment.md").write_text(markdown_table(t) + "\n")
    return t


if __name__ == "__main__":
    print(main().to_string())
