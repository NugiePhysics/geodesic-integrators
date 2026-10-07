"""Helpers shared by the figure scripts in this directory."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from geoint.plotting import style

ROOT = Path(__file__).resolve().parents[1]
SUMMARY = ROOT / "results" / "summary"


def save_table(df: pd.DataFrame, name: str, columns=None) -> Path:
    """Write a small CSV to ``results/summary/`` (committed; the raw cache is not)."""
    SUMMARY.mkdir(parents=True, exist_ok=True)
    path = SUMMARY / f"{name}.csv"
    (df[columns] if columns else df).to_csv(path, index=False, float_format="%.6g")
    return path


def markdown_table(df: pd.DataFrame, floatfmt: str = ".3g") -> str:
    """A GitHub-flavoured Markdown table (no extra dependency)."""

    def fmt(v):
        if isinstance(v, float):
            return format(v, floatfmt) if pd.notna(v) else "–"
        return str(v)

    head = "| " + " | ".join(map(str, df.columns)) + " |"
    rule = "|" + "|".join("---" for _ in df.columns) + "|"
    rows = ["| " + " | ".join(fmt(v) for v in row) + " |" for row in df.itertuples(index=False)]
    return "\n".join([head, rule, *rows])


style.apply()
