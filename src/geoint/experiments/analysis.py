"""Fits and interpolations shared by the figure and table scripts."""

from __future__ import annotations

import numpy as np


def fit_slope(x, y, *, lo: float = 1e-12, hi: float = 1e-3):
    """Least-squares slope of ``log y`` against ``log x`` over ``lo < y < hi``.

    The window drops the round-off floor (below) and the pre-asymptotic regime (above).
    Returns ``(slope, standard_error, n_points)``.
    """
    x, y = np.asarray(x, float), np.asarray(y, float)
    keep = np.isfinite(x) & np.isfinite(y) & (y > lo) & (y < hi)
    if keep.sum() < 3:
        return np.nan, np.nan, int(keep.sum())
    X, Y = np.log(x[keep]), np.log(y[keep])
    A = np.vstack([X, np.ones_like(X)]).T
    coef, residual, *_ = np.linalg.lstsq(A, Y, rcond=None)
    n = len(X)
    sigma2 = float(residual[0]) / (n - 2) if residual.size and n > 2 else 0.0
    se = np.sqrt(sigma2 / np.sum((X - X.mean()) ** 2))
    return float(coef[0]), float(se), n


def pareto_front(cost, error):
    """Points of a work-precision cloud that no cheaper point beats, sorted by cost."""
    cost, error = np.asarray(cost, float), np.asarray(error, float)
    ok = np.isfinite(cost) & np.isfinite(error) & (error > 0)
    order = np.argsort(cost[ok])
    c, e = cost[ok][order], error[ok][order]
    keep = e < np.minimum.accumulate(np.concatenate([[np.inf], e[:-1]]))
    return c[keep], e[keep]


def cost_at_error(cost, error, target: float) -> float:
    """Cost at which the work-precision front first reaches ``target`` (log-log interpolation)."""
    c, e = pareto_front(cost, error)
    if c.size == 0 or e.min() > target:
        return np.nan
    if e[0] <= target:
        return float(c[0])
    k = int(np.argmax(e <= target))
    t = (np.log(target) - np.log(e[k - 1])) / (np.log(e[k]) - np.log(e[k - 1]))
    return float(np.exp(np.log(c[k - 1]) + t * (np.log(c[k]) - np.log(c[k - 1]))))


def error_at_cost(cost, error, budget: float) -> float:
    """Error of the work-precision front at cost ``budget`` (log-log interpolation)."""
    c, e = pareto_front(cost, error)
    if c.size == 0 or budget < c[0]:
        return np.nan
    if budget >= c[-1]:
        return float(e[-1])
    k = int(np.searchsorted(c, budget))
    t = (np.log(budget) - np.log(c[k - 1])) / (np.log(c[k]) - np.log(c[k - 1]))
    return float(np.exp(np.log(e[k - 1]) + t * (np.log(e[k]) - np.log(e[k - 1]))))


def envelope(x, y, width: float):
    """Maximum of ``y`` in consecutive windows of ``x`` of the given width (e.g. one orbit).

    A trailing window with fewer than half the typical number of samples (typically the single
    end point) is dropped, since its maximum is not comparable with the others.
    """
    x, y = np.asarray(x, float), np.asarray(y, float)
    bins = np.floor((x - x[0]) / width).astype(int)
    edges = np.flatnonzero(np.diff(bins)) + 1
    groups = np.split(np.arange(x.size), edges)
    if len(groups) > 2 and len(groups[-1]) < 0.5 * np.median([len(g) for g in groups]):
        groups = groups[:-1]
    return np.array([x[g[-1]] for g in groups]), np.array([np.max(y[g]) for g in groups])


GROWTH_CLASSES = {0.0: "bounded", 0.5: "random walk", 1.0: "linear", 2.0: "quadratic"}


def classify_growth(exponent: float) -> str:
    """Nearest of the laws 0, 1/2, 1, 2 (roadmap §1.4)."""
    if not np.isfinite(exponent):
        return "n/a"
    nearest = min(GROWTH_CLASSES, key=lambda p: abs(p - exponent))
    return GROWTH_CLASSES[nearest]
