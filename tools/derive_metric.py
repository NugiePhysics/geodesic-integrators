"""Derive everything a metric kernel provides from ``g_{mu nu}`` alone, with SymPy.

From the covariant metric this computes ``g^{mu nu}``, ``d_a g_{mu nu}``, ``d_a g^{mu nu}`` and
``Gamma^mu_{a b}`` in the index conventions of :mod:`geoint.metrics.base`. The tests compare the
hand-written kernels against it; for Kerr it is meant to grow into a code generator.

Usage::

    uv run python tools/derive_metric.py schwarzschild   # print the non-zero components
"""

from __future__ import annotations

import argparse
from collections.abc import Callable
from dataclasses import dataclass

import sympy as sp

N = 4


@dataclass(frozen=True)
class SymbolicMetric:
    """A metric ``g_{mu nu}(x)`` and the numerical values ``geoint`` uses for its parameters."""

    name: str
    coords: tuple[sp.Symbol, ...]
    g: sp.Matrix
    params: dict[sp.Symbol, sp.Expr]


@dataclass(frozen=True)
class DerivedMetric:
    """Derived quantities, each an ``sp.Array`` indexed like the matching kernel."""

    g: sp.Array  # [mu, nu] = g_{mu nu}
    g_inv: sp.Array  # [mu, nu] = g^{mu nu}
    dg: sp.Array  # [a, mu, nu] = d_a g_{mu nu}
    dg_inv: sp.Array  # [a, mu, nu] = d_a g^{mu nu}
    christoffel: sp.Array  # [mu, a, b] = Gamma^mu_{a b}

    def arrays(self) -> dict[str, sp.Array]:
        return {
            "g": self.g,
            "g_inv": self.g_inv,
            "dg": self.dg,
            "dg_inv": self.dg_inv,
            "christoffel": self.christoffel,
        }


def schwarzschild(M: float = 1.0) -> SymbolicMetric:
    """Schwarzschild metric in Schwarzschild coordinates (theory eq. 1)."""
    t, r, theta, phi = sp.symbols("t r theta phi", real=True)
    m = sp.Symbol("M", positive=True)
    f = 1 - 2 * m / r
    g = sp.diag(-f, 1 / f, r**2, r**2 * sp.sin(theta) ** 2)
    # Rational(M) is the exact binary value of the float, so no precision is lost here.
    return SymbolicMetric("schwarzschild", (t, r, theta, phi), g, {m: sp.Rational(M)})


SPACETIMES: dict[str, Callable[..., SymbolicMetric]] = {"schwarzschild": schwarzschild}


def derive(metric: SymbolicMetric) -> DerivedMetric:
    """All derived arrays of ``metric``, simplified."""
    x, g = metric.coords, metric.g
    g_inv = sp.simplify(g.inv())
    dg = [[[sp.diff(g[m, n], x[a]) for n in range(N)] for m in range(N)] for a in range(N)]
    dg_inv = [
        [[sp.simplify(sp.diff(g_inv[m, n], x[a])) for n in range(N)] for m in range(N)]
        for a in range(N)
    ]
    gamma = [
        [
            [
                sp.simplify(
                    sum(g_inv[s, m] * (dg[a][m][b] + dg[b][m][a] - dg[m][a][b]) for m in range(N))
                    / 2
                )
                for b in range(N)
            ]
            for a in range(N)
        ]
        for s in range(N)
    ]
    return DerivedMetric(
        g=sp.Array(g),
        g_inv=sp.Array(g_inv),
        dg=sp.Array(dg),
        dg_inv=sp.Array(dg_inv),
        christoffel=sp.Array(gamma),
    )


def lambdify(
    metric: SymbolicMetric, array: sp.Array, modules: str = "mpmath"
) -> Callable[..., list]:
    """Numerical function ``x -> nested list`` of ``array``, with the metric parameters inserted.

    With ``modules="mpmath"`` the result is exact to the working precision of ``mpmath.mp``,
    which makes it a reference for float64 kernels.
    """
    expr = array.subs(metric.params).tolist()
    return sp.lambdify([metric.coords], expr, modules=modules)


def describe(metric: SymbolicMetric, derived: DerivedMetric) -> str:
    """The non-zero components, one per line; symmetric pairs are listed once."""
    c = metric.coords
    sections = [
        ("g_{mu nu}", lambda m, n: f"g_{{{c[m]} {c[n]}}}", derived.g),
        ("g^{mu nu}", lambda m, n: f"g^{{{c[m]} {c[n]}}}", derived.g_inv),
    ]
    lines = [f"{metric.name}, parameters {metric.params}"]
    for title, label, array in sections:
        lines.append(f"\n{title}:")
        lines += [
            f"  {label(m, n)} = {array[m, n]}"
            for m in range(N)
            for n in range(m, N)
            if array[m, n] != 0
        ]
    lines.append("\nd_a g^{mu nu}:")
    lines += [
        f"  d_{c[a]} g^{{{c[m]} {c[n]}}} = {derived.dg_inv[a, m, n]}"
        for a in range(N)
        for m in range(N)
        for n in range(m, N)
        if derived.dg_inv[a, m, n] != 0
    ]
    lines.append("\nGamma^mu_{a b}:")
    lines += [
        f"  Gamma^{c[m]}_{{{c[a]} {c[b]}}} = {derived.christoffel[m, a, b]}"
        for m in range(N)
        for a in range(N)
        for b in range(a, N)
        if derived.christoffel[m, a, b] != 0
    ]
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("spacetime", choices=sorted(SPACETIMES))
    args = parser.parse_args(argv)
    metric = SPACETIMES[args.spacetime]()
    print(describe(metric, derive(metric)))


if __name__ == "__main__":
    main()
