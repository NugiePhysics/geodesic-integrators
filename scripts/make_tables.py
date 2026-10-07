"""Turn the summary tables in ``results/summary/`` into LaTeX for the report (``make tables``).

Every table and every number quoted in ``report/`` comes from here: tables go to
``report/tables/*.tex`` and headline numbers to ``report/tables/numbers.tex`` as macros, so
that rerunning the experiments and then this script updates the report without hand edits.

    uv run python scripts/make_tables.py
"""

from __future__ import annotations

import math
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SUMMARY = ROOT / "results" / "summary"
OUT = ROOT / "report" / "tables"
METHODS = ["RK4", "DP5", "DOP853", "GL1", "GL2", "GL3", "Tao2", "Tao4"]
UNICODE = {
    "Δ": r"$\Delta$",
    "φ": r"$\phi$",
    "θ": r"$\theta$",
    "ω": r"$\omega$",
    "°": r"$^\circ$",
    "²": r"$^2$",
}


def tex(text) -> str:
    """Escape plain text for LaTeX and map the few non-ASCII symbols the tables use."""
    text = str(text)
    for char, repl in (("\\", r"\textbackslash{}"), ("&", r"\&"), ("%", r"\%"), ("#", r"\#")):
        text = text.replace(char, repl)
    text = text.replace("_", r"\_")
    for char in "|<>":
        text = text.replace(char, f"${char}$")
    for char, repl in UNICODE.items():
        text = text.replace(char, repl)
    return text


def sci(x, digits: int = 2, math_mode: bool = True) -> str:
    """``1.4\\times10^{3}``; ``--`` for a missing value."""
    if x is None or (isinstance(x, float) and not math.isfinite(x)) or pd.isna(x):
        return "--"
    x = float(x)
    if x == 0:
        return "$0$" if math_mode else "0"
    exponent = math.floor(math.log10(abs(x)))
    mantissa = x / 10**exponent
    if round(mantissa, digits - 1) >= 10:
        mantissa /= 10
        exponent += 1
    if -2 <= exponent <= 2 and digits >= 2:
        body = f"{x:.{max(0, digits - 1 - exponent)}f}"
    else:
        body = f"{mantissa:.{digits - 1}f}\\times10^{{{exponent}}}"
    return f"${body}$" if math_mode else body


def fixed(x, places: int = 2) -> str:
    if pd.isna(x):
        return "--"
    text = f"{float(x):.{places}f}"
    return text[1:] if text.startswith("-") and float(text) == 0 else text


def tabular(columns: str, header: list[str], rows: list[list[str]], groups=None) -> str:
    """A booktabs tabular; ``groups`` maps row indices to midrules inserted before them."""
    groups = set(groups or ())
    lines = [rf"\begin{{tabular}}{{{columns}}}", r"\toprule", " & ".join(header) + r" \\"]
    lines.append(r"\midrule")
    for k, row in enumerate(rows):
        if k in groups:
            lines.append(r"\midrule")
        lines.append(" & ".join(row) + r" \\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    return "\n".join(lines) + "\n"


def read(name: str) -> pd.DataFrame:
    return pd.read_csv(SUMMARY / f"{name}.csv")


def t1() -> str:
    df = read("t1_methods")
    rows = [
        [
            r["method"],
            tex(r["structure"]),
            str(r["order"]),
            r["symplectic"],
            r["symmetric"],
            r["adaptive"],
            tex(r["evals_per_step"]),
            tex(r["formulations"]),
        ]
        for _, r in df.iterrows()
    ]
    header = ["Method", "Structure", "Order", "Sympl.", "Symm.", "Adapt.", "Evaluations per step"]
    return tabular("llcccclc", [*header, "Form."], rows)


def t2() -> str:
    df = read("t2_validation")
    rows = [
        [tex(r["check"]), f"{r['numerical']:.12g}", f"{r['analytical']:.12g}", sci(r["difference"])]
        for _, r in df.iterrows()
    ]
    return tabular("lrrr", ["Check", "Numerical", "Analytical", "Difference"], rows)


def t3() -> str:
    df = read("t3_orders")
    cols = [("deflection", "a"), ("deflection", "b"), ("orbit", "a"), ("orbit", "b")]
    rows = []
    for method in ["RK4", "GL1", "GL2", "GL3", "Tao2", "Tao4"]:
        g = df[df["method"] == method]
        theory = int(g["theory"].iloc[0])
        row = [method, str(theory)]
        for case, form in cols:
            hit = g[(g["case"] == case) & (g["formulation"] == form)]
            if hit.empty:
                row.append("--")
            else:
                r = hit.iloc[0]
                row.append(f"${r['measured']:.2f} \\pm {r['fit_se']:.2f}$")
        rows.append(row)
    header = ["Method", "Theory", "Deflection (a)", "Deflection (b)", "Orbit (a)", "Orbit (b)"]
    return tabular("lccccc", header, rows)


def t4() -> str:
    df = read("t4_growth").set_index(["method", "formulation"])
    rows, groups = [], []
    for form in ("b", "a"):
        groups.append(len(rows))
        for method in METHODS:
            if (method, form) not in df.index:
                continue
            r = df.loc[(method, form)]
            rows.append(
                [
                    f"{method} ({form})",
                    sci(r["dH_last_orbit"]),
                    fixed(r["dH_exponent"]),
                    tex(r["dH_growth"]),
                    sci(r["phase_err_end"]),
                    fixed(r["phase_exponent"]),
                    tex(r["phase_growth"]),
                    sci(r["phase_c1"]),
                    sci(r["phase_c2"]),
                ]
            )
    header = [
        "Method",
        r"$|\delta H|$, last period",
        "Exp.",
        "Law",
        "Phase error",
        "Exp.",
        "Law",
        "$c_1$",
        "$c_2$",
    ]
    return tabular("lrrlrrlrr", header, rows, groups[1:])


def t5() -> str:
    df = read("t5_cost_to_target").set_index(["case", "formulation", "method"])
    cases = ["deflection b = 6", "orbit (20, 0.5), 3 periods", "inclined orbit, 60°, 2 periods"]
    rows = []
    for method in METHODS:
        row = [method]
        for case in cases:
            for target in ("1e-06", "1e-10"):
                cell = []
                for form in ("b", "a"):
                    key = (case, form, method)
                    value = df.loc[key, f"nfev@{target}"] if key in df.index else float("nan")
                    cell.append(sci(value, math_mode=False))
                b, a = cell
                row.append(f"${b}$ ({'$' + a + '$' if a != '--' else a})" if b != "--" else "--")
        rows.append(row)
    header = [
        "Method",
        r"Deflection, $10^{-6}$",
        r"$10^{-10}$",
        r"Orbit, $10^{-6}$",
        r"$10^{-10}$",
        r"Inclined, $10^{-6}$",
        r"$10^{-10}$",
    ]
    return tabular("lrrrrrr", header, rows)


def t5b() -> str:
    df = read("t5b_cost_long").set_index(["formulation", "method"])
    cols = ["phase 0.0001", "phase 1e-08", "dH 1e-10"]
    rows = []
    for method in METHODS:
        row = [method]
        for col in cols:
            b = df.loc[("b", method), col] if ("b", method) in df.index else float("nan")
            a = df.loc[("a", method), col] if ("a", method) in df.index else float("nan")
            row.append(f"{sci(b)} ({sci(a)})" if not pd.isna(a) else sci(b))
        rows.append(row)
    header = ["Method", r"Phase $10^{-4}$", r"Phase $10^{-8}$", r"$|\delta H|$ $10^{-10}$"]
    return tabular("lrrr", header, rows)


def t6() -> str:
    df = read("t6_photon_sphere")
    rows = []
    for _, r in df.iterrows():
        rows.append(
            [
                f"{r['method']} ({r['formulation']})",
                tex(r["setting"]),
                fixed(r["orbits_per_decade"]),
                fixed(r["predicted_per_decade"]),
                fixed(r["orbits_tightest"]),
            ]
        )
    header = ["Method", "Varied", "Orbits per decade", "Predicted", "Orbits, finest"]
    return tabular("llrrr", header, rows)


def t8() -> str:
    df = read("t8_environment")
    rows = [[tex(r["item"]), tex(r["value"])] for _, r in df.iterrows()]
    return tabular("ll", ["Item", "Value"], rows)


def numbers() -> str:
    """Headline numbers quoted in the text, as macros."""
    t4 = read("t4_growth").set_index(["method", "formulation"])
    t5 = read("t5_cost_to_target").set_index(["case", "formulation", "method"])
    t5b = read("t5b_cost_long").set_index(["formulation", "method"])
    t6 = read("t6_photon_sphere").set_index(["method", "formulation", "setting"])
    f11 = read("f11_isco_scaling").set_index("series")
    f19 = read("f19_roundoff").set_index(["quantity", "method", "compensated"])
    ps = read("f19_photon_sphere_summation")
    orbit = "orbit (20, 0.5), 3 periods"

    def ps_orbits(method, comp, spo):
        g = ps[(ps["method"] == method) & (ps["opt_compensated"] == comp)]
        return g.iloc[(g["steps_per_orbit"] - spo).abs().argmin()]["orbits"]

    gl3_sat = t6.loc[("GL3", "b", "step size"), "orbits_tightest"]
    values = {
        "dHRKfourLast": sci(t4.loc[("RK4", "b"), "dH_last_orbit"]),
        "dHDOPLast": sci(t4.loc[("DOP853", "b"), "dH_last_orbit"]),
        "dHGLtwo": sci(t4.loc[("GL2", "b"), "dH_last_orbit"]),
        "dHGLthreeLast": sci(t4.loc[("GL3", "b"), "dH_last_orbit"]),
        "expGLthree": fixed(t4.loc[("GL3", "b"), "dH_exponent"]),
        "phaseDOPend": sci(t4.loc[("DOP853", "b"), "phase_err_end"]),
        "phaseGLthreeEnd": sci(t4.loc[("GL3", "b"), "phase_err_end"]),
        "crossRKfour": f"{t4.loc[('RK4', 'b'), 'phase_crossover']:.0f}",
        "shortDOP": sci(t5.loc[(orbit, "b", "DOP853"), "nfev@1e-10"] / 3),
        "shortGLthree": sci(t5.loc[(orbit, "b", "GL3"), "nfev@1e-10"] / 3),
        "longDOPfour": sci(t5b.loc[("b", "DOP853"), "phase 0.0001"]),
        "longGLthreeFour": sci(t5b.loc[("b", "GL3"), "phase 0.0001"]),
        "longGLthreeEight": sci(t5b.loc[("b", "GL3"), "phase 1e-08"]),
        "psDOPdecade": fixed(t6.loc[("DOP853", "b", "tolerance"), "orbits_per_decade"]),
        "psRKfourDecade": fixed(t6.loc[("RK4", "b", "step size"), "orbits_per_decade"]),
        "psGLthreeDecade": fixed(t6.loc[("GL3", "b", "step size"), "orbits_per_decade"]),
        "psGLthreeSaturation": fixed(gl3_sat, 1),
        "psPlainRKfour": fixed(ps_orbits("RK4", False, 8000), 1),
        "psCompRKfour": fixed(ps_orbits("RK4", True, 8000), 1),
        "iscoDPfive": fixed(f11.loc["DP5 (b)", "slope"], 3),
        "iscoGLtwo": fixed(f11.loc["GL2 (b)", "slope"]),
        "iscoGLthree": fixed(f11.loc["GL3 (b)", "slope"]),
        "circPlain": sci(f19.loc[("circular phase", "RK4", False), "final"]),
        "circComp": sci(f19.loc[("circular phase", "RK4", True), "final"]),
        "circPlainExp": fixed(f19.loc[("circular phase", "RK4", False), "exponent"]),
    }
    lines = ["% Generated by scripts/make_tables.py from results/summary/. Do not edit."]
    lines += [rf"\newcommand{{\num{key}}}{{{value}}}" for key, value in values.items()]
    return "\n".join(lines) + "\n"


TABLES = {"t1": t1, "t2": t2, "t3": t3, "t4": t4, "t5": t5, "t5b": t5b, "t6": t6, "t8": t8}


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    header = "% Generated by scripts/make_tables.py from results/summary/. Do not edit.\n"
    for name, build in TABLES.items():
        (OUT / f"{name}.tex").write_text(header + build())
    (OUT / "numbers.tex").write_text(numbers())
    print(f"wrote {len(TABLES)} tables and numbers.tex to {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
