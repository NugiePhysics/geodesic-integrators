# Reproducibility

## One command

```bash
uv sync --all-extras    # the locked environment (uv.lock)
make all                # every figure and table, then the report
```

`make all` runs `make figures`, `make tables` and `make report`:

- `make figures` runs the scripts in [`experiments/`](../experiments/) in dependency order ([`scripts/make_figures.py`](../scripts/make_figures.py)). Each integration is described by a specification (test case, formulation, method, options) whose hash keys a cache in `results/raw/` (git-ignored). A second run only redraws. `make clean-cache` forces a full recomputation, which takes about 25 minutes on four cores.
- Small summary tables go to [`results/summary/`](../results/summary/) (committed), figures to [`figures/`](../figures/) as PDF and PNG. The PDFs carry no creation date, so a rerun produces byte-identical files.
- `make tables` turns the summaries into LaTeX tables and number macros for the report ([`scripts/make_tables.py`](../scripts/make_tables.py)). No number in the report is typed by hand.
- `make report` compiles [`report/main.tex`](../report/main.tex) and the one-page [`report/summary.tex`](../report/summary.tex) with [Tectonic](https://tectonic-typesetting.github.io/).

Running everything from an empty cache reproduced every summary table and figure byte for byte. The only differences were wall times (median ratio 0.98, 90 % within 0.75–1.15). A scheduled CI job ([`reproduce.yml`](../.github/workflows/reproduce.yml)) runs `make all` from a clean checkout and uploads the figures and the report.

## Provenance

Every cached result row carries the git commit (and whether the tree was dirty), the versions of Python, NumPy, SciPy and Numba, the CPU, and a timestamp. The environment of the published results:

<!-- include: results/summary/t8_environment.md -->

*Table T8.*

## Verification levels

| Level | What is checked | Where |
|---|---|---|
| 1. Units | $`g\,g^{-1} = I`$; symmetry of $`\Gamma`$; $`\Gamma`$ and $`\partial g^{\mu\nu}`$ against a SymPy derivation at 50 digits; (a)–(b) consistency through $`p = g u`$; initial data on the mass shell (property tests) | `tests/unit` |
| 2. Integrators on toy problems | orders, symplecticity $`\Psi'^\mathsf{T}J\Psi' = J`$, symmetry $`\Psi_{-h}\circ\Psi_h = \mathrm{id}`$, bounded energy over $`10^6`$ steps | `tests/integrators` |
| 3. Cross-implementation | our DP5 and DOP853 against SciPy, step for step; compiled against interpreted (`NUMBA_DISABLE_JIT=1`); formulation (a) against (b) | `tests/integrators`, `tests/physics` |
| 4. Physics | every test case against its exact reference | `tests/physics`, table T2 |
| 5. Regression | the validation table of the predecessor project | `tests/regression` |

`make test` runs the fast suite (about 30 s once compiled); `make test-slow` the long integrations; `make test-nojit` the same tests as plain Python.

## Timing

Costs are counted in vector-field evaluations ([ADR 0003](decisions/0003-nfev-as-cost-metric.md)). Wall times compare only our own compiled implementations, with `NUMBA_NUM_THREADS = OMP_NUM_THREADS = 1`, after a warm-up run, as the median of the repeats. Timing runs execute one at a time; other sweeps use up to four worker processes.
