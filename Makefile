.PHONY: sync lint format test test-slow test-nojit figures tables report all docs docs-serve clean-cache

sync:            ## Install the locked environment (all extras + dev group)
	uv sync --locked --all-extras

lint:            ## Same checks as CI
	uv run ruff check .
	uv run ruff format --check .

format:
	uv run ruff check --fix .
	uv run ruff format .

test:            ## Fast tests (CI)
	uv run pytest

test-slow:       ## Long integrations only
	uv run pytest -m slow

test-nojit:      ## Same tests as pure Python, for debugging
	NUMBA_DISABLE_JIT=1 uv run pytest

figures:         ## Every figure and summary table (cached integrations in results/raw/)
	NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 uv run python scripts/make_figures.py

tables:          ## LaTeX tables and number macros for the report, from results/summary/
	uv run python scripts/make_tables.py

# A fixed date in the PDF metadata makes the report byte-identical between builds.
REPORT_EPOCH := 1791331200

report: tables   ## The report and the one-page summary (needs tectonic)
	cd report && export SOURCE_DATE_EPOCH=$(REPORT_EPOCH) \
		&& tectonic -X compile main.tex && tectonic -X compile summary.tex

all: figures report  ## Every figure, table and the report (integrations are cached in results/raw/)

docs:            ## Documentation site into site/ (fails on broken links)
	uv run --group docs mkdocs build --strict

docs-serve:      ## Documentation site with live reload
	uv run --group docs mkdocs serve

clean-cache:     ## Delete the integration cache, so that the next make all starts from scratch
	rm -rf results/raw
