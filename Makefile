.PHONY: sync lint format test test-slow test-nojit figures docs docs-serve

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

docs:            ## Documentation site into site/ (fails on broken links)
	uv run --group docs mkdocs build --strict

docs-serve:      ## Documentation site with live reload
	uv run --group docs mkdocs serve
