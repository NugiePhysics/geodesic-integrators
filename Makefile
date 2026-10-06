.PHONY: sync lint format test test-slow test-nojit

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
