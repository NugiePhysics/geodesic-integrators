"""The experiment harness: rows, caching, provenance, parallel runs."""

import json

import pytest

from geoint.experiments import runner
from geoint.testcases.schwarzschild import Circular, Eccentric


@pytest.fixture
def cache_dir(tmp_path, monkeypatch):
    monkeypatch.setattr(runner, "RAW", tmp_path)
    return tmp_path


def test_row_contents_and_cache(cache_dir):
    spec = runner.RunSpec.make(Circular(10.0, 1.0), "b", "RK4", n_steps=200)
    row = runner.run(spec)
    assert row["case"] == "TC3-r_c10-n_orbits1" and row["method"] == "RK4"
    assert row["opt_n_steps"] == 200 and row["n_steps"] == 200 and row["nfev"] == 800
    assert row["status"] == "completed" and row["error"] < 1e-6
    for key in ("git_hash", "numpy", "numba", "cpu", "wall_time", "dH_max"):
        assert key in row
    path = cache_dir / f"{spec.key()}.json"
    assert json.loads(path.read_text()) == row
    row2 = runner.run(spec)  # from the cache
    assert row2 == row


def test_key_depends_on_every_part_of_the_spec():
    base = runner.RunSpec.make(Circular(10.0, 1.0), "b", "RK4", n_steps=200)
    others = [
        runner.RunSpec.make(Circular(10.0, 2.0), "b", "RK4", n_steps=200),
        runner.RunSpec.make(Circular(10.0, 1.0), "a", "RK4", n_steps=200),
        runner.RunSpec.make(Circular(10.0, 1.0), "b", "GL2", n_steps=200),
        runner.RunSpec.make(Circular(10.0, 1.0), "b", "RK4", n_steps=201),
    ]
    assert len({base.key(), *(s.key() for s in others)}) == 5


def test_tao_couples_only_the_non_cyclic_coordinates(cache_dir):
    spec = runner.RunSpec.make(Eccentric(20.0, 0.5, 1.0), "b", "Tao4", n_steps=2000, omega=0.01)
    solution, _ = runner.solve(spec)
    assert solution.options["coupling"] == [False, True, True, False]
    z = solution.extended
    assert (z[:, 4] == z[0, 4]).all() and (z[:, 12] == z[0, 4]).all()  # E exact in both copies


def test_run_many_in_parallel(cache_dir):
    specs = [
        runner.RunSpec.make(Circular(10.0, 1.0), form, "RK4", n_steps=n)
        for form in ("a", "b")
        for n in (100, 200)
    ]
    table = runner.run_many(specs, processes=2)
    assert list(table["opt_n_steps"]) == [100, 200, 100, 200]
    assert len(list(cache_dir.glob("*.json"))) == 4
