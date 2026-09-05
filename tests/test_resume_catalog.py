import gzip
import importlib.util
import json
from pathlib import Path
import sys

import pytest


def test_resume_preserves_source_indices_and_config_and_does_not_repeat_saved_work(tmp_path, monkeypatch):
    spec = importlib.util.spec_from_file_location("resume_catalog_source",
        Path(__file__).parents[1] / "bench/resume_catalog_source.py")
    runner = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(runner)
    prior, output = tmp_path / "prior", tmp_path / "resumed"
    (prior / "parts").mkdir(parents=True)
    (prior / "progress_audit.json").write_text(json.dumps({"unfinished_sources": [
        {"source_id": "source", "smiles": "C", "row_index": 572}]}))
    config = {"target_smiles": "CO", "config": {"branch_limit": 100},
              "minimum_target_coverage_fraction": None}
    (prior / "parts/part_0.summary.json").write_text(json.dumps(config))
    calls = []
    monkeypatch.setattr(runner, "_worker_init", lambda *args: calls.append(args))
    def search(row, *, seed_workers):
        calls.append((row, seed_workers))
        return {"rows": 1}, {"row_index": row[0], "source_id": row[2]}
    monkeypatch.setattr(runner, "_search_one", search)
    monkeypatch.setattr(sys, "argv", ["resume", "--prior-run", str(prior),
        "--output-dir", str(output), "--index", "0", "--workers", "48"])
    runner.main()
    assert calls == [("CO", {"branch_limit": 100}, None, True, None), ((572, "C", "source"), 47)]
    with gzip.open(output / "parts/part_0.jsonl.gz", "rt") as stream:
        assert json.loads(stream.read()) == {"row_index": 572, "source_id": "source"}
    with pytest.raises(FileExistsError):
        runner.main()
    assert len(calls) == 2
