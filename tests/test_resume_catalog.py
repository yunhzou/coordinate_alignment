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
        "--output-dir", str(output), "--index", "0", "--workers", "48", "--worker"])
    runner.main()
    assert calls == [("CO", {"branch_limit": 100}, None, True, None), ((572, "C", "source"), 47)]
    with gzip.open(output / "parts/part_0.jsonl.gz", "rt") as stream:
        assert json.loads(stream.read()) == {"row_index": 572, "source_id": "source"}
    with pytest.raises(FileExistsError):
        runner.main()
    assert len(calls) == 2


def test_watchdog_kills_entire_session_and_records_timeout(tmp_path, monkeypatch, capsys):
    spec = importlib.util.spec_from_file_location("resume_watchdog",
        Path(__file__).parents[1] / "bench/resume_catalog_source.py")
    runner = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(runner)
    calls = []
    class Process:
        pid = 12345
        def __enter__(self): return self
        def __exit__(self, *args): pass
        def wait(self, timeout=None):
            calls.append(timeout)
            if timeout is not None:
                raise runner.subprocess.TimeoutExpired("detector", timeout)
            return -9
    def launch(command, *, start_new_session):
        assert start_new_session
        assert command[-1] == "--worker"
        return Process()
    killed = []
    monkeypatch.setattr(runner.subprocess, "Popen", launch)
    monkeypatch.setattr(runner.os, "killpg", lambda pid, sig: killed.append((pid, sig)))
    with pytest.raises(SystemExit) as error:
        runner.supervise([], tmp_path, 2)
    assert error.value.code == 124
    assert calls == [300, 300, None]
    assert killed == [(12345, runner.signal.SIGKILL)]
    assert json.loads((tmp_path / "parts/part_2.watchdog.json").read_text())["result_complete"] is False
    assert json.loads(capsys.readouterr().out)["event"] == "slow_run"
