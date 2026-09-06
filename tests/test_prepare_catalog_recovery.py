"""Restart selection preserves bank identity and never selects saved sources."""
import gzip
import importlib.util
import json
from pathlib import Path
import sys


def test_prepare_only_unfinished_rows_in_selected_shards(tmp_path, monkeypatch, capsys):
    spec = importlib.util.spec_from_file_location(
        "prepare_catalog_recovery",
        Path(__file__).parents[1] / "bench/prepare_catalog_recovery.py",
    )
    runner = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(runner)
    run = tmp_path / "run"
    for directory in ("parts", "logs", "checkpoints"):
        (run / directory).mkdir(parents=True)
    config = {"shard_count": 2, "config": {"iso_tolerance": 1.0, "branch_limit": 100}}
    (run / "parts/part_0.summary.json").write_text(json.dumps(config))
    (run / "logs/42_1.out").write_text('{"precursor_id": "B"}\n{"partial":')
    (run / "checkpoints/3.detection.pkl.gz").touch()
    catalog = tmp_path / "bank.csv.gz"
    with gzip.open(catalog, "wt") as stream:
        stream.write("Bank ID,SMILES\nA,C\nB,N\nC,O\nD,F\n")
    output = tmp_path / "recovery"
    monkeypatch.setattr(sys, "argv", ["prepare", "--run", str(run),
        "--catalog", str(catalog), "--id-column", "Bank ID", "--shards", "1",
        "--output-dir", str(output)])
    runner.main()
    audit = json.loads((output / "progress_audit.json").read_text())
    assert audit["unfinished_sources"] == [{"source_id": "D", "smiles": "F", "row_index": 3}]
    assert json.loads((output / "parts/config.summary.json").read_text()) == config
    assert (output / "checkpoints").resolve() == run / "checkpoints"
    assert json.loads(capsys.readouterr().out)["typed_checkpoints"] == [True]
