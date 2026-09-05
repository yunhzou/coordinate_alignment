import importlib.util
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import threading
import gzip
import io
import json


spec = importlib.util.spec_from_file_location("search_catalog", Path(__file__).parents[1] / "tools/search_mcule_retro.py")
catalog = importlib.util.module_from_spec(spec)
spec.loader.exec_module(catalog)


def test_budgeted_scheduler_reuses_slots_before_slow_job_finishes():
    slow_started, small_finished = threading.Event(), threading.Event()
    lock = threading.Lock()
    active = 0
    seen = []

    def work(row, budget):
        nonlocal active
        with lock:
            active += budget
            assert active <= 4
            seen.append(row)
        if row == "slow":
            slow_started.set()
            assert small_finished.wait(5), "phase barrier blocked a runnable small job"
        elif row == "fast":
            assert slow_started.wait(5)
        else:
            small_finished.set()
        with lock:
            active -= budget
        return row

    with ThreadPoolExecutor(4) as executor:
        results = list(catalog._budgeted_results(executor,
            [("slow", 3), ("fast", 1), ("small", 1)], 4, work))
    assert sorted(results) == ["fast", "slow", "small"]
    assert sorted(seen) == sorted(results)


def test_outlier_budgets_fit_allocation():
    rows = [(0, "CCCCCCCC", "a"), (1, "CCC", "b"), (2, "C", "c")]
    ordinary, outliers, _ = catalog._adaptive_partition(rows, 4)
    budgets = catalog._outlier_worker_budgets(outliers, 4)
    assert sum(budgets) <= 4
    assert all(b >= 1 for b in budgets)
    assert sorted(ordinary + outliers) == rows


def test_worker_gzip_members_are_one_complete_jsonl_stream():
    records = [{"id": i, "atoms": list(range(i)), "name": "α"} for i in range(4)]
    data = catalog._encode_records(()) + b"".join(catalog._encode_records([r]) for r in records)
    with gzip.open(io.BytesIO(data), "rt") as stream:
        assert [json.loads(line) for line in stream] == records


def test_progress_audit_counts_flushed_rows_without_claiming_complete_scan(tmp_path):
    audit_spec = importlib.util.spec_from_file_location("catalog_progress",
        Path(__file__).parents[1] / "bench/catalog_scan_progress.py")
    progress = importlib.util.module_from_spec(audit_spec)
    audit_spec.loader.exec_module(progress)
    (tmp_path / "logs").mkdir()
    (tmp_path / "parts").mkdir()
    inventory = tmp_path / "bank.csv.gz"
    with gzip.open(inventory, "wt") as stream:
        stream.write("SMILES,Inventory ID\nC,a\nO,b\n")
    row = {"precursor_id": "a", "pair_elapsed_seconds": 1.0}
    (tmp_path / "logs/1_0.out").write_text(json.dumps(row) + '\n{\n  "counts": {}\n}\n')
    report = progress.audit(tmp_path, inventory, "Inventory ID")
    assert report["saved_rows"] == 1
    assert not report["scan_complete"]
    assert report["unfinished_sources"] == [{"source_id": "b", "smiles": "O", "row_index": 1}]
