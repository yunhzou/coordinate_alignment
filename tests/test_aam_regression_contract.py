import json
from pathlib import Path
import runpy
from types import SimpleNamespace

from rxn_core.benchmark_regression import evaluate_record, load_contract


CONTRACT = Path(__file__).parents[1] / "benchmarks" / "aam_regression_contract.json"
BENCHMARK_TOOL = Path(__file__).parents[1] / "tools" / "benchmark_aam_versions.py"
CASE = "pr15.Fe_crosscoupling_ACScat2023_TS8_step2_reductive_elimination_mult2"


def _healthy_record():
    return {
        "case": CASE,
        "atom_count": 82,
        "selected_mechanism_count": 1,
        "aam_seconds": 300.0,
        "post_aam_seconds": 500.0,
        "total_seconds": 900.0,
        "peak_process_tree_rss_kb": 2_000_000,
        "mechanisms": [{
            "broken_bonds_R": [[15, 12], [37, 12]],
            "formed_bonds_R": [[12, 10], [12, 11], [37, 15]],
            "event_certificate_digest": (
                "68a5b155eac53fef4cf6c614aa84d8879a3a9d6f0937ccbaa2ffb430faf6e792"),
            "index_chirality_violations": 0,
            "degeneracy_groups": [],
        }],
        "regression_metrics": {
            "configured_max_branches": 100,
            "max_live_branches": 100,
            "post_aam_call_counts": {
                "complete_chosen_automorphism_groups": 1,
            },
            "pipeline_post_aam_metrics": {
                "completed_candidate_group_requests": 12,
                "completed_candidate_group_calculations": 7,
                "completed_candidate_group_cache_hits": 5,
            },
        },
    }


def test_contract_is_versioned_and_healthy_record_passes():
    contract = load_contract(CONTRACT)
    assert contract["schema_version"] == "rxn_core.aam_regression_contract/v1"
    report = evaluate_record(
        _healthy_record(), contract, performance_profile="cpunodes_16")
    assert report["passed"] is True
    assert report["failures"] == []
    assert report["advisories"] == []


def test_contract_detects_chemistry_branch_memory_and_time_regressions():
    record = _healthy_record()
    record["mechanisms"][0]["formed_bonds_R"] = [[0, 1]]
    record["mechanisms"][0]["event_certificate_digest"] = "wrong"
    record["mechanisms"][0]["index_chirality_violations"] = 1
    record["regression_metrics"]["max_live_branches"] = 101
    record["total_seconds"] = 1101
    record["peak_process_tree_rss_kb"] = 110_000_000

    report = evaluate_record(
        record, load_contract(CONTRACT), performance_profile="cpunodes_16")
    text = "\n".join(report["failures"])
    assert "mechanism event certificates" in text
    assert "chirality violations" in text
    assert "max_live_branches exceeded" in text
    assert "total_seconds exceeded" in text
    assert "peak_process_tree_rss_kb exceeded" in text


def test_contract_detects_missing_confirmed_degeneracy_group():
    contract = load_contract(CONTRACT)
    contract["cases"][CASE]["required_degeneracy_groups"] = [{
        "center_R": 9,
        "r_atoms": [7, 10],
    }]
    report = evaluate_record(_healthy_record(), contract)
    assert report["passed"] is False
    assert "required degeneracy groups are missing" in report["failures"][0]


def test_duplicate_degeneracy_work_is_a_hard_regression():
    record = _healthy_record()
    record["regression_metrics"]["post_aam_call_counts"][
        "complete_chosen_automorphism_groups"] = 2

    report = evaluate_record(record, load_contract(CONTRACT))
    assert report["passed"] is False
    assert report["failures"] == [
        "complete_chosen_automorphism_groups call count exceeded: 2>1"]


def test_contract_json_contains_no_environment_specific_paths():
    raw = json.loads(CONTRACT.read_text(encoding="utf-8"))
    assert "/h/" not in json.dumps(raw)


def test_benchmark_instrumentation_counts_calls_and_pool_structure():
    module = runpy.run_path(str(BENCHMARK_TOOL))
    ledger = module["_CallLedger"]()
    target = SimpleNamespace(work=lambda value: value + 1)
    original = target.work
    ledger.instrument(target, "work")
    assert target.work(1) == 2
    assert target.work(3) == 4
    ledger.restore()
    assert ledger.calls == {"work": 2}
    assert ledger.seconds["work"] >= 0.0
    assert target.work is original

    metrics = module["_pool_metrics"]({
        "event": {
            "branches": [{
                "branch_symmetry": {
                    "fragments": [{"symmetry": {"blocks": [{}, {}]}}]
                }
            }]
        }
    })
    assert metrics == {
        "pool_branch_count": 1,
        "max_branches_per_mechanism": 1,
        "pool_fragment_record_count": 1,
        "pool_symmetry_block_count": 2,
    }
