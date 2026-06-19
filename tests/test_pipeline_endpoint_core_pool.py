import json
import sys

import pytest
import numpy as np
from pathlib import Path

import rxn_core.pipeline as pipeline
from rxn_core.pipeline import (
    _merge_endpoint_core_pools,
    _product_core_pool_to_reactant,
)


def _write_graph_json(path, nodes, edges):
    n = len(nodes)
    weights = np.zeros((n, n))
    for i, j, w in edges:
        weights[i, j] = weights[j, i] = w
    Path(path).write_text(json.dumps({
        "nodes": nodes,
        "weights": weights.tolist(),
    }))


def test_product_core_pool_pulls_back_and_merges_with_reactant_pool():
    core_R = [0, 1]
    mapping_RP = {0: 10, 1: 11}

    r_pool = {
        (((0, 100), (1, 101)), ()): {
            "mapping": {0: 100, 1: 101},
            "dedup_count": 1,
        },
    }
    p_pool_native = {
        (((10, 100), (11, 101)), ()): {
            "mapping": {10: 100, 11: 101},
            "dedup_count": 3,
        },
        (((10, 102), (11, 103)), ()): {
            "mapping": {10: 102, 11: 103},
            "dedup_count": 2,
        },
    }

    p_pool_as_r = _product_core_pool_to_reactant(
        p_pool_native, mapping_RP, core_R)
    merged = _merge_endpoint_core_pools(core_R, r_pool, p_pool_as_r)

    assert len(merged) == 2
    first = merged[(((0, 100), (1, 101)), ())]
    assert first["mapping"] == {0: 100, 1: 101}
    assert first["sources"] == {"R", "P"}
    assert first["dedup_count"] == 4

    second = merged[(((0, 102), (1, 103)), ())]
    assert second["mapping"] == {0: 102, 1: 103}
    assert second["sources"] == {"P"}
    assert second["dedup_count"] == 2


def test_rp_cut_chunks_merge_matches_direct_stage():
    wbo = np.zeros((2, 2))
    wbo[0, 1] = wbo[1, 0] = 0.9
    inputs = pipeline.step_inputs_from_arrays(
        "h2",
        ["H", "H"],
        np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.75]]),
        wbo,
        ["H", "H"],
        np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.75]]),
        wbo,
    )
    cfg = pipeline.rp_stage_config()
    cfg["n_seeds"] = 1

    direct = pipeline.run_rp_stage(inputs, config=cfg, inner_workers=1)
    chunks = [
        pipeline.run_rp_cut_chunk(inputs, [cut], config=cfg)
        for cut in pipeline.rp_cut_work_items(inputs, config=cfg)
    ]
    merged = pipeline.merge_rp_cut_chunks(inputs, chunks, config=cfg)

    assert len(merged["mechanisms"]) == len(direct["mechanisms"])
    assert merged["mechanisms"][0]["mapping_RP"] == direct["mechanisms"][0]["mapping_RP"]


def test_merge_ts_stage_chunks_recomputes_global_top_flags():
    chunk_a = {
        "stage": "ts",
        "step": "toy",
        "config": {},
        "mechanisms": [
            {"id": 1, "gt": None, "igs": [
                {"label": "a", "S": 1.0},
                {"label": "b", "S": 5.0},
            ]},
        ],
        "endpoint_results": [{"target_label": "a"}],
    }
    chunk_b = {
        "stage": "ts",
        "step": "toy",
        "config": {},
        "mechanisms": [
            {"id": 1, "gt": None, "igs": [
                {"label": "c", "S": 4.0},
            ]},
        ],
        "endpoint_results": [{"target_label": "c"}],
    }

    merged = pipeline.merge_ts_stage_chunks([chunk_a, chunk_b])

    igs = merged["mechanisms"][0]["igs"]
    assert [ig["label"] for ig in igs] == ["a", "b", "c"]
    assert [ig["is_top2"] for ig in igs] == [False, True, True]
    assert merged["union_top_labels"] == ["b", "c"]
    assert len(merged["endpoint_results"]) == 2


def test_best_under_mech_scores_all_endpoint_allowed_core_maps(monkeypatch):
    rt_pool = {
        "r_only": {
            "mapping": {0: 0, 1: 1},
            "sources": {"R"},
            "dedup_count": 1,
        },
        "both": {
            "mapping": {0: 1, 1: 0},
            "sources": {"R", "P"},
            "dedup_count": 2,
        },
    }

    def fake_score_one(*args, **kwargs):
        mapping = args[4]
        if mapping == {0: 0, 1: 1}:
            return {"S": 10.0, "core_map": {"0": 0, "1": 1}}
        return {"S": 5.0, "core_map": {"0": 1, "1": 0}}

    monkeypatch.setattr(pipeline, "score_one", fake_score_one)

    best = pipeline.best_under_mech_using_pool(
        ["H", "H"],
        np.zeros((2, 3)),
        np.zeros((2, 2)),
        np.zeros((2, 2)),
        ["H", "H"],
        np.zeros((2, 3)),
        np.zeros((2, 2)),
        np.zeros(1),
        np.zeros((1, 2, 3)),
        rt_pool,
        {0: 0, 1: 1},
        [],
        [],
        [0, 1],
        prefer_endpoint_consensus=False,
    )

    assert best["core_map"] == {"0": 0, "1": 1}
    assert best["core_sources"] == ["R"]
    assert best["endpoint_consensus"] is False


def test_best_under_mech_can_prefer_endpoint_consensus(monkeypatch):
    rt_pool = {
        "r_only": {
            "mapping": {0: 0, 1: 1},
            "sources": {"R"},
            "dedup_count": 1,
        },
        "both": {
            "mapping": {0: 1, 1: 0},
            "sources": {"R", "P"},
            "dedup_count": 2,
        },
    }

    def fake_score_one(*args, **kwargs):
        mapping = args[4]
        if mapping == {0: 0, 1: 1}:
            return {"S": 10.0, "core_map": {"0": 0, "1": 1}}
        return {"S": 5.0, "core_map": {"0": 1, "1": 0}}

    monkeypatch.setattr(pipeline, "score_one", fake_score_one)

    best = pipeline.best_under_mech_using_pool(
        ["H", "H"],
        np.zeros((2, 3)),
        np.zeros((2, 2)),
        np.zeros((2, 2)),
        ["H", "H"],
        np.zeros((2, 3)),
        np.zeros((2, 2)),
        np.zeros(1),
        np.zeros((1, 2, 3)),
        rt_pool,
        {0: 0, 1: 1},
        [],
        [],
        [0, 1],
        prefer_endpoint_consensus=True,
    )

    assert best["core_map"] == {"0": 1, "1": 0}
    assert best["core_sources"] == ["P", "R"]
    assert best["endpoint_consensus"] is True


def test_branch_pool_variants_use_strict_ts_automorphisms(monkeypatch):
    branch_pool = {
        "state": {
            "sources": {"R", "P"},
            "dedup_count": 1,
            "records": [
                {
                    "mapping": {0: 0, 1: 1},
                    "blocks": [{"r_atoms": [0, 1], "p_atoms": [0, 1, 2]}],
                    "source": "R",
                    "dedup_count": 1,
                },
                {
                    "mapping": {0: 0, 1: 1},
                    "blocks": [{"r_atoms": [0, 1], "p_atoms": [0, 1, 2]}],
                    "source": "P",
                    "dedup_count": 1,
                },
            ],
        },
    }
    wbo_r = np.zeros((2, 2))
    wbo_r[0, 1] = wbo_r[1, 0] = 1.0
    wbo_p = np.zeros((2, 2))
    wbo_t = np.zeros((3, 3))
    wbo_t[0, 1] = wbo_t[1, 0] = 1.0

    scored = []

    def fake_score_one(*args, **kwargs):
        mapping = args[4]
        scored.append(dict(mapping))
        # This would win if raw block permutation were still allowed.
        if mapping == {0: 0, 1: 2}:
            return {"S": 100.0, "core_map": {"0": 0, "1": 2}}
        return {
            "S": 10.0 if mapping == {0: 0, 1: 1} else 1.0,
            "core_map": {str(k): v for k, v in mapping.items()},
        }

    monkeypatch.setattr(pipeline, "score_one", fake_score_one)

    best = pipeline.best_under_mech_using_branch_pool(
        ["C", "H"],
        np.zeros((2, 3)),
        wbo_r,
        wbo_p,
        ["C", "H", "H"],
        np.zeros((3, 3)),
        wbo_t,
        np.zeros(1),
        np.zeros((1, 3, 3)),
        branch_pool,
        {0: 0, 1: 1},
        [],
        [],
        [0, 1],
        prefer_endpoint_consensus=True,
        symmetry_wbo_tol=0.2,
    )

    assert {0: 0, 1: 2} not in scored
    assert best["core_map"] == {"0": 0, "1": 1}
    assert best["endpoint_consensus"] is True


def test_sp_cache_cache_only_rejects_missing_wbo(tmp_path):
    cache = tmp_path / "R"
    cache.mkdir()
    (cache / "r.xyz").write_text("1\nR\nH 0 0 0\n")

    with pytest.raises(RuntimeError, match="cache-only"):
        pipeline._ensure_sp_cache(cache, "R", xtb_mode="cache-only")


def test_sp_cache_auto_runs_xtb_for_missing_wbo(tmp_path, monkeypatch):
    cache = tmp_path / "R"
    cache.mkdir()
    (cache / "r.xyz").write_text("2\nR\nH 0 0 0\nH 0 0 1\n")
    calls = []

    def fake_run_xtb(xyz_path, workdir, charge=0, uhf=0, omp_threads=1):
        calls.append((xyz_path, workdir, charge, uhf, omp_threads))
        (workdir / "wbo").write_text("1 2 0.900\n")

    monkeypatch.setattr(pipeline, "_xtb_available", lambda: True)
    monkeypatch.setattr(pipeline, "run_xtb", fake_run_xtb)

    pipeline._ensure_sp_cache(cache, "R", xtb_mode="auto")

    assert calls == [(cache / "r.xyz", cache, 0, 0, 8)]
    assert (cache / "wbo").exists()


def test_sp_cache_auto_uses_charge_and_multiplicity(tmp_path, monkeypatch):
    cache = tmp_path / "R"
    cache.mkdir()
    (cache / "r.xyz").write_text("1\nR\nH 0 0 0\n")
    calls = []

    def fake_run_xtb(xyz_path, workdir, charge=0, uhf=0, omp_threads=1):
        calls.append((charge, uhf))
        (workdir / "wbo").write_text("")

    monkeypatch.setattr(pipeline, "_xtb_available", lambda: True)
    monkeypatch.setattr(pipeline, "run_xtb", fake_run_xtb)

    pipeline._ensure_sp_cache(
        cache, "R", xtb_mode="auto", charge=-1, multiplicity=4)

    assert calls == [(-1, 3)]


def test_hess_cache_auto_uses_charge_and_multiplicity(tmp_path, monkeypatch):
    hess = tmp_path / "hess_iter1"
    hess.mkdir()
    (hess / "ts.xyz").write_text("1\nTS\nH 0 0 0\n")
    calls = []

    def fake_run_xtb_hess(xyz_path, workdir, charge=0, uhf=0, omp_threads=1):
        calls.append((charge, uhf))
        (workdir / "g98.out").write_text("fake")
        (workdir / "wbo").write_text("")

    monkeypatch.setattr(pipeline, "_xtb_available", lambda: True)
    monkeypatch.setattr(pipeline, "run_xtb_hess", fake_run_xtb_hess)
    monkeypatch.setattr(pipeline, "parse_g98_modes", lambda path: ([], []))

    pipeline._ensure_hess_cache(
        hess, "iter1", xtb_mode="auto", charge=2, multiplicity=2)

    assert calls == [(2, 1)]


def test_xtb_worker_auto_respects_thread_budget(monkeypatch):
    monkeypatch.setenv("TSDISCO_SLURM_CPUS", "48")
    monkeypatch.setattr(pipeline, "XTB_OMP_THREADS", "4")
    monkeypatch.setattr(pipeline, "XTB_MAX_THREADS", 4)
    monkeypatch.setattr(pipeline, "XTB_WORKERS", "auto")

    assert pipeline._resolve_xtb_workers(inner_workers=24) == 12


def test_xtb_worker_explicit_override(monkeypatch):
    monkeypatch.setenv("TSDISCO_SLURM_CPUS", "48")
    monkeypatch.setattr(pipeline, "XTB_OMP_THREADS", "4")
    monkeypatch.setattr(pipeline, "XTB_MAX_THREADS", 4)

    assert pipeline._resolve_xtb_workers(inner_workers=24, workers="8") == 8


def test_run_subgraph_cli_matches_by_feature_and_anchor(tmp_path):
    query = tmp_path / "query.json"
    target = tmp_path / "target.json"
    _write_graph_json(
        query,
        [
            {"features": {"outer_shell": 6}},
            {"features": {"outer_shell": 6}},
        ],
        [(0, 1, 1.0)],
    )
    _write_graph_json(
        target,
        [
            {"features": {"outer_shell": 6}},
            {"features": {"outer_shell": 6}},
            {"features": {"outer_shell": 6}},
            {"features": {"outer_shell": 6}},
        ],
        [(0, 1, 1.0), (2, 3, 1.0)],
    )

    result = pipeline.run_subgraph_cli(
        query,
        target,
        node_policy_fields=["outer_shell"],
        anchor_values=["0:2"],
        iso_tol=0.1,
        orbit_dedup=False,
    )

    assert result["n_matches"] == 1
    assert result["matches"][0]["mapping"] == {"0": 2, "1": 3}


def test_main_subgraph_mode_writes_output_json(tmp_path, monkeypatch):
    query = tmp_path / "query.json"
    target = tmp_path / "target.json"
    out = tmp_path / "subgraph_out.json"
    _write_graph_json(
        query,
        [
            {"features": {"outer_shell": 6}},
            {"features": {"outer_shell": 6}},
        ],
        [(0, 1, 1.0)],
    )
    _write_graph_json(
        target,
        [
            {"features": {"outer_shell": 6}},
            {"features": {"outer_shell": 6}},
        ],
        [(0, 1, 1.0)],
    )
    monkeypatch.setattr(sys, "argv", [
        "rxn-core",
        "--subgraph-query-json", str(query),
        "--subgraph-target-json", str(target),
        "--subgraph-node-policy", "outer_shell",
        "--subgraph-anchor", "0:0",
        "--subgraph-output", str(out),
        "--iso-tol", "0.1",
    ])

    pipeline.main()

    data = json.loads(out.read_text())
    assert data["mode"] == "subgraph"
    assert data["n_matches"] == 1
    assert data["matches"][0]["mapping"] == {"0": 0, "1": 1}


def test_alignment_inputs_from_xyz_uses_explicit_cache_dirs(tmp_path, monkeypatch):
    r_xyz = tmp_path / "reactant.xyz"
    p_xyz = tmp_path / "product.xyz"
    r_xyz.write_text("1\nR\nH 0 0 0\n")
    p_xyz.write_text("1\nP\nH 1 0 0\n")
    r_cache = tmp_path / "cache_a"
    p_cache = tmp_path / "cache_b"
    calls = []

    def fake_run_xtb(xyz_path, workdir, charge=0, uhf=0, omp_threads=1):
        calls.append((Path(xyz_path).name, workdir.name, charge, uhf))
        (workdir / "wbo").write_text("")

    monkeypatch.setattr(pipeline, "_xtb_available", lambda: True)
    monkeypatch.setattr(pipeline, "run_xtb", fake_run_xtb)

    inputs = pipeline.alignment_inputs_from_xyz(
        r_xyz, p_xyz, name="direct",
        reactant_workdir=r_cache, product_workdir=p_cache,
        charge=-1, multiplicity=3, xtb_mode="auto")

    assert inputs.step_name == "direct"
    assert inputs.elR == ["H"]
    assert inputs.elP == ["H"]
    assert calls == [
        ("reactant.xyz", "cache_a", -1, 2),
        ("product.xyz", "cache_b", -1, 2),
    ]


def test_discover_mechanisms_from_xyz_can_return_inputs(tmp_path, monkeypatch):
    r_xyz = tmp_path / "R.xyz"
    p_xyz = tmp_path / "P.xyz"
    r_xyz.write_text("2\nR\nH 0 0 0\nH 0 0 0.75\n")
    p_xyz.write_text("2\nP\nH 1 0 0\nH 1 0 0.75\n")

    def fake_run_xtb(_xyz_path, workdir, charge=0, uhf=0, omp_threads=1):
        (workdir / "wbo").write_text("1 2 0.900\n")

    monkeypatch.setattr(pipeline, "_xtb_available", lambda: True)
    monkeypatch.setattr(pipeline, "run_xtb", fake_run_xtb)
    cfg = pipeline.rp_stage_config()
    cfg["n_seeds"] = 1

    inputs, result = pipeline.discover_mechanisms_from_xyz(
        r_xyz, p_xyz, workdir=tmp_path / "cache", name="h2",
        xtb_mode="auto", config=cfg, return_inputs=True)

    assert inputs.step_name == "h2"
    assert result["stage"] == "rp"
    assert result["mechanisms"]


def test_process_xyz_stage_runs_rp_without_step_schema(tmp_path, monkeypatch):
    r_xyz = tmp_path / "R.xyz"
    p_xyz = tmp_path / "P.xyz"
    r_xyz.write_text("2\nR\nH 0 0 0\nH 0 0 0.75\n")
    p_xyz.write_text("2\nP\nH 1 0 0\nH 1 0 0.75\n")

    def fake_run_xtb(_xyz_path, workdir, charge=0, uhf=0, omp_threads=1):
        (workdir / "wbo").write_text("1 2 0.900\n")

    monkeypatch.setattr(pipeline, "_xtb_available", lambda: True)
    monkeypatch.setattr(pipeline, "run_xtb", fake_run_xtb)
    monkeypatch.setattr(pipeline, "OUT_ROOT", tmp_path / "views")
    monkeypatch.setattr(pipeline, "STAGE_ROOT", tmp_path / "stages")
    monkeypatch.setattr(pipeline, "ALIGNMENT_OUT_ROOT", tmp_path / "alignments")
    cfg = pipeline.rp_stage_config()
    cfg["n_seeds"] = 1
    cfg["dwbo_threshold"] = 0.7

    rec = pipeline.process_xyz_stage(
        "h2_direct", r_xyz, p_xyz,
        workdir=tmp_path / "work", stage="rp",
        inner_workers=0, charge=0, multiplicity=1, xtb_mode="auto",
        rp_config=cfg)

    assert rec["slim"]["n_mechs"] >= 1
    assert rec["rp"]["config"]["dwbo_threshold"] == 0.7
    assert (tmp_path / "stages" / "h2_direct" / "rp_stage.json").exists()
    assert (tmp_path / "views" / "h2_direct" / "view.html").exists()


def test_smiles_to_formal_wbo_preserves_explicit_h_and_atom_maps():
    pytest.importorskip("rdkit")

    endpoint = pipeline.smiles_to_formal_wbo("[CH3:1][O:2][H:3]")

    assert endpoint.elements == ["C", "O", "H"]
    assert endpoint.atom_maps == {0: 1, 1: 2, 2: 3}
    assert endpoint.wbo[0, 1] == pytest.approx(1.0)
    assert endpoint.wbo[1, 2] == pytest.approx(1.0)
    assert endpoint.coords.shape == (3, 3)


def test_smiles_to_formal_wbo_can_expand_hydrogens():
    pytest.importorskip("rdkit")

    endpoint = pipeline.smiles_to_formal_wbo(
        "[CH3:1][O:2]", expand_hydrogens=True)

    assert endpoint.elements == ["C", "O", "H", "H", "H"]
    assert endpoint.atom_maps == {0: 1, 1: 2}
    assert endpoint.hydrogen_policy == "expand_hydrogens"
    assert endpoint.coords.shape == (5, 3)
    assert endpoint.wbo[0, 1] == pytest.approx(1.0)
    assert sorted(
        j for j, value in enumerate(endpoint.wbo[0])
        if value == pytest.approx(1.0)
    ) == [1, 2, 3, 4]


def test_process_smiles_stage_runs_rp_from_formal_bond_orders(
        tmp_path, monkeypatch):
    pytest.importorskip("rdkit")

    monkeypatch.setattr(pipeline, "OUT_ROOT", tmp_path / "views")
    monkeypatch.setattr(pipeline, "STAGE_ROOT", tmp_path / "stages")
    monkeypatch.setattr(pipeline, "ALIGNMENT_OUT_ROOT", tmp_path / "alignments")
    cfg = pipeline.rp_stage_config()
    cfg["n_seeds"] = 1

    rec = pipeline.process_smiles_stage(
        "smiles_deprot",
        "[O:1][H:2]",
        "[O-:1].[H+:2]",
        workdir=tmp_path / "work",
        stage="rp",
        inner_workers=0,
        rp_config=cfg)

    assert rec["slim"]["n_mechs"] >= 1
    assert rec["rp"]["config"]["n_seeds"] == 1
    broken = [list(pair) for pair in rec["rp"]["mechanisms"][0]["broken_bonds_R"]]
    assert broken == [[0, 1]]
    assert (tmp_path / "work" / "source" / "R_wbo.json").exists()
    assert (tmp_path / "work" / "source" / "R_graph.json").exists()
    assert (tmp_path / "stages" / "smiles_deprot" / "rp_stage.json").exists()
    assert (tmp_path / "views" / "smiles_deprot" / "view.html").exists()


def test_main_smiles_mode_writes_view_and_eval(tmp_path, monkeypatch):
    pytest.importorskip("rdkit")

    out = tmp_path / "eval.json"
    monkeypatch.setattr(pipeline, "OUT_ROOT", tmp_path / "views")
    monkeypatch.setattr(pipeline, "EVAL_JSON", out)
    monkeypatch.setattr(sys, "argv", [
        "rxn-core",
        "--stage", "rp",
        "--name", "smiles_cli",
        "--reactant-smiles", "[O:1][H:2]",
        "--product-smiles", "[O-:1].[H+:2]",
        "--workdir", str(tmp_path / "work"),
        "--stage-root", str(tmp_path / "stages"),
        "--workers", "1",
    ])

    pipeline.main()

    data = json.loads(out.read_text())
    assert data[0]["step"] == "smiles_cli"
    assert data[0]["n_mechs"] >= 1
    assert (tmp_path / "stages" / "smiles_cli" / "rp_stage.json").exists()
    assert (tmp_path / "views" / "smiles_cli" / "view.html").exists()


def test_process_xyz_stage_can_resume_rp_for_collective_validation(
        tmp_path, monkeypatch):
    r_xyz = tmp_path / "R.xyz"
    p_xyz = tmp_path / "P.xyz"
    r_xyz.write_text("2\nR\nH 0 0 0\nH 0 0 0.75\n")
    p_xyz.write_text("2\nP\nH 1 0 0\nH 1 0 0.75\n")

    def fake_run_xtb(_xyz_path, workdir, charge=0, uhf=0, omp_threads=1):
        (workdir / "wbo").write_text("1 2 0.900\n")

    monkeypatch.setattr(pipeline, "_xtb_available", lambda: True)
    monkeypatch.setattr(pipeline, "run_xtb", fake_run_xtb)
    monkeypatch.setattr(pipeline, "OUT_ROOT", tmp_path / "views")
    monkeypatch.setattr(pipeline, "STAGE_ROOT", tmp_path / "stages")
    monkeypatch.setattr(pipeline, "ALIGNMENT_OUT_ROOT", tmp_path / "alignments")
    cfg = pipeline.rp_stage_config()
    cfg["n_seeds"] = 1

    pipeline.process_xyz_stage(
        "h2_direct", r_xyz, p_xyz,
        workdir=tmp_path / "work", stage="rp",
        inner_workers=0, charge=0, multiplicity=1, xtb_mode="auto",
        rp_config=cfg)

    def fail_run_rp_stage(*_args, **_kwargs):
        raise AssertionError("R-P stage should not rerun")

    monkeypatch.setattr(pipeline, "run_rp_stage", fail_run_rp_stage)
    rec = pipeline.process_xyz_stage(
        "h2_direct", r_xyz, p_xyz,
        workdir=tmp_path / "work", stage="full", resume_rp=True,
        target_specs=[], inner_workers=0, charge=0, multiplicity=1,
        xtb_mode="auto")

    assert rec["ts"]["stage"] == "ts"
    assert rec["slim"]["n_mechs"] >= 1
    assert (tmp_path / "stages" / "h2_direct" / "ts_stage.json").exists()


def test_ts_target_from_xyz_uses_explicit_cache_dirs(tmp_path, monkeypatch):
    ts_xyz = tmp_path / "guess.xyz"
    ts_xyz.write_text("1\nenergy: -1.234567 gnorm: 0.0\nH 0 0 0\n")
    sp_cache = tmp_path / "single_point"
    hess_cache = tmp_path / "hessian"
    calls = []

    def fake_run_xtb(xyz_path, workdir, charge=0, uhf=0, omp_threads=1):
        calls.append(("sp", Path(xyz_path).name, workdir.name, charge, uhf))
        (workdir / "wbo").write_text("")

    def fake_run_xtb_hess(xyz_path, workdir, charge=0, uhf=0, omp_threads=1):
        calls.append(("hess", Path(xyz_path).name, workdir.name, charge, uhf))
        (workdir / "g98.out").write_text("fake")

    monkeypatch.setattr(pipeline, "_xtb_available", lambda: True)
    monkeypatch.setattr(pipeline, "run_xtb", fake_run_xtb)
    monkeypatch.setattr(pipeline, "run_xtb_hess", fake_run_xtb_hess)
    monkeypatch.setattr(
        pipeline, "parse_g98_modes",
        lambda _path: (np.array([-500.0]), np.zeros((1, 1, 3))))

    target = pipeline.ts_target_from_xyz(
        "ig", "guess", ts_xyz,
        sp_workdir=sp_cache, hess_workdir=hess_cache,
        charge=1, multiplicity=2, xtb_mode="auto")

    assert target.kind == "ig"
    assert target.label == "guess"
    assert target.freqs.tolist() == [-500.0]
    assert target.energy == -1.234567
    assert calls == [
        ("sp", "guess.xyz", "single_point", 1, 1),
        ("hess", "guess.xyz", "hessian", 1, 1),
    ]


def test_multiplicity_must_be_positive():
    with pytest.raises(ValueError, match="multiplicity"):
        pipeline._xtb_charge_uhf(0, 0)


def test_sp_cache_can_copy_xyz_fallback_without_xtb(tmp_path):
    sp = tmp_path / "sp_iter1"
    hess = tmp_path / "hess_iter1"
    hess.mkdir()
    fallback = hess / "xtbhess.xyz"
    fallback.write_text("1\nTS\nH 0 0 0\n")
    sp.mkdir()
    (sp / "wbo").write_text("")

    pipeline._ensure_sp_cache(
        sp, "iter1", xyz_fallback=fallback, xtb_mode="cache-only")

    assert (sp / "iter1.xyz").read_text() == fallback.read_text()


def test_xtb_threads_auto_is_capped_per_molecule(monkeypatch):
    monkeypatch.setattr(pipeline.os, "cpu_count", lambda: 64)

    assert pipeline._resolve_xtb_threads("auto") == 8
    assert pipeline._resolve_xtb_threads("32") == 8
    assert pipeline._resolve_xtb_threads("4") == 4
    assert pipeline._resolve_xtb_threads("32", max_threads=12) == 12


def test_default_worker_count_is_not_the_xtb_thread_cap(monkeypatch):
    monkeypatch.setattr(pipeline.os, "cpu_count", lambda: 64)

    assert pipeline._default_worker_count() == 63


def test_viewer_uses_step_level_download_button():
    html = pipeline.HTML.format(
        title="x",
        data_json=(
            '{"step":"s","n_atoms":1,'
            '"reactant":{"elements":["H"],"coords":[[0,0,0]]},'
            '"product":{"elements":["H"],"coords":[[0,0,0]]},'
            '"mechanisms":[],"default_mech_id":null,'
            '"include_gt":false,'
            '"score_config":{"EVENT_WEIGHT_POWER":1,"WBO_PROGRESS_POWER":1}}'
        ),
    )

    assert 'id="downloadAllBtn">Download</button>' in html
    assert 'id="r_meta"' in html
    assert 'id="p_meta"' in html
    assert 'id="showAtomIndices"' in html
    assert 'id="showDegeneracy"' in html
    assert 'id="rOrdered"' in html
    assert 'TS R-order' not in html
    assert '<option value="endpoint">Endpoint</option>' not in html
    assert 'addAtomLabels' in html
    assert 'applyDegeneracyStyles' in html
    assert 'branchBlocks' in html
    assert 'branchColorGroups' in html
    assert 'productAtomLabels' in html
    assert '"R"+i+"/P"+mapping[String(i)]' not in html
    assert "return els.map((_, i) => String(i));" in html
    assert 'const pAligned = !!(rOrdered && mech.product_xyz_in_R)' in html
    assert '"native P"' in html
    assert 'id="zipBtn"' not in html
    assert 'mechanism.json' in html
    assert 'energies_frequencies.json' in html
    assert 'buildEnergyFrequencySummary' in html
    assert 'viewer_data.json' in html
    assert 'view_html:"view.html"' in html
    assert 'root+"/view.html"' in html
    assert 'mechanisms/mechanism_' in html


def test_view_data_carries_endpoint_energy_metadata():
    wbo = np.zeros((1, 1))
    inputs = pipeline.step_inputs_from_arrays(
        "h",
        ["H"], np.array([[0.0, 0.0, 0.0]]), wbo,
        ["H"], np.array([[1.0, 0.0, 0.0]]), wbo,
        energy_R=-1.0, energy_P=-0.9,
    )

    data = pipeline.build_view_data(inputs, {"mechanisms": []})

    assert data["reactant"]["metadata"]["energy_hartree"] == -1.0
    assert data["product"]["metadata"]["energy_hartree"] == -0.9
    assert data["metadata_units"] == {
        "energy": "hartree",
        "frequency": "cm^-1",
    }


def test_view_data_does_not_emit_endpoint_degeneracy_orbits():
    wbo = np.zeros((2, 2))
    wbo[0, 1] = wbo[1, 0] = 0.9
    inputs = pipeline.step_inputs_from_arrays(
        "h2",
        ["H", "H"],
        np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.75]]),
        wbo,
        ["H", "H"],
        np.array([[1.0, 0.0, 0.0], [1.0, 0.0, 0.75]]),
        wbo,
    )

    data = pipeline.build_view_data(inputs, {"mechanisms": []})

    assert "degeneracy" not in data["reactant"]
    assert "degeneracy" not in data["product"]


def test_rp_stage_carries_branch_symmetry_to_mechanism():
    wbo = np.zeros((2, 2))
    wbo[0, 1] = wbo[1, 0] = 0.9
    inputs = pipeline.step_inputs_from_arrays(
        "h2",
        ["H", "H"],
        np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.75]]),
        wbo,
        ["H", "H"],
        np.array([[1.0, 0.0, 0.0], [1.0, 0.0, 0.75]]),
        wbo,
    )
    cfg = pipeline.rp_stage_config()
    cfg["n_seeds"] = 1
    cfg["symmetry_repair"] = False

    result = pipeline.run_rp_stage(inputs, config=cfg, inner_workers=0)

    branch_symmetry = result["mechanisms"][0]["branch_symmetry"]
    assert branch_symmetry["rule"] == "mechanism_dedup_branch_symmetry"
    assert branch_symmetry["dedup_witness_count"] >= 1
    assert any(
        block["r_atoms"] == [0, 1] and block["p_atoms"] == [0, 1]
        for block in branch_symmetry["blocks"]
    )


def test_array_based_mechanism_discovery_returns_aligned_product():
    el = ["H", "H"]
    xyz_r = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.75]])
    xyz_p = np.array([[1.0, 0.0, 0.0], [1.0, 0.0, 0.75]])
    wbo = np.array([[0.0, 0.9], [0.9, 0.0]])
    cfg = pipeline.rp_stage_config()
    cfg["n_seeds"] = 1

    result = pipeline.discover_mechanisms_from_arrays(
        el, xyz_r, wbo, el, xyz_p, wbo,
        step_name="h2", config=cfg, inner_workers=0)

    assert result["stage"] == "rp"
    assert result["mechanisms"]
    mech = result["mechanisms"][0]
    assert set(mech["mapping_RP"]) == {0, 1}
    assert len(mech["product_xyz_in_R"]) == 2


def test_rp_alignment_file_export_is_neb_ready(tmp_path):
    el = ["H", "H"]
    xyz_r = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.75]])
    xyz_p = np.array([[1.0, 0.0, 0.0], [1.0, 0.0, 0.75]])
    wbo = np.array([[0.0, 0.9], [0.9, 0.0]])
    inputs = pipeline.step_inputs_from_arrays(
        "h2", el, xyz_r, wbo, el, xyz_p, wbo)
    cfg = pipeline.rp_stage_config()
    cfg["n_seeds"] = 1
    rp = pipeline.run_rp_stage(inputs, config=cfg)

    out = pipeline.write_rp_alignment_files(inputs, rp, tmp_path)

    assert out["n_mechanisms"] == len(rp["mechanisms"])
    assert (tmp_path / "manifest.json").exists()
    assert (tmp_path / "R.xyz").exists()
    assert (tmp_path / "P_original.xyz").exists()
    mdir = tmp_path / "mechanisms" / "mechanism_001"
    assert (mdir / "R.xyz").exists()
    assert (mdir / "P_aligned.xyz").exists()
    assert (mdir / "neb_endpoints.xyz").read_text().count("\n2\n") == 1
    assert "R_index,R_element,P_index,P_element" in (
        mdir / "mapping_R_to_P.csv").read_text()


def test_ts_alignment_file_export_writes_best_score_core_frame(tmp_path):
    el = ["H", "H"]
    xyz_r = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.75]])
    xyz_p = np.array([[1.0, 0.0, 0.0], [1.0, 0.0, 0.75]])
    wbo = np.array([[0.0, 0.9], [0.9, 0.0]])
    inputs = pipeline.step_inputs_from_arrays(
        "h2", el, xyz_r, wbo, el, xyz_p, wbo)
    ts_result = {
        "stage": "ts",
        "step": "h2",
        "mechanisms": [{
            "id": 1,
            "label": "#1",
            "broken_bonds_R": [[0, 1]],
            "formed_bonds_R": [],
            "core_atoms": [0, 1],
            "gt": {
                "S": 1.2,
                "beta": 0.8,
                "wbo_progress": 0.9,
                "wbo_progress_factor": 0.9,
                "freq": -500.0,
                "k": 0,
                "energy_hartree": -2.5,
                "frequency_units": "cm^-1",
                "frequencies_cm1": [-500.0, 100.0],
                "frequency_summary": {
                    "n_modes": 2,
                    "n_imaginary": 1,
                    "imaginary_cm1": [-500.0],
                    "lowest_cm1": -500.0,
                    "highest_cm1": 100.0,
                },
                "core_map": {"0": 1, "1": 0},
                "core_sources": ["R"],
                "core_pool_dedup_count": 1,
                "elements": el,
                "xyz": xyz_p.tolist(),
                "xyz_in_R": [[1.0, 0.0, 0.75], [1.0, 0.0, 0.0]],
                "picked_disp_R": [[0.1, 0.0, 0.0], [-0.1, 0.0, 0.0]],
                "broken_bonds_T": [[1, 0]],
                "formed_bonds_T": [],
            },
            "igs": [],
        }],
    }

    out = pipeline.write_ts_alignment_files(inputs, ts_result, tmp_path)

    assert out["n_targets"] == 1
    tdir = tmp_path / "mechanisms" / "mechanism_001" / "gt_GT"
    assert (tdir / "TS_native.xyz").exists()
    assert (tdir / "TS_core_aligned_R_frame.xyz").exists()
    assert (tdir / "picked_mode_R_frame.xyz").exists()
    score_doc = pipeline.read_stage_json(tdir / "score.json")
    assert score_doc["energy_hartree"] == -2.5
    assert score_doc["frequencies_cm1"] == [-500.0, 100.0]
    assert score_doc["frequency_summary"]["n_imaginary"] == 1
    assert "spectator bijections" in (tdir / "score.json").read_text()
