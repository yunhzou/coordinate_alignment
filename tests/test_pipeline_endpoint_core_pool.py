import pytest
import numpy as np
from pathlib import Path

import rxn_core.pipeline as pipeline
from rxn_core.pipeline import (
    _merge_endpoint_core_pools,
    _product_core_pool_to_reactant,
)


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
    ts_xyz.write_text("1\nTS\nH 0 0 0\n")
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
            '"include_gt":false,"score_config":{"W_RXN":1,"W_CORE":0.2,"IMAG_PEN":0.3}}'
        ),
    )

    assert 'id="downloadAllBtn">Download</button>' in html
    assert 'id="showAtomIndices"' in html
    assert 'addAtomLabels' in html
    assert 'id="zipBtn"' not in html
    assert 'mechanism.json' in html
    assert 'viewer_data.json' in html
    assert 'mechanisms/mechanism_' in html


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
                "rho": 0.3,
                "kappa": 0.2,
                "freq": -500.0,
                "k": 0,
                "n_imag": 1,
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
    assert "spectator bijections" in (tdir / "score.json").read_text()
