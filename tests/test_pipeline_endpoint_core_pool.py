import pytest

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

    def fake_run_xtb(xyz_path, workdir):
        calls.append((xyz_path, workdir))
        (workdir / "wbo").write_text("1 2 0.900\n")

    monkeypatch.setattr(pipeline, "_xtb_available", lambda: True)
    monkeypatch.setattr(pipeline, "run_xtb", fake_run_xtb)

    pipeline._ensure_sp_cache(cache, "R", xtb_mode="auto")

    assert calls == [(cache / "r.xyz", cache)]
    assert (cache / "wbo").exists()


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
