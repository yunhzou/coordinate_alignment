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
