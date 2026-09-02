from rxn_core._native import paired_mapping_invariant


def test_native_paired_mapping_invariant_refines_sparse_relations():
    result = paired_mapping_invariant(
        ((1, 1), (1, 1), (1, 1)),
        (0, 0),
        ((0, 1, 5, 5),),
    )

    assert result == (
        (((1, 1), 3),),
        ((0, 2), (1, 1)),
        (
            ((0, 0, (5, 5)), 1),
            ((0, 1, (0, 0)), 2),
        ),
    )
