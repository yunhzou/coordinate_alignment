from rxn_core.fragment_matching import FragmentCandidate
from rxn_core.retrosynthesis.enumeration import (
    CoverageEnumerationConfig,
    enumerate_coverage_patterns,
)
from rxn_core.fragment_matching.serialization import fragment_candidate_to_record


def test_fragment_candidate_serialization_preserves_fragment_units():
    candidate = FragmentCandidate(
        source_id="CO2",
        mapping=((0, 0), (1, 1), (2, 2)),
        retained_atoms=(0, 1, 2),
        covered_target_atoms=(0, 1, 2),
        leftover_fragments=(),
        boundary_bonds=((1, 2),),
        attachment_atoms_source=(1, 2),
        attachment_atoms_target=(1, 2),
        copied_residual_placements=(),
        augmented_target_atom_count=4,
        retained_fragments=((0, 1), (2,)),
    )

    record = fragment_candidate_to_record(candidate)

    assert record["retained_fragments"] == [[0, 1], [2]]
    assert record["mapping"] == [[0, 0], [1, 1], [2, 2]]


def test_all_coverage_modes_find_the_same_simple_exact_covers():
    masks = (0b0011, 0b1100, 0b0101, 0b1010)

    def rank_pattern(pattern, _covered_atom_count):
        return tuple(pattern)

    expected = {(0b0011, 0b1100), (0b0101, 0b1010)}
    for mode in ("exhaustive", "modular", "recommendation"):
        result = enumerate_coverage_patterns(
            masks,
            4,
            rank_pattern,
            config=CoverageEnumerationConfig(
                maximum_precursors=2,
                mode=mode,
                beam_width=20,
                patterns_per_coverage=4,
                state_limit=100,
            ),
        )
        assert set(result.patterns) == expected
        assert result.complete


def test_overlapping_fragments_cannot_manufacture_an_exact_cover():
    result = enumerate_coverage_patterns(
        (0b011, 0b110),
        3,
        lambda pattern, _covered: tuple(pattern),
        config=CoverageEnumerationConfig(
            maximum_precursors=2,
            mode="exhaustive",
        ),
    )

    assert result.patterns == ()
