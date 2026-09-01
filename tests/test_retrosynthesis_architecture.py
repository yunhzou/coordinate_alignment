from rxn_core.fragment_matching import FragmentCandidate
from rxn_core.retrosynthesis.enumeration import (
    CoverageEnumerationConfig,
    enumerate_coverage_patterns,
)
from rxn_core.fragment_matching.serialization import fragment_candidate_to_record
from rxn_core.retrosynthesis.ownership import resolve_overlapping_ownership


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


def test_overlap_masks_are_deferred_to_joint_ownership():
    masks = (0b011, 0b110)

    result = enumerate_coverage_patterns(
        masks,
        3,
        lambda pattern, _covered: tuple(pattern),
        config=CoverageEnumerationConfig(
            maximum_precursors=2,
            mode="exhaustive",
            allow_overlaps=True,
        ),
    )

    assert result.patterns == ((0b011, 0b110),)


def test_joint_ownership_rebuilds_overlapping_precursor_entries():
    def entry(identifier, smiles, mapping):
        atom_count = 8 if smiles == "CC" else 6
        return {
            "precursor_id": identifier,
            "smiles": smiles,
            "row_index": 0,
            "complete": True,
            "status": "matched",
            "best_fragment_size": 2,
            "covered_target_atoms": sorted(target for _, target in mapping),
            "retained_atoms": sorted(source for source, _ in mapping),
            "leftover_fragments": [],
            "boundary_bonds": [],
            "attachment_atoms_target": [],
            "mapping": [list(pair) for pair in mapping],
            "retained_fragments": [sorted(source for source, _ in mapping)],
            "leftover_atom_count": atom_count - 2,
            "structure_key": smiles,
            "retained_heavy_atoms": 2,
            "total_heavy_atoms": 2,
            "heavy_atom_retention": 1.0,
            "retained_atom_count": 2,
            "total_atom_count": atom_count,
            "atom_retention": 2 / atom_count,
            "attachment_trimmed_target_atoms": [],
            "chirality_violations": 0,
            "chirality_violation_target_atoms": [],
        }

    result = resolve_overlapping_ownership(
        (
            entry("left", "CC", ((0, 0), (1, 1))),
            entry("right", "CO", ((0, 1), (1, 2))),
        ),
        3,
        ((0, 1), (1, 2)),
        beam_width=20,
        assembly_limit=10,
    )

    assert result.assemblies
    for assembly in result.assemblies:
        regions = [
            set(precursor["covered_target_atoms"])
            for precursor in assembly["precursors"]
        ]
        assert set.union(*regions) == {0, 1, 2}
        assert sum(map(len, regions)) == 3
