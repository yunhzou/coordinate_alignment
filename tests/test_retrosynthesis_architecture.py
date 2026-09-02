from rxn_core.alignment.post_aam import (
    AAMHierarchy,
    AtomPermutation,
    FragmentMatch,
    SymmetryDomain,
)
from rxn_core.fragment_matching import FragmentCandidate
from rxn_core.retrosynthesis.enumeration import (
    CoverageEnumerationConfig,
    enumerate_coverage_patterns,
)
from rxn_core.retrosynthesis.compressed_coverage import (
    CoverageRecommendationConfig,
    assign_domain_signatures,
    candidate_target_domains,
    coverage_signature,
    recommend_compressed_coverage_patterns,
)
from rxn_core.fragment_matching.serialization import (
    fragment_candidate_from_record,
    fragment_candidate_to_record,
)


def test_fragment_candidate_serialization_preserves_fragment_units():
    hierarchy = AAMHierarchy((FragmentMatch(
        fragment_index=2,
        island_index=4,
        r_atoms=(0, 1, 2),
        deferred_edges=((1, 2),),
        symmetry_domains=(SymmetryDomain(
            r_atoms=(1, 2),
            p_atoms=(5, 6),
            source="sym_block",
            extendable=True,
        ),),
        representative_assignments=((0, 3), (1, 5), (2, 6)),
        exact_fixed=(0,),
        multiplicity=2,
        target_generators=(AtomPermutation((0, 1, 2, 3, 4, 6, 5)),),
    ),))
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
        aam_hierarchy=hierarchy,
    )

    record = fragment_candidate_to_record(candidate)
    restored = fragment_candidate_from_record(record)

    assert record["retained_fragments"] == [[0, 1], [2]]
    assert record["mapping"] == [[0, 0], [1, 1], [2, 2]]
    assert restored.aam_hierarchy == hierarchy


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


def test_repeated_symmetric_ligand_assembles_without_bijection_expansion():
    ligand = FragmentCandidate(
        source_id="ligand",
        mapping=((0, 1),),
        retained_atoms=(0,),
        covered_target_atoms=(1,),
        leftover_fragments=(),
        boundary_bonds=(),
        attachment_atoms_source=(),
        attachment_atoms_target=(),
        copied_residual_placements=(),
        augmented_target_atom_count=4,
        retained_fragments=((0,),),
        aam_hierarchy=AAMHierarchy((FragmentMatch(
            fragment_index=0,
            island_index=0,
            r_atoms=(0,),
            symmetry_domains=(SymmetryDomain(
                r_atoms=(0,),
                p_atoms=(1, 2, 3),
                source="sym_block",
            ),),
            representative_assignments=((0, 1),),
        ),)),
    )
    ligand_signature = coverage_signature(
        candidate_target_domains(ligand))
    center_signature = ((0,),)

    result = recommend_compressed_coverage_patterns(
        (center_signature, ligand_signature),
        4,
        lambda pattern, covered: (-covered, tuple(map(str, pattern))),
        result_limit=4,
        config=CoverageRecommendationConfig(
            maximum_precursors=4,
            frontier_limit=20,
        ),
    )

    assert result.patterns == ((
        center_signature,
        ligand_signature,
        ligand_signature,
        ligand_signature,
    ),)
    witness = assign_domain_signatures(result.patterns[0], 4)
    assert witness[0] == (0,)
    assert {targets[0] for targets in witness[1:]} == {1, 2, 3}
