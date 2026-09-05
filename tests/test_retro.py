from collections import Counter

import numpy as np
import pytest

from rxn_core import (
    WeightedGraph,
)
from rxn_core.fragment_matching import (
    FragmentCandidate,
    FragmentDetectionConfig,
    FragmentDetectionExecution,
    detect_fragments,
    detect_fragments_parallel,
    materialize_target_coverage_orbit,
    prepare_fragment_target,
    progressive_fragment_matching,
)
from rxn_core.fragment_matching.augmentation import match_augmented_residuals
from rxn_core.retrosynthesis import assemble_fragment_cover
from rxn_core.smiles import smiles_to_weighted_graph


def test_prepared_fragment_target_matches_direct_detection():
    config = FragmentDetectionConfig(candidate_limit=100)
    source = smiles_to_weighted_graph(
        "CCBr", expand_hydrogens=True)
    target = smiles_to_weighted_graph(
        "CCO", expand_hydrogens=True)

    direct = detect_fragments(source, target, config=config)
    prepared = prepare_fragment_target(target, config=config)
    reused = detect_fragments(source, prepared, config=config)

    assert reused == direct


def test_progressive_matching_recovers_all_t05_ground_truth_atoms():
    target = smiles_to_weighted_graph(
        "C/C(NC1=C(C(C)C)C=CC=C1C(C)C)=C/C(C)=N/"
        "C2=NC3=C(C=C(C#N)C=C3)S2",
        expand_hydrogens=True,
    )
    sources = tuple(
        (source_id, smiles_to_weighted_graph(smiles, expand_hydrogens=True))
        for source_id, smiles in (
            ("acetylacetone", "CC(=O)CC(C)=O"),
            ("aniline", "CC(C)c1cccc(C(C)C)c1N"),
            ("benzothiazole", "N#Cc1ccc2sc(N)nc2c1"),
        )
    )

    result = progressive_fragment_matching(
        sources,
        target,
        config=FragmentDetectionConfig(
            seed_mode="all",
            branch_limit=100,
            candidate_limit=512,
            iso_tolerance=0.5,
        ),
    )

    assert result.uncovered_target_atoms == ()
    by_id = {placement.source_id: placement for placement in result.placements}
    acetylacetone = dict(by_id["acetylacetone"].mapping)
    assert {acetylacetone[index] for index in (0, 1, 3, 4, 5)} == {
        0, 1, 15, 16, 17,
    }
    assert 2 not in acetylacetone
    assert 6 not in acetylacetone
    benzothiazole = dict(by_id["benzothiazole"].mapping)
    assert 0 in benzothiazole
    assert 1 in benzothiazole
    occupied = [
        target_atom
        for placement in result.placements
        for _source_atom, target_atom in placement.mapping
    ]
    assert len(occupied) == len(set(occupied)) == len(target.nodes)
    assert result.selections
    for selection in result.selections:
        assert selection.candidate.derivations
        original_pairs = {(selection.source_atoms[a], selection.target_atoms[b])
                          for a, b in selection.candidate.mapping}
        assert original_pairs.issubset(set(result.placements[selection.source_index].mapping))


def test_fragment_detection_seed_limit_is_reported_as_incomplete():
    config = FragmentDetectionConfig(seed_limit=1)
    result = detect_fragments(
        smiles_to_weighted_graph("CCBr", expand_hydrogens=True),
        smiles_to_weighted_graph("CCO", expand_hydrogens=True),
        config=config,
    )

    assert not result.complete
    assert result.status == "seed_limited"


def _matrix(size, edges):
    matrix = np.zeros((size, size), dtype=float)
    for left, right, weight in edges:
        matrix[left, right] = matrix[right, left] = weight
    return matrix


def test_noncompetitive_singleton_is_recorded_without_augmented_copy():
    precursor = WeightedGraph(
        ["C", "Br"],
        _matrix(2, [(0, 1, 1.0)]),
    )
    target = WeightedGraph(
        ["C", "C"],
        _matrix(2, [(0, 1, 1.0)]),
    )

    result = detect_fragments(
        precursor,
        target,
        source_id="carbon-bromide",
        config=FragmentDetectionConfig(
            minimum_fragment_size=1,
            branch_limit=8,
            maximum_boundary_bonds=1,
            maximum_leftover_fragments=1,
        ),
    )

    assert result.status == "matched"
    assert result.complete
    assert result.best_fragment_size == 1
    assert result.candidates
    candidate = result.candidates[0]
    assert candidate.leftover_fragments == ((1,),)
    assert candidate.boundary_bonds == ((0, 1),)
    assert candidate.copied_residual_placements == ()
    assert candidate.augmented_target_atom_count == 2


def test_augmented_copy_baseline_survives_invalid_greedy_target_mapping():
    source = WeightedGraph(
        ["C", "O"],
        _matrix(2, [(0, 1, 1.0)]),
    ).to_networkx()
    target = WeightedGraph(
        ["C", "O"],
        _matrix(2, []),
    ).to_networkx()

    result = match_augmented_residuals(
        source,
        target,
        {0: 0},
        {1},
        ((0, 1),),
        graph_floor=0.2,
        iso_tolerance=0.5,
        branch_limit=100,
    )

    assert not result.capped
    assert result.augmented_target_atom_count == 2
    assert tuple(item.mapping for item in result.placements) == (((0, 0),),)
    assert result.placements[0].hierarchy.fragments


def test_augmented_leftover_competes_for_unused_target_atoms():
    precursor = WeightedGraph(
        ["O", "C", "O"],
        _matrix(3, [(0, 1, 2.0), (1, 2, 2.0)]),
    )
    target = WeightedGraph(
        ["O", "C", "O"],
        _matrix(3, [(0, 1, 2.0), (1, 2, 1.0)]),
    )

    result = detect_fragments(
        precursor,
        target,
        source_id="carbon-dioxide",
        config=FragmentDetectionConfig(
            minimum_fragment_size=1,
            iso_tolerance=0.5,
            branch_limit=100,
            candidate_limit=100,
        ),
    )

    assert result.status == "matched"
    assert result.complete
    assert result.best_fragment_size == 3
    assert all(candidate.covered_target_atoms == (0, 1, 2)
               for candidate in result.candidates)
    full = [candidate for candidate in result.candidates
            if candidate.covered_target_atoms == (0, 1, 2)]
    assert full
    assert full[0].retained_atoms == (0, 1, 2)
    assert full[0].retained_fragments in (
        ((0,), (1, 2)),
        ((0, 1), (2,)),
    )
    assert full[0].leftover_fragments == ()
    assert len(full[0].boundary_bonds) == 1
    assert full[0].copied_residual_placements == ()


def test_target_region_recovers_small_donor_despite_larger_incidental_match():
    target = smiles_to_weighted_graph(
        "CNC(=O)CO", expand_hydrogens=True)
    donor = smiles_to_weighted_graph(
        "CSCC(=O)CO", expand_hydrogens=True)
    target_graph = target.to_networkx()
    methyl_carbon = 0
    methyl_region = {methyl_carbon} | {
        neighbor for neighbor in target_graph.neighbors(methyl_carbon)
        if target_graph.nodes[neighbor]["element"] == "H"
    }
    config = FragmentDetectionConfig(
        minimum_fragment_size=1,
        branch_limit=100,
        candidate_limit=512,
    )

    untargeted = detect_fragments(donor, target, config=config)
    targeted = detect_fragments(
        donor,
        target,
        config=config,
        target_region_atoms=methyl_region,
    )

    assert any(
        methyl_region.issubset(candidate.covered_target_atoms)
        for candidate in untargeted.candidates
    )
    assert any(
        not methyl_region.intersection(candidate.covered_target_atoms)
        for candidate in untargeted.candidates
    )
    assert any(
        set(candidate.covered_target_atoms) == methyl_region
        for candidate in targeted.candidates
    )
    assert all(
        methyl_region.intersection(candidate.covered_target_atoms)
        for candidate in targeted.candidates
    )


def test_branch_cap_is_reported_as_incomplete():
    precursor = WeightedGraph(
        ["C", "C"],
        _matrix(2, [(0, 1, 1.0)]),
    )
    target = WeightedGraph(
        ["C"] * 6,
        _matrix(6, [(index, index + 1, 1.0) for index in range(5)]),
    )

    result = detect_fragments(
        precursor,
        target,
        source_id="ambiguous-CC",
        config=FragmentDetectionConfig(branch_limit=2),
    )

    assert result.status == "capped"
    assert not result.complete
    assert result.capped_seed_count > 0
    assert result.maximum_branch_count > result.branch_limit


def test_detection_keeps_target_automorphism_family_compressed_until_assembly():
    precursor = WeightedGraph(["N"], _matrix(1, []))
    target = WeightedGraph(
        ["C", "N", "N", "N"],
        _matrix(4, [(0, 1, 1.0), (0, 2, 1.0), (0, 3, 1.0)]),
    )

    result = detect_fragments(
        precursor,
        target,
        source_id="symmetric-ligand",
        config=FragmentDetectionConfig(candidate_limit=100),
    )

    assert len(result.candidates) == 1
    variants = materialize_target_coverage_orbit(
        result.candidates[0], target)
    assert {candidate.covered_target_atoms for candidate in variants} == {
        (1,), (2,), (3,),
    }
    assert result.candidates[0].aam_hierarchy.fragments


def test_rough_seed_mode_stops_after_majority_source_fragment():
    precursor = WeightedGraph(
        ["C", "C", "C", "C"],
        _matrix(4, [(0, 1, 1.0), (1, 2, 1.0), (2, 3, 1.0)]),
    )
    target = WeightedGraph(
        ["C", "C", "C", "C"],
        _matrix(4, [(0, 1, 1.0), (1, 2, 1.0), (2, 3, 1.0)]),
    )

    result = detect_fragments(
        precursor,
        target,
        config=FragmentDetectionConfig(
            seed_mode="fragment_cover",
            rough_retention_threshold=0.5,
        ),
    )

    assert result.best_fragment_size == 4
    assert result.seed_attempt_count == 1
    assert result.seed_pruned_count == 3
    assert result.rough_stop_hit
    assert result.status == "rough"
    assert not result.complete


def test_orbit_representative_seed_mode_is_explicitly_rough():
    precursor = WeightedGraph(
        ["C", "C", "C", "C"],
        _matrix(4, [(0, 1, 1.0), (1, 2, 1.0), (2, 3, 1.0)]),
    )
    target = WeightedGraph(
        ["C", "C", "C", "C"],
        _matrix(4, [(0, 1, 1.0), (1, 2, 1.0), (2, 3, 1.0)]),
    )

    result = detect_fragments(
        precursor,
        target,
        config=FragmentDetectionConfig(seed_mode="orbit_representatives"),
    )

    assert result.best_fragment_size == 4
    assert result.seed_attempt_count == 2
    assert result.seed_pruned_count == 2
    assert result.status == "rough"
    assert not result.complete


def test_parallel_seed_execution_preserves_exact_detection_result():
    precursor = WeightedGraph(
        ["C", "C", "O", "H", "H"],
        _matrix(5, [
            (0, 1, 1.0), (1, 2, 1.0), (0, 3, 1.0), (1, 4, 1.0),
        ]),
    )
    target = WeightedGraph(
        ["C", "C", "O", "H", "H", "N"],
        _matrix(6, [
            (0, 1, 1.0), (1, 2, 1.0), (0, 3, 1.0), (1, 4, 1.0),
            (2, 5, 1.0),
        ]),
    )
    config = FragmentDetectionConfig(seed_mode="all")

    sequential = detect_fragments(precursor, target, config=config)
    parallel = detect_fragments_parallel(
        precursor,
        target,
        config=config,
        execution=FragmentDetectionExecution(seed_workers=2),
    )

    assert parallel == sequential


def test_balanced_williamson_reaction_is_recovered_with_hidden_side_product():
    pytest.importorskip("rdkit")
    from rdkit import Chem
    from rxn_core.smiles import smiles_to_weighted_graph

    reactant_smiles = ("CCBr", "C[O-].[Na+]")
    product_smiles = ("CCOC", "[Na+].[Br-]")

    def composition(smiles_group):
        return Counter(
            atom.GetSymbol()
            for smiles in smiles_group
            for atom in Chem.MolFromSmiles(smiles).GetAtoms()
        )

    assert composition(reactant_smiles) == composition(product_smiles)

    target = smiles_to_weighted_graph("CCOC", expand_hydrogens=False)
    config = FragmentDetectionConfig(
        minimum_fragment_size=2,
        branch_limit=100,
        candidate_limit=100,
        maximum_boundary_bonds=1,
        maximum_leftover_fragments=1,
    )
    bromoethane = detect_fragments(
        smiles_to_weighted_graph("CCBr", expand_hydrogens=False),
        target,
        source_id="bromoethane",
        config=config,
    )
    sodium_methoxide = detect_fragments(
        smiles_to_weighted_graph("C[O-].[Na+]", expand_hydrogens=False),
        target,
        source_id="sodium-methoxide",
        config=config,
    )

    assert bromoethane.complete and bromoethane.best_fragment_size == 2
    assert sodium_methoxide.complete and sodium_methoxide.best_fragment_size == 2
    assert all(candidate.leftover_fragments == ((2,),)
               for candidate in bromoethane.candidates)
    assert all(candidate.leftover_fragments == ((2,),)
               for candidate in sodium_methoxide.candidates)

    assembly_result = assemble_fragment_cover(
        target,
        bromoethane.candidates + sodium_methoxide.candidates,
        maximum_precursors=2,
        require_attachment_bonds=True,
    )

    assert assembly_result.status == "matched"
    assert assembly_result.complete
    assert assembly_result.assemblies
    assembly = assembly_result.assemblies[0]
    assert set(assembly.precursor_ids) == {"bromoethane", "sodium-methoxide"}
    assert len(assembly.formed_bonds) == 1
    assert len(assembly.broken_bonds) == 1


def test_symmetric_cover_can_reuse_one_precursor_three_times():
    target = WeightedGraph(
        ["C", "N", "N", "N"],
        _matrix(4, [(0, 1, 1.0), (0, 2, 1.0), (0, 3, 1.0)]),
    )

    def candidate(precursor_id, source, image):
        return FragmentCandidate(
            source_id=precursor_id,
            mapping=((source, image),),
            retained_atoms=(source,),
            covered_target_atoms=(image,),
            leftover_fragments=(),
            boundary_bonds=(),
            attachment_atoms_source=(source,),
            attachment_atoms_target=(image,),
            copied_residual_placements=(),
            augmented_target_atom_count=4,
        )

    candidates = (
        candidate("carbon-core", 0, 0),
        candidate("same-ligand", 0, 1),
        candidate("same-ligand", 0, 2),
        candidate("same-ligand", 0, 3),
    )
    repeated = assemble_fragment_cover(
        target, candidates, maximum_precursors=4,
        allow_repeated_precursors=True)
    distinct_only = assemble_fragment_cover(
        target, candidates, maximum_precursors=4,
        allow_repeated_precursors=False)

    assert repeated.status == "matched"
    assert repeated.assemblies[0].precursor_ids.count("same-ligand") == 3
    assert len(repeated.assemblies[0].formed_bonds) == 3
    assert distinct_only.status == "no_cover"


def test_two_step_triphenylamine_route_materializes_symmetric_placements():
    pytest.importorskip("rdkit")
    from rxn_core.smiles import smiles_to_weighted_graph

    def graph(smiles):
        return smiles_to_weighted_graph(smiles, expand_hydrogens=True)
    config = FragmentDetectionConfig(
        minimum_fragment_size=1,
        iso_tolerance=0.5,
        branch_limit=100,
        candidate_limit=512,
    )
    bromobenzene_smiles = "Brc1ccccc1"
    aniline_smiles = "Nc1ccccc1"
    diphenylamine_smiles = "c1ccc(Nc2ccccc2)cc1"
    triphenylamine_smiles = "c1ccc(N(c2ccccc2)c2ccccc2)cc1"

    def cover(target_smiles, sources, maximum_precursors):
        target = graph(target_smiles)
        candidates = tuple(
            candidate
            for source_id, source_smiles in sources
            for candidate in detect_fragments(
                graph(source_smiles),
                target,
                source_id=source_id,
                config=config,
            ).candidates
        )
        return assemble_fragment_cover(
            target,
            candidates,
            maximum_precursors=maximum_precursors,
            allow_repeated_precursors=True,
            require_attachment_bonds=False,
        )

    step_one = cover(diphenylamine_smiles, (
        ("aniline", aniline_smiles),
        ("bromobenzene", bromobenzene_smiles),
    ), 2)
    step_two = cover(triphenylamine_smiles, (
        ("diphenylamine", diphenylamine_smiles),
        ("bromobenzene", bromobenzene_smiles),
    ), 2)
    direct = cover(triphenylamine_smiles, (
        ("aniline", aniline_smiles),
        ("bromobenzene", bromobenzene_smiles),
    ), 3)

    assert step_one.status == "matched"
    assert Counter(step_one.assemblies[0].precursor_ids) == {
        "aniline": 1,
        "bromobenzene": 1,
    }
    assert step_two.status == "matched"
    assert Counter(step_two.assemblies[0].precursor_ids) == {
        "diphenylamine": 1,
        "bromobenzene": 1,
    }
    assert direct.status == "matched"
    assert Counter(direct.assemblies[0].precursor_ids) == {
        "aniline": 1,
        "bromobenzene": 2,
    }


def test_large_symmetric_star_uses_three_repeated_biphenyl_arms():
    pytest.importorskip("rdkit")
    from rxn_core.smiles import smiles_to_weighted_graph

    def graph(smiles):
        return smiles_to_weighted_graph(smiles, expand_hydrogens=True)

    target_smiles = (
        "c1ccc(-c2ccc(-c3cc(-c4ccc(-c5ccccc5)cc4)cc"
        "(-c4ccc(-c5ccccc5)cc4)c3)cc2)cc1"
    )
    sources = (
        ("1,3,5-tribromobenzene", "Brc1cc(Br)cc(Br)c1"),
        (
            "4-biphenylboronic acid MIDA ester",
            "CN1CC(=O)OB(c2ccc(-c3ccccc3)cc2)OC(=O)C1",
        ),
    )
    target = graph(target_smiles)
    config = FragmentDetectionConfig(
        minimum_fragment_size=1,
        iso_tolerance=0.5,
        branch_limit=100,
        candidate_limit=512,
    )
    detections = [
        detect_fragments(
            graph(source_smiles),
            target,
            source_id=source_id,
            config=config,
        )
        for source_id, source_smiles in sources
    ]
    result = assemble_fragment_cover(
        target,
        tuple(
            candidate
            for detection in detections
            for candidate in detection.candidates
        ),
        maximum_precursors=4,
        allow_repeated_precursors=True,
        require_attachment_bonds=False,
    )

    assert all(detection.complete for detection in detections)
    assert result.status == "matched"
    assembly = result.assemblies[0]
    assert Counter(assembly.precursor_ids) == {
        "1,3,5-tribromobenzene": 1,
        "4-biphenylboronic acid MIDA ester": 3,
    }
    assert len(assembly.formed_bonds) == 3
    assert len(assembly.broken_bonds) == 6


def test_chloroform_plus_three_repeated_pyrazoles_covers_target():
    pytest.importorskip("rdkit")
    from rxn_core.smiles import smiles_to_weighted_graph

    target = smiles_to_weighted_graph(
        "N1(C(N2C=CC=N2)N3N=CC=C3)N=CC=C1",
        expand_hydrogens=True)
    config = FragmentDetectionConfig(
        minimum_fragment_size=1, branch_limit=100, candidate_limit=1000)
    core = detect_fragments(
        smiles_to_weighted_graph("ClC(Cl)Cl", expand_hydrogens=True),
        target, source_id="CHCl3", config=config)
    ligand = detect_fragments(
        smiles_to_weighted_graph("c1cn[nH]c1", expand_hydrogens=True),
        target, source_id="pyrazole", config=config)

    result = assemble_fragment_cover(
        target, core.candidates + ligand.candidates,
        maximum_precursors=4, allow_repeated_precursors=True,
        require_attachment_bonds=False)
    strict_first_pass = assemble_fragment_cover(
        target, core.candidates + ligand.candidates,
        maximum_precursors=4, allow_repeated_precursors=True,
        require_attachment_bonds=True)

    assert result.status == "matched"
    assembly = result.assemblies[0]
    assert Counter(assembly.precursor_ids) == {
        "CHCl3": 1, "pyrazole": 3}
    assert len(assembly.formed_bonds) == 3
    assert len(assembly.broken_bonds) == 6
    assert strict_first_pass.status == "matched"
    assert len(strict_first_pass.assemblies[0].formed_bonds) == 3


def test_known_mcule_suzuki_precursors_cover_ortho_chlorobiphenyl():
    pytest.importorskip("rdkit")
    from rxn_core.smiles import smiles_to_weighted_graph

    # Both exact structures occur in the downloaded Mcule building-block file:
    #   BrC1=CC=CC=C1                 MCULE-5539191636
    #   C1(=CC=CC=C1Cl)B(O)O         MCULE-6011753091
    target = smiles_to_weighted_graph(
        "Clc1ccccc1-c1ccccc1", expand_hydrogens=False)
    bromobenzene = smiles_to_weighted_graph(
        "BrC1=CC=CC=C1", expand_hydrogens=False)
    chlorophenylboronic_acid = smiles_to_weighted_graph(
        "C1(=CC=CC=C1Cl)B(O)O", expand_hydrogens=False)
    config = FragmentDetectionConfig(
        branch_limit=64,
        maximum_boundary_bonds=1,
        maximum_leftover_fragments=1,
    )

    bromine_result = detect_fragments(
        bromobenzene,
        target,
        source_id="MCULE-5539191636",
        config=config,
    )
    boron_result = detect_fragments(
        chlorophenylboronic_acid,
        target,
        source_id="MCULE-6011753091",
        config=config,
    )

    assert bromine_result.complete and bromine_result.best_fragment_size == 6
    assert boron_result.complete and boron_result.best_fragment_size == 7
    assert any(candidate.leftover_fragments == ((0,),)
               for candidate in bromine_result.candidates)
    assert any(len(candidate.leftover_fragments[0]) == 3
               for candidate in boron_result.candidates)

    assembly_result = assemble_fragment_cover(
        target,
        bromine_result.candidates + boron_result.candidates,
        maximum_precursors=2,
        require_attachment_bonds=True,
    )

    assert assembly_result.status == "matched"
    assert assembly_result.complete
    assert assembly_result.assemblies
    assembly = assembly_result.assemblies[0]
    assert set(assembly.precursor_ids) == {
        "MCULE-5539191636",
        "MCULE-6011753091",
    }
    assert len(assembly.formed_bonds) == 1
    assert len(assembly.broken_bonds) == 2
