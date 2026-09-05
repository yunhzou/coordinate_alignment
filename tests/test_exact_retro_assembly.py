"""Small exhaustive oracles for the recommendation layer, not known reactions."""
from itertools import combinations
from types import SimpleNamespace

import numpy as np
from rdkit import Chem

from rxn_core import WeightedGraph
from rxn_core.retrosynthesis.decision_graph import CoverageDecisionGraph
from rxn_core.retrosynthesis.assembly import AssemblyProblem
from rxn_core.retrosynthesis.ranking import assembly_rank, validate_atom_ownership
from rxn_core.retrosynthesis.catalog_index import exact_source_copy_capacity
from rxn_core.fragment_matching import detect_fragments, FragmentDetectionConfig
from rxn_core.fragment_matching import FragmentCandidate, materialize_target_coverage_orbit
from rxn_core.fragment_matching.graph_ops import fragment_equivalence_classes
from rxn_core.alignment.post_aam import AAMHierarchy, FragmentMatch, AtomPermutation


def test_decision_graph_equals_exhaustive_union_covers_including_redundancy():
    masks = (1, 2, 3, 4, 5, 6, 7)
    expected = set()
    for size in range(1, len(masks) + 1):
        for selected in combinations(masks, size):
            union = 0
            for mask in selected:
                union |= mask
            if union == 7:
                expected.add(frozenset(selected))
    graph = CoverageDecisionGraph.build(masks, 3)
    assert {frozenset(c) for c in graph.covers()} == expected
    assert frozenset((3, 6)) in expected  # Overlap is valid support.
    assert frozenset((1, 7)) in expected  # No invented minimal-cover constraint.
    assert len(graph.nodes) < 2 ** len(masks)


def test_explicit_copy_limit_and_duplicate_slots_equal_exhaustive_oracle():
    masks = (1, 3, 3, 5, 6, 7)
    for maximum in range(1, len(masks) + 1):
        expected = set()
        for size in range(1, maximum + 1):
            for slots in combinations(range(len(masks)), size):
                if CoverageDecisionGraph._union(masks[i] for i in slots) == 7:
                    expected.add(slots)
        graph = CoverageDecisionGraph.build(masks, 3, maximum_regions=maximum)
        assert set(graph.paths()) == expected


def test_pattern_quotient_preserves_fragment_hierarchy_and_target_symmetry():
    from rxn_core.retrosynthesis.assembly import construction_pattern
    labels = ((6, 0, 0),) * 3
    bonds = ((0, 1, 1.), (1, 2, 1.))
    def part(atoms):
        return {"target_fragment_atoms": (atoms,), "preserved_target_bonds": ()}
    left = construction_pattern((part((0,)), part((1, 2))), labels, bonds)
    right = construction_pattern((part((0, 1)), part((2,))), labels, bonds)
    assert left == right  # End-to-end target reflection and source-copy reordering.
    split = {"target_fragment_atoms": ((0,), (1, 2)), "preserved_target_bonds": ()}
    assert construction_pattern((split,), labels, bonds) != left


def item(name, atoms, *, total=8, fragments=None):
    atoms = tuple(atoms)
    fragments = fragments or (tuple(range(len(atoms))),)
    mapping = tuple(enumerate(atoms))
    occupation = {"covered_target_atoms": atoms, "mapping": mapping,
        "retained_fragments": fragments, "attachment_atoms_target": (),
        "target_fragment_atoms": tuple(tuple(atoms[a] for a in f) for f in fragments)}
    return {"precursor_id": name, "structure_key": name, "target_occupations": (occupation,),
        "retained_fragments": fragments, "retained_atom_count": len(atoms),
        "total_atom_count": total, "symmetry_retained_atom_count": len(atoms),
        "retained_heavy_atoms": len(atoms), "total_heavy_atoms": total,
        "symmetry_retained_heavy_atoms": len(atoms), "complete": True,
        "chirality_violations": 0, "boundary_bonds": (), "leftover_atom_count": total - len(atoms),
        "preserved_source_bonds": ()}


def test_best_first_order_equals_full_set_ranking_not_local_pool_ranking():
    target = Chem.MolFromSmiles("CCC")
    items = [item("common", (0,), total=30), item("common", (1, 2), total=30),
             item("cheap-A", (0,), total=1), item("cheap-B", (1, 2), total=2),
             item("overlap", (0, 1)), item("overlap", (1, 2))]
    problem = AssemblyProblem.from_index(SimpleNamespace(groups={0: items}), target)
    exhaustive = sorted(problem.assemblies(), key=assembly_rank)
    ranked = list(problem.ranked_assemblies())
    assert [assembly_rank(a) for a in ranked] == [assembly_rank(a) for a in exhaustive]
    assert set(ranked[0]["precursor_stoichiometry"]) == {"overlap"}
    assert any(a["precursor_stoichiometry"] == {"common": 2} for a in ranked)


def test_equal_coverage_keeps_different_fragment_partitions():
    problem = AssemblyProblem.from_index(SimpleNamespace(groups={0: [
        item("same", (0, 1), fragments=((0, 1),)),
        item("same", (0, 1), fragments=((0,), (1,))),
    ]}), Chem.MolFromSmiles("CC"))
    results = list(problem.ranked_assemblies())
    assert len(results) == 3  # Either relation alone, or both overlapping relations.
    assert results[0]["pattern_key"] != results[1]["pattern_key"]


def test_overlap_does_not_create_an_order_dependent_bond_edit():
    a = {"covered_target_atoms": (0, 1), "preserved_target_bonds": ((0, 1),)}
    b = {"covered_target_atoms": (1, 2), "preserved_target_bonds": ((1, 2),)}
    assert validate_atom_ownership((a, b), ((0, 1), (1, 2))) == []
    assert validate_atom_ownership((b, a), ((0, 1), (1, 2))) == []


def test_typed_assembly_consumes_fragment_bonds_without_reinterpreting_aam():
    from rxn_core.retrosynthesis import assemble_fragment_cover
    target = WeightedGraph(["C", "O"], np.array([[0., 1.], [1., 0.]]))
    candidate = FragmentCandidate("CO", ((0, 0), (1, 1)), (0, 1), (0, 1),
        (), (), (), (), (), 2, retained_fragments=((0, 1),),
        preserved_source_bonds=((0, 1),))
    result = assemble_fragment_cover(target, (candidate,))
    assert result.status == "matched"
    assert result.assemblies[0].formed_bonds == ()


def test_whole_fragment_packing_rejects_an_orbit_count_upper_bound():
    # Opposite pairs of a 6-cycle: the orbit has three disjoint pairs.
    rotation = dict(enumerate((1, 2, 3, 4, 5, 0)))
    capacity, atoms = exact_source_copy_capacity((0, 3), (rotation,))
    assert capacity == 3 and atoms == list(range(6))
    # All rotations of {0,1,3} intersect: atom capacity says 2, exact answer 1.
    capacity, _atoms = exact_source_copy_capacity((0, 1, 3), (rotation,))
    assert capacity == 1


def test_residual_singletons_use_saved_augmented_aam_not_attachment_rules():
    source = WeightedGraph(["C", "O"], np.array([[0., 1.], [1., 0.]]))
    target = WeightedGraph(["C", "O"], np.zeros((2, 2)))
    result = detect_fragments(source, target, source_id="CO", config=FragmentDetectionConfig())
    assert any(c.covered_target_atoms == (0, 1) for c in result.candidates)
    assert result.search_graphs
    assert any(len(g.contexts[0].target_atoms) > 2 for g in result.search_graphs)
    assert all(c.derivations for c in result.candidates)
    assert FragmentDetectionConfig().candidate_limit is None


def test_equivalent_hydrogen_fragments_do_not_expand_factorial_bijections():
    count = 10
    graph = WeightedGraph(["H"] * count, np.zeros((count, count)))
    fragments = tuple((a,) for a in range(count))
    generators = []
    for a in range(count - 1):
        permutation = list(range(count))
        permutation[a], permutation[a + 1] = permutation[a + 1], permutation[a]
        generators.append(AtomPermutation(tuple(permutation)))
    hierarchy = AAMHierarchy((FragmentMatch(0, 0, tuple(range(count)),
        representative_assignments=tuple((a, a) for a in range(count)),
        target_generators=tuple(generators)),))
    candidate = FragmentCandidate("H-units", tuple((a, a) for a in range(count)),
        tuple(range(count)), tuple(range(count)), (), (), (), (), (), count,
        retained_fragments=fragments, aam_hierarchy=hierarchy,
        fragment_classes=fragment_equivalence_classes(graph.to_networkx(), (), fragments, 0.5))
    variants = materialize_target_coverage_orbit(candidate, graph, generators=())
    assert len(variants) == 1  # 10! atom permutations, one occupied fragment relation.


def test_all_reported_co2_fragment_bonds_are_supported_by_target():
    source = WeightedGraph(["O", "C", "O"], np.array([[0., 2., 0.], [2., 0., 2.], [0., 2., 0.]]))
    target = WeightedGraph(["O", "C", "O"], np.array([[0., 2., 0.], [2., 0., 1.], [0., 1., 0.]]))
    detected = detect_fragments(source, target, config=FragmentDetectionConfig())
    assert any(c.covered_target_atoms == (0, 1, 2) for c in detected.candidates)
    for candidate in detected.candidates:
        for variant in materialize_target_coverage_orbit(candidate, target):
            mapping = dict(variant.mapping)
            for a, b in variant.preserved_source_bonds:
                assert abs(source.weights[a, b] - target.weights[mapping[a], mapping[b]]) <= 0.5
