import numpy as np
import pytest

from rxn_core import build_graph
from rxn_core.alignment.index_chirality import (
    compile_analytical_mapping_family,
)
from rxn_core.alignment.post_aam import (
    AffineChiralityConstraint,
    AtomBijection,
    AtomPermutation,
    ConstraintEvaluation,
    FixedMappingRMSD,
    OrientedSimplex,
    PermutationGroup,
    PostAAMMechanism,
    PostAAMSelectionProblem,
)


def test_post_aam_mechanism_ignores_concrete_witnesses():
    entry = {
        "mapping": {0: 0, 1: 1, 2: 2},
        "dedup_count": 17,
        "has_no_cut": True,
        "branch_symmetry": {
            "dedup_witness_count": 9,
            "witnesses": [
                {"mapping": {0: 1, 1: 0, 2: 2}},
            ],
            "fragments": [{
                "fragment_index": 0,
                "island_idx": 4,
                "fragment": [0, 1, 2],
                "deferred_edges": [[1, 2]],
                "symmetry": {"blocks": [{
                    "r_atoms": [0, 1],
                    "p_atoms": [0, 1],
                    "source": "exact_automorph_group",
                }]},
            }],
        },
    }

    mechanism = PostAAMMechanism.from_pool_entry(((), ()), entry)

    assert mechanism.representative == AtomBijection((0, 1, 2))
    assert mechanism.raw_branch_count == 17
    assert len(mechanism.hierarchy.fragments) == 1
    assert mechanism.hierarchy.fragments[0].symmetry_domains[0].r_atoms == (0, 1)
    assert not hasattr(mechanism, "witnesses")


def test_bijection_group_action_is_target_o_mapping_o_source_inverse():
    mapping = AtomBijection((0, 1, 2, 3))
    source = AtomPermutation((1, 0, 2, 3))
    target = AtomPermutation((0, 1, 3, 2))

    transformed = mapping.act(source=source, target=target)

    assert transformed.images == (1, 0, 3, 2)


def test_permutation_group_orbits_come_from_generators_not_samples():
    group = PermutationGroup(6, (
        AtomPermutation((1, 0, 2, 3, 4, 5)),
        AtomPermutation((0, 2, 1, 3, 4, 5)),
        AtomPermutation((0, 1, 2, 4, 3, 5)),
    ))

    assert group.orbits() == ((0, 1, 2), (3, 4), (5,))


def test_fixed_mapping_rmsd_never_changes_correspondence():
    reactant = np.array([
        [0.0, 0.0, 0.0],
        [1.2, 0.0, 0.0],
        [0.0, 1.7, 0.0],
        [0.0, 0.0, 2.3],
    ])
    mechanism = PostAAMMechanism(
        mechanism_key=((), ()),
        representative=AtomBijection((0, 1, 2, 3)),
        hierarchy=PostAAMMechanism.from_pool_entry(
            ((), ()), {"mapping": {0: 0, 1: 1, 2: 2, 3: 3}}
        ).hierarchy,
        endpoint_source_symmetry=PermutationGroup.trivial(4),
        endpoint_target_symmetry=PermutationGroup.trivial(4),
    )
    objective = FixedMappingRMSD(reactant, reactant)

    identity = objective.score(mechanism, AtomBijection((0, 1, 2, 3)))
    swapped = objective.score(mechanism, AtomBijection((1, 0, 2, 3)))

    assert identity < 1e-12
    assert swapped > 0.1


def test_constraints_and_objective_are_separate_postprocessing_objects():
    class KeepFirstAtom:
        name = "keep_first_atom"

        def evaluate(self, mechanism, mapping):
            valid = mapping.images[0] == mechanism.representative.images[0]
            return ConstraintEvaluation(self.name, valid, int(not valid))

    coords = np.array([
        [0.0, 0.0, 0.0], [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0], [0.0, 0.0, 1.0],
    ])
    mechanism = PostAAMMechanism.from_pool_entry(
        ((), ()), {"mapping": {0: 0, 1: 1, 2: 2, 3: 3}})
    problem = PostAAMSelectionProblem(
        mechanism, (KeepFirstAtom(),), FixedMappingRMSD(coords, coords))

    accepted = problem.evaluate(AtomBijection((0, 1, 2, 3)))
    rejected = problem.evaluate(AtomBijection((1, 0, 2, 3)))

    assert accepted.valid
    assert not rejected.valid
    assert accepted.objective_name == "fixed_mapping_proper_fit_rmsd"


def test_endpoint_groups_are_computed_from_graphs_not_witness_samples():
    pytest.importorskip("pynauty")
    elements = ["C", "H", "H", "H"]
    wbo = np.zeros((4, 4))
    for atom in (1, 2, 3):
        wbo[0, atom] = wbo[atom, 0] = 1.0
    graph = build_graph(elements, wbo, bond_cut=0.2)
    base = {
        "mapping": {atom: atom for atom in range(4)},
        "dedup_count": 999,
        "branch_symmetry": {
            "witnesses": [{"mapping": {0: 0, 1: 2, 2: 1, 3: 3}}],
            "fragments": [],
        },
    }

    first = PostAAMMechanism.from_aam_graphs(
        ((), ()), base, graph, graph, symmetry_wbo_tolerance=1.0)
    base["branch_symmetry"]["witnesses"] = []
    second = PostAAMMechanism.from_aam_graphs(
        ((), ()), base, graph, graph, symmetry_wbo_tolerance=1.0)

    assert first.endpoint_source_symmetry == second.endpoint_source_symmetry
    assert first.endpoint_target_symmetry == second.endpoint_target_symmetry
    assert first.endpoint_source_symmetry.orbits() == ((0,), (1, 2, 3))
    assert "witnesses" not in first.symmetry_record()


def test_chirality_is_an_independent_mapping_constraint():
    coords = np.array([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ])
    mechanism = PostAAMMechanism.from_pool_entry(
        ((), ()), {"mapping": {0: 0, 1: 1, 2: 2, 3: 3}})
    constraint = AffineChiralityConstraint(
        coords, (OrientedSimplex(0, (1, 2, 3), 1),))

    preserved = constraint.evaluate(mechanism, AtomBijection((0, 1, 2, 3)))
    reversed_result = constraint.evaluate(
        mechanism, AtomBijection((0, 2, 1, 3)))

    assert preserved.valid
    assert not reversed_result.valid
    assert reversed_result.violations == 1


def test_analytical_mapping_family_dedupes_equal_cosets_not_distinct_cosets():
    elements = ["C"] * 4
    wbo = np.zeros((4, 4))
    hierarchy = {
        "fragments": [
            {"fragment_index": 0, "fragment": [0, 1]},
            {"fragment_index": 1, "fragment": [2, 3]},
        ],
    }
    identity = {0: 0, 1: 1, 2: 2, 3: 3}
    within_block = {0: 1, 1: 0, 2: 2, 3: 3}
    cross_block = {0: 0, 1: 2, 2: 1, 3: 3}

    first = compile_analytical_mapping_family(
        identity, hierarchy, elements, wbo, elements, wbo)
    same = compile_analytical_mapping_family(
        within_block, hierarchy, elements, wbo, elements, wbo)
    distinct = compile_analytical_mapping_family(
        cross_block, hierarchy, elements, wbo, elements, wbo)
    unrestricted = compile_analytical_mapping_family(
        identity, {"fragments": [{"fragment": [0, 1, 2, 3]}]},
        elements, wbo, elements, wbo)

    assert first.equivalent(same)
    assert same.equivalent(first)
    assert not first.equivalent(distinct)
    assert first.contains(identity)
    assert first.contains(within_block)
    assert not first.contains(cross_block)
    assert first.is_subset_of(unrestricted)
    assert distinct.is_subset_of(unrestricted)
    assert not unrestricted.is_subset_of(first)
