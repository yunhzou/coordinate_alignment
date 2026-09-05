from itertools import permutations

from rxn_core.alignment.post_aam import (
    AAMHierarchy, AAMHierarchyView, AAMHierarchyChain, AtomPermutation, FragmentMatch)
from rxn_core.fragment_matching import FragmentCandidate
from rxn_core.fragment_matching.serialization import (
    fragment_candidate_to_record, fragment_candidate_from_record,
    fragment_archive_from_record, repack_fragment_detection_v4,
    repack_fragment_detection_v6,
)


def base():
    generator = AtomPermutation((1, 0, 2, 3))
    return AAMHierarchy(tuple(FragmentMatch(i, i, (i,),
        representative_assignments=((i, i),), target_generators=(generator,)) for i in (0, 1)))


def candidate(hierarchy):
    return FragmentCandidate("source", ((0, 0), (1, 1)), (0, 1), (0, 1),
        (), (), (), (), (), 4, retained_fragments=((0,), (1,)), aam_hierarchy=hierarchy)


def test_lazy_transport_equals_naive_conjugation_for_every_four_atom_action():
    hierarchy = base()
    for images in permutations(range(4)):
        view = hierarchy.relabel_target(dict(enumerate(images)))
        assert isinstance(view, AAMHierarchyView)
        assert "materialized" not in view.__dict__
        expected = [0] * 4
        original = hierarchy.fragments[0].target_generators[0].images
        for a in range(4):
            expected[images[a]] = images[original[a]]
        assert view.fragments[0].target_generators[0].images == tuple(expected)
        assert view.fragments[0].target_generators[0] is view.fragments[1].target_generators[0]
        assert view.fragments[0].representative_assignments == ((0, images[0]),)


def test_composed_views_flatten_and_pad_augmentation_frames():
    hierarchy = base()
    first = dict(enumerate((4, 1, 2, 3, 0)))
    second = dict(enumerate((0, 4, 2, 3, 1)))
    view = hierarchy.relabel_target(first).relabel_target(second)
    assert view.base is hierarchy
    assert "materialized" not in view.__dict__
    combined = {a: second[first[a]] for a in first}
    assert view.to_record() == hierarchy._materialize_target(combined).to_record()
    assert view.fragments[0].target_generators[0].degree == 5


def test_archive_does_not_materialize_views_and_interns_generators():
    view = base().relabel_target({0: 1, 1: 0})
    record = fragment_candidate_to_record(candidate(view))
    assert "materialized" not in view.__dict__
    assert len(record["hierarchy_fragments"]) == 2
    assert len(record["generators"]) == 1
    restored = fragment_candidate_from_record(record)
    assert restored.aam_hierarchy == view
    assert restored.aam_hierarchy.to_record() == view.to_record()
    assert restored.aam_hierarchy.base.fragments[0].target_generators[0] is restored.aam_hierarchy.base.fragments[1].target_generators[0]


def test_v4_repacking_preserves_full_hierarchy_without_search():
    original = candidate(base())
    record = fragment_candidate_to_record(original)
    legacy_candidate = {k: v for k, v in record.items()
                        if k not in ("search_graphs", "generators", "hierarchy_fragments")}
    legacy_candidate["aam_hierarchy"] = original.aam_hierarchy.to_record()
    legacy = {"schema": "rxn_core.fragment_detection/v4", "search_graphs": [],
              "candidates": [legacy_candidate], "source_id": "source"}
    upgraded = repack_fragment_detection_v4(legacy)
    graphs, fragments = fragment_archive_from_record(upgraded)
    restored = fragment_candidate_from_record(upgraded["candidates"][0],
        search_graphs=graphs, hierarchy_fragments=fragments)
    assert restored.aam_hierarchy == original.aam_hierarchy
    assert legacy["candidates"][0]["aam_hierarchy"] == original.aam_hierarchy.to_record()


def test_chain_archive_preserves_independent_frames_without_materialization():
    prefix = base()
    residual = base().relabel_target({0: 1, 1: 0})
    chain = AAMHierarchyChain((prefix, residual))
    record = fragment_candidate_to_record(candidate(chain))
    assert "fragments" not in chain.__dict__
    assert "materialized" not in residual.__dict__
    assert len(record["generators"]) == 1
    assert len(record["hierarchy_fragments"]) == 2
    restored = fragment_candidate_from_record(record).aam_hierarchy
    assert restored == chain
    expected = AAMHierarchy(prefix.fragments + residual.fragments)
    assert restored.to_record() == expected.to_record()
    action = dict(enumerate((4, 1, 2, 3, 0)))
    assert chain.relabel_target(action).to_record() == expected.relabel_target(action).to_record()


def test_v6_repacking_changes_only_reference_envelope():
    original = candidate(base().relabel_target({0: 1, 1: 0}))
    record = fragment_candidate_to_record(original)
    legacy_candidate = {k: v for k, v in record.items()
                        if k not in ("search_graphs", "generators", "hierarchy_fragments")}
    legacy_candidate["aam_hierarchy"] = record["aam_hierarchy"]["segments"][0]
    legacy = {"schema": "rxn_core.fragment_detection/v6", "candidates": [legacy_candidate],
              **{k: record[k] for k in ("search_graphs", "generators", "hierarchy_fragments")}}
    migrated = repack_fragment_detection_v6(legacy)
    graphs, fragments = fragment_archive_from_record(migrated)
    restored = fragment_candidate_from_record(migrated["candidates"][0],
        search_graphs=graphs, hierarchy_fragments=fragments)
    assert restored.aam_hierarchy.to_record() == original.aam_hierarchy.to_record()
