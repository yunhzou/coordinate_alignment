"""Search-graph contracts independent of optional event/geometry selection."""
import itertools
import json
import random
from dataclasses import replace

import numpy as np
import pytest

from rxn_core import AAMProblem, AAMSearchConfig, MolecularEndpoint, search_aam
from rxn_core.aam import finalize_graph_symmetry
from rxn_core.artifacts import aam_from_record, aam_record
from rxn_core.alignment.branch import _Branch, find_islands
from rxn_core.analytical import _payload_key, compile_mapping_families
from rxn_core.frag import build_graph
from rxn_core.fragment import FragmentPlacement, FragmentMatchConfig, FragmentMatchContext, match_fragment
from rxn_core.growth.result import IslandBranchLimitExceeded
from rxn_core.matcher import _nauty_orbits
from rxn_core.search_graph import AAMSearchGraph, SearchContext, SearchGraphBuilder


def network(elements, bonds=()):
    weights = np.zeros((len(elements), len(elements)))
    for left, right in bonds:
        weights[left, right] = weights[right, left] = 1
    return build_graph(elements, weights, bond_cut=0.2)


def iso(mapping):
    return FragmentPlacement(dict(mapping), frozenset(mapping), frozenset(),
                             {'witness': mapping, 'blocks': []}, ())


def problem():
    weights = np.array([[0., 1.], [1., 0.]])
    endpoint = MolecularEndpoint(('H', 'H'), [[0, 0, 0], [1, 0, 0]], weights)
    return AAMProblem(endpoint, endpoint, 'hydrogen')


def test_prefix_fork_reconvergence_does_not_copy_or_cross_contexts():
    graph = network(['C'] * 4)
    recorder = SearchGraphBuilder(SearchContext(tuple(graph), tuple(graph), (0, 1, 2, 3)))
    prefix = _Branch(recorder)
    prefix.commit(iso({0: 0}))
    left, right = prefix.fork(), prefix.fork()
    left.commit(iso({1: 1}))
    left.commit(iso({2: 2}))
    left.commit(iso({1: 1, 2: 2}))
    right.commit(iso({1: 1, 2: 2}))
    assert left.islands_R == right.islands_R
    left.merge_exact_paths(right)
    left.commit(iso({3: 3}))
    recorder.stop(left, 'objective_met')
    result = recorder.finish()
    paths = tuple(result.paths())
    assert len(paths) == 2
    assert all(p.transitions[0] == 0 for p in paths)
    assert paths[0].transitions[-1] == paths[1].transitions[-1]
    assert len(result.branches()) == 2  # distinct fragment histories
    assert all(edge.source < edge.target for edge in result.transitions)
    combined = AAMSearchGraph.combine((result, result))
    assert len(tuple(combined.paths())) == 4  # not a Cartesian product
    assert len(combined.branches()) == 2
    assert all(len(b.paths) == 2 for b in combined.branches())


def test_online_admission_reuses_equal_continuation(monkeypatch):
    import rxn_core.fragment as fragment_module
    calls = []
    def grow(source, target, seed, mapping, **kwargs):
        calls.append((seed, tuple(mapping)))
        if seed == 0:
            return [iso({0: 0}), iso({0: 0})]
        return [iso({1: 1})]
    monkeypatch.setattr(fragment_module, 'grow_island', grow)
    graph = network(['C', 'O'])
    result = find_islands(graph, graph, [0, 1], max_branches=1)
    assert calls == [(0, ()), (1, (0,))]
    assert not result.capped
    assert len(tuple(result.paths())) == 2
    assert len(result.branches()) == 1


def test_sibling_commits_share_snapshots_without_mutating_the_parent():
    parent = _Branch(anchor_map={0: 10, 2: 12})
    before = parent.state_key()
    left, right = parent.fork(), parent.fork()
    assert left.islands_R is right.islands_R is parent.islands_R
    left.commit(iso({1: 11}))
    right.commit(iso({1: 13}))
    assert parent.state_key() == before
    assert parent.mapping == {0: 10, 2: 12}
    assert parent.islands_R == {0: 1, 2: 2}
    # Only the source partition is shared: different target witnesses survive.
    assert left.state_key()[1] is right.state_key()[1]
    assert left.state_key() != right.state_key()
    assert left.state_key()[0][0] is before[0][0]
    old_node = left.node
    left.commit(iso({1: 11}))
    assert left.node == old_node  # a no-op doesn't invent a new decision
    assert len(parent.graph.states) == 3


@pytest.mark.parametrize('seed', range(20))
def test_cached_partitions_and_events_match_literal_commit(seed):
    """Compare the optimized commit against the pre-optimization definition."""
    rng = random.Random(seed)
    atoms = list(range(0, 36, 3))
    assignment = dict(zip(atoms, rng.sample(range(50, 150), len(atoms))))
    anchors = dict((a, assignment[a]) for a in rng.sample(atoms, 3))
    branch = _Branch(anchor_map=anchors)
    mapping, islands = dict(anchors), dict(branch.islands_R)
    cuts, next_iid = set(), branch.next_iid
    for _ in range(30):
        selected = rng.sample(atoms, rng.randrange(1, 6))
        pairs = {a: assignment[a] for a in selected}
        new_cuts = frozenset({tuple(sorted(rng.sample(atoms, 2)))})
        match = FragmentPlacement(pairs, frozenset(pairs), new_cuts,
                                  {'witness': pairs, 'blocks': []}, ())
        touched = {islands[a] for a in pairs if a in islands}
        iid = min(touched) if touched else next_iid
        added, relabeled = [], []
        for a, b in pairs.items():
            if a not in mapping:
                mapping[a] = b
                added.append((a, b))
            elif islands[a] != iid:
                relabeled.append((a, islands[a]))
            islands[a] = iid
        for a, label in list(islands.items()):
            if label in touched and label != iid:
                relabeled.append((a, label))
                islands[a] = iid
        groups = {}
        for a, label in islands.items():
            groups.setdefault(label, []).append(a)
        groups = sorted(tuple(sorted(group)) for group in groups.values())
        islands = {a: label for label, group in enumerate(groups, 1) for a in group}
        next_iid = len(groups) + 1
        cuts.update(new_cuts)
        expected_event = dict(type='island_locked', island_idx=iid, pairs=added,
                              merged_with=sorted(touched - {iid}), relabeled=relabeled,
                              mapped_total=len(mapping))
        events = []
        branch.commit(match, events=events)
        assert events == [expected_event]
        assert branch.mapping == mapping
        assert list(branch.islands_R.items()) == list(islands.items())
        assert branch.islands_P == {mapping[a]: label for a, label in islands.items()}
        assert branch.next_iid == next_iid
        assert branch.state_key() == (tuple(sorted(mapping.items())),
                                      tuple(sorted(islands.items())), tuple(sorted(cuts)))


@pytest.mark.parametrize('seed', range(10))
def test_fragment_bonds_match_whole_graph_definition(seed):
    from rxn_core.growth.result import _IsoResult
    rng = random.Random(seed)
    graph = network(['C'] * 18)
    for a in graph:
        for b in range(a, 18):
            if rng.random() < .2:
                graph.add_edge(a, b)
    atoms = frozenset(rng.sample(list(graph), 9))
    cuts = frozenset(tuple(sorted(edge)) for edge in rng.sample(list(graph.edges()), 4))
    raw = _IsoResult({a: a + 100 for a in atoms}, fragment=atoms, deferred_edges=cuts)
    expected = tuple(sorted(tuple(sorted((a, b))) for a, b in graph.edges()
                            if a in atoms and b in atoms and tuple(sorted((a, b))) not in cuts))
    result = FragmentPlacement.from_match(raw, graph)
    assert result.preserved_bonds == expected


def test_symmetry_finalization_preserves_history_but_only_evaluates_result_ancestry():
    target = network(['C'] * 3)
    recorder = SearchGraphBuilder(SearchContext(tuple(target), tuple(target), (0, 1, 2)))
    prefix = _Branch(recorder)
    prefix.commit(iso({0: 0}))
    live, capped = prefix.fork(), prefix.fork()
    live.commit(iso({1: 1}))
    capped.commit(iso({2: 2}))
    recorder.stop(live, 'objective_met')
    recorder.stop(capped, 'capped')
    raw = recorder.finish()
    selected = raw.ancestor_transitions(raw.terminals)
    assert len(selected) == 2
    result, metrics = finalize_graph_symmetry(raw, target, iso_tolerance=.5)
    assert metrics['completed_candidate_group_requests'] == 2
    assert result.states == raw.states and result.stops == raw.stops
    assert len(result.transitions) == len(raw.transitions)
    assert result.capped
    for edge in result.transitions:
        assert (edge.match['symmetry'].get('automorph_group_source')
                == 'conditioned_search_transition') == (edge.id in selected)
    complete, metrics = finalize_graph_symmetry(result, target, iso_tolerance=.5,
                                               states=range(len(result.states)))
    assert metrics['completed_candidate_group_requests'] == 1
    for edge_id in selected:
        assert complete.transitions[edge_id].match == result.transitions[edge_id].match
    again, metrics = finalize_graph_symmetry(complete, target, iso_tolerance=.5)
    assert metrics['completed_candidate_group_requests'] == 0
    assert again == complete
    restored = AAMSearchGraph.from_record(json.loads(json.dumps(result.to_record())))
    assert next(restored.paths()).hierarchy == next(result.paths()).hierarchy


def test_graph_shares_generator_values_and_typed_prefix_groups(monkeypatch):
    import rxn_core.search_symmetry as symmetry
    target = network(['C'] * 3)
    recorder = SearchGraphBuilder(SearchContext(tuple(target), tuple(target), (0, 1, 2)))
    branch = _Branch(recorder)
    branch.commit(iso({0: 0}))
    branch.commit(iso({1: 1}))
    recorder.stop(branch, 'objective_met')
    # Fresh equal arrays from separate conditioned calculations must be interned.
    monkeypatch.setattr(symmetry._CandidateAutomorphismCanonicalizer,
                        'atom_generators', lambda *a, **k: [list(range(3))])
    graph, _ = finalize_graph_symmetry(recorder.finish(), target, iso_tolerance=.5)
    first, second = graph.transitions
    assert first.match['symmetry']['automorph_generators'] is second.match['symmetry']['automorph_generators']
    paths = [next(graph.paths()), next(graph.paths())]
    a, b = paths[0].hierarchy.fragments
    assert a.target_generators[0] is b.target_generators[0]
    assert a.target_generators is paths[1].hierarchy.fragments[0].target_generators
    assert [f.fragment_index for f in paths[0].hierarchy.fragments] == [0, 1]


@pytest.mark.parametrize('stage', ['fragment_growth', 'combined_live_leaves'])
def test_cap_is_recorded_without_claiming_a_match(monkeypatch, stage):
    import rxn_core.fragment as fragment_module
    def grow(source, target, seed, mapping, **kwargs):
        if stage == 'fragment_growth':
            raise IslandBranchLimitExceeded(2, 1, seed=seed)
        return [iso({0: 0}), iso({0: 1})]
    monkeypatch.setattr(fragment_module, 'grow_island', grow)
    graph = network(['C', 'C'])
    result = find_islands(graph, graph, [0], max_branches=1)
    assert not result.terminals
    assert result.capped
    assert all(s.stage == stage and s.count == 2 and s.limit == 1 for s in result.stops)


def test_partial_fragment_api_and_sparse_original_atom_indices():
    source = network(['C', 'O', 'H'], [(0, 1), (0, 2)]).subgraph([0, 2]).copy()
    target = network(['C', 'O', 'H', 'C', 'H'], [(0, 2), (3, 4)]).subgraph([0, 2, 3, 4]).copy()
    config = FragmentMatchConfig()
    context = FragmentMatchContext(source_orbits=_nauty_orbits(source, wbo_tol=.5),
                                   target_orbits=_nauty_orbits(target, wbo_tol=.5))
    result = match_fragment(source, target, seed=0, context=context, config=config)
    assert result.matches and not result.capped
    assert all(set(m) == {0, 2} for m in result.matches)
    graph = AAMSearchGraph.initial_fragment_search(source, target, 0, result, config)
    finalized, _ = finalize_graph_symmetry(graph, target, iso_tolerance=.5)
    restored = AAMSearchGraph.from_record(json.loads(json.dumps(finalized.to_record())))
    for path in restored.paths():
        realization = path.sample(random.Random(7))
        mapping = dict(realization.mapping)
        assert set(mapping) == {0, 2}
        assert target.has_edge(mapping[0], mapping[2])


def test_tiny_exact_symmetry_matches_exhaustive_ring_automorphisms():
    bonds = [(0, 1), (1, 2), (2, 3), (3, 0)]
    graph = network(['C'] * 4, bonds)
    result = find_islands(graph, graph, [0], iso_tol=.5)
    result, _ = finalize_graph_symmetry(result, graph, iso_tolerance=.5)
    path = next(result.paths())
    edge = next(e for e in result.transitions if e.match is not None)
    generators = edge.placement.target_generators
    orbit = {tuple(path.mapping[a] for a in range(4))}
    queue = list(orbit)
    for mapping in queue:
        for generator in generators:
            moved = tuple(generator.images[a] for a in mapping)
            if moved not in orbit:
                orbit.add(moved)
                queue.append(moved)
    exact = {p for p in itertools.permutations(range(4))
             if all(graph.has_edge(p[a], p[b]) for a, b in bonds)}
    assert orbit == exact
    assert len(exact) == 8
    for seed in range(10):
        sample = path.sample(random.Random(seed), steps_per_fragment=4)
        assert tuple(dict(sample.mapping)[a] for a in range(4)) in exact


def test_saved_raw_aam_is_independent_of_postprocessing(tmp_path, monkeypatch):
    import rxn_core.mechanisms as mechanisms
    def forbidden(*args, **kwargs):
        raise AssertionError('raw search must not score mechanisms')
    monkeypatch.setattr(mechanisms, '_score_branch_mapping', forbidden)
    base = problem()
    enriched = AAMProblem(replace(base.reactant, energy=-1.2, metadata={'origin': 'saved-test'}),
                          base.product)
    result = search_aam(enriched, AAMSearchConfig(seed_count=2), intermediate_dir=tmp_path)
    restored = aam_from_record(json.loads((tmp_path / 'aam.json').read_text()))
    assert restored.problem.reactant.energy == -1.2
    assert restored.problem.reactant.metadata['origin'] == 'saved-test'
    assert restored.graph.to_record() == aam_from_record(
        json.loads(json.dumps(aam_record(result)))).graph.to_record()
    assert list(tmp_path.glob('cut_*.json'))
    assert len(restored.branches) == len(result.branches)
    compiled = compile_mapping_families(restored)
    assert compiled.branches
    assert compiled.branches[0].aam_branch.path_provenance[0]['search_paths']


def test_same_source_atoms_are_not_an_exact_relation_dedup_key():
    fragment = {'fragment': [0], 'deferred_edges': [],
                'symmetry': {'witness': {0: 0}, 'blocks': []}}
    left = {'mapping': {0: 0}, 'hierarchy': {'fragments': [fragment]}}
    right = {'mapping': {0: 0}, 'hierarchy': {'fragments': [
        {**fragment, 'symmetry': {**fragment['symmetry'],
          'blocks': [{'r_atoms': [0], 'p_atoms': [0, 1]}]}}]}}
    assert _payload_key(left) != _payload_key(right)


def test_serial_and_workers_produce_identical_search_graphs():
    config = AAMSearchConfig(seed_count=2)
    serial = search_aam(problem(), config, workers=1)
    parallel = search_aam(problem(), config, workers=2)
    assert serial.graph.to_record() == parallel.graph.to_record()


def test_later_generator_is_transported_with_earlier_fragment():
    graph = network(['C'] * 3)
    recorder = SearchGraphBuilder(SearchContext(tuple(graph), tuple(graph), (0, 1, 2)))
    branch = _Branch(recorder)
    for atom, generators in [(0, [[1, 0, 2]]), (1, [[0, 2, 1]]), (2, [])]:
        placement = iso({atom: atom})
        placement.symmetry['automorph_generators'] = generators
        branch.commit(placement)
    recorder.stop(branch, 'objective_met')
    path = next(recorder.finish().paths())
    sample = path.realize({0: (0,), 1: (0,)})
    assert dict(sample.mapping) == {0: 1, 1: 2, 2: 0}
    # The second choice must not move the already selected first assignment.
    assert dict(sample.hierarchy.fragments[0].representative_assignments) == {0: 1}


def test_capped_sibling_does_not_discard_a_successful_path(monkeypatch):
    import rxn_core.fragment as fragment_module
    def grow(source, target, seed, mapping, **kwargs):
        if seed == 0:
            return [iso({0: 0}), iso({0: 1})]
        if mapping[0] == 0:
            raise IslandBranchLimitExceeded(3, 2, seed=seed)
        return [iso({1: 0})]
    monkeypatch.setattr(fragment_module, 'grow_island', grow)
    graph = network(['C', 'C'])
    result = find_islands(graph, graph, [0, 1], max_branches=2)
    assert result.capped
    assert [p.mapping for p in result.paths()] == [{0: 1, 1: 0}]


def test_fragment_archive_and_target_action_preserve_evidence():
    from rxn_core.smiles import smiles_to_weighted_graph
    from rxn_core.fragment_matching import detect_fragments, materialize_target_coverage_orbit
    from rxn_core.fragment_matching.serialization import (
        fragment_candidate_from_record, fragment_candidate_to_record,
        fragment_detection_to_record)
    source = smiles_to_weighted_graph('CO')
    target = smiles_to_weighted_graph('COC')
    result = detect_fragments(source, target)
    candidate = result.candidates[0]
    variants = materialize_target_coverage_orbit(candidate, target)
    assert {dict(v.mapping)[0] for v in variants} == {0, 2}
    for variant in variants:
        assignments = {a: b for fragment in variant.aam_hierarchy.fragments
                        for a, b in fragment.representative_assignments if a in variant.retained_atoms}
        assert assignments == dict(variant.mapping)
        record = json.loads(json.dumps(fragment_candidate_to_record(variant)))
        restored = fragment_candidate_from_record(record)
        assert restored.mapping == variant.mapping
        assert restored.aam_hierarchy == variant.aam_hierarchy
        derivation = restored.derivations[0]
        action = dict(derivation.target_action)
        witness = {a: b for path in derivation.initial_paths[:1] + derivation.residual_paths
                   for a, b in path.mapping.items() if a in variant.retained_atoms}
        assert {a: action.get(b, b) for a, b in witness.items()} == assignments
    record = json.loads(json.dumps(fragment_detection_to_record(
        result, row_index=0, representation='CO')))
    from rxn_core.fragment_matching.serialization import fragment_archive_from_record
    graphs, fragments = fragment_archive_from_record(record)
    assert all('search_graphs' not in c for c in record['candidates'])
    restored = fragment_candidate_from_record(record['candidates'][0], search_graphs=graphs,
                                              hierarchy_fragments=fragments)
    graph_id = record['candidates'][0]['derivations'][0]['initial_paths'][0]['graph']
    assert restored.derivations[0].initial_paths[0].graph is graphs[graph_id]


def test_offline_viewer_uses_saved_graph_without_new_matching(tmp_path, monkeypatch):
    import re
    import shutil
    import subprocess
    import rxn_core.fragment as fragment_module
    from rxn_core.artifacts import write_aam_bundle
    result = search_aam(problem(), AAMSearchConfig(seed_count=1))
    def forbidden(*args, **kwargs):
        raise AssertionError('viewer must not invoke AAM')
    monkeypatch.setattr(fragment_module, 'grow_island', forbidden)
    write_aam_bundle(result, tmp_path)
    page = (tmp_path / 'search.html').read_text()
    assert '__DATA__' not in page and '__LIBRARY__' not in page
    assert '<script src=' not in page
    assert 'Next recorded history' in page
    if shutil.which('node'):
        script = tmp_path / 'viewer.js'
        script.write_text(re.findall(r'<script>(.*?)</script>', page, re.S)[-1])
        subprocess.run(['node', '--check', str(script)], check=True, capture_output=True)
