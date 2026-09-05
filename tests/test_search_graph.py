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
from rxn_core.fragment import FragmentMatchConfig, FragmentMatchContext, match_fragment
from rxn_core.growth.result import _IsoResult, IslandBranchLimitExceeded
from rxn_core.matcher import _nauty_orbits
from rxn_core.search_graph import AAMSearchGraph, SearchContext, SearchGraphBuilder


def network(elements, bonds=()):
    weights = np.zeros((len(elements), len(elements)))
    for left, right in bonds:
        weights[left, right] = weights[right, left] = 1
    return build_graph(elements, weights, bond_cut=0.2)


def iso(mapping):
    return _IsoResult(mapping, fragment=mapping,
                      symmetry={'witness': mapping, 'blocks': []})


def problem():
    weights = np.array([[0., 1.], [1., 0.]])
    endpoint = MolecularEndpoint(('H', 'H'), [[0, 0, 0], [1, 0, 0]], weights)
    return AAMProblem(endpoint, endpoint, 'hydrogen')


def test_prefix_fork_reconvergence_does_not_copy_or_cross_contexts():
    graph = network(['C'] * 4)
    recorder = SearchGraphBuilder(SearchContext(tuple(graph), tuple(graph), (0, 1, 2, 3)))
    prefix = _Branch(recorder)
    prefix.commit(iso({0: 0}), graph)
    left, right = prefix.fork(), prefix.fork()
    left.commit(iso({1: 1}), graph)
    left.commit(iso({2: 2}), graph)
    left.commit(iso({1: 1, 2: 2}), graph)
    right.commit(iso({1: 1, 2: 2}), graph)
    assert left.islands_R == right.islands_R
    left.merge_exact_paths(right)
    left.commit(iso({3: 3}), graph)
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
        branch.commit(placement, graph)
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
    graphs = tuple(AAMSearchGraph.from_record(g) for g in record['search_graphs'])
    assert all('search_graphs' not in c for c in record['candidates'])
    restored = fragment_candidate_from_record(record['candidates'][0], search_graphs=graphs)
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
