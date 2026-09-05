"""Finalize exact symmetry on recorded fragment transitions, without events."""
from dataclasses import replace

from .matcher.state import candidate_from_record
from .matcher.canonical import _CandidateAutomorphismCanonicalizer
from .search_graph import frozen_value


def finalize_graph_symmetry(graph, target, *, iso_tolerance):
    """Attach generators once per conditioned transition, independently of events."""
    cache, base_cache, edges = {}, {}, []
    requests = 0
    for edge in graph.transitions:
        if edge.match is None:
            edges.append(edge)
            continue
        requests += 1
        state = edge.match['symmetry']
        locked = graph.states[edge.source].mapping
        key = (locked, frozen_value(state))
        if key not in cache:
            candidate = candidate_from_record(state)
            canonicalizer = _CandidateAutomorphismCanonicalizer(
                target, locked_mapping=dict(locked), wbo_tol=iso_tolerance,
                base_cache=base_cache)
            cache[key] = canonicalizer.atom_generators(candidate)
        symmetry = {**state, 'automorph_generators': [list(g) for g in cache[key]],
                    'automorph_group_source': 'conditioned_search_transition'}
        edges.append(replace(edge, match={**edge.match, 'symmetry': symmetry}))
    return replace(graph, transitions=tuple(edges)), {
        'completed_candidate_group_requests': requests,
        'completed_candidate_group_calculations': len(cache),
        'completed_candidate_group_cache_hits': requests - len(cache),
    }
