"""Finalize exact symmetry on recorded fragment transitions, without events."""
from dataclasses import replace
from collections import defaultdict

from .matcher.state import candidate_from_record
from .matcher.canonical import _CandidateAutomorphismCanonicalizer
from .search_graph import frozen_value


def finalize_graph_symmetry(graph, target, *, iso_tolerance):
    """Attach generators once per conditioned transition, independently of events."""
    cache, coloring_cache, edges = {}, {}, []
    # One graph topology for this target, recolored sequentially for each exact
    # conditioned transition. No graph object escapes this finalization pass.
    canonicalizer = None
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
            if canonicalizer is None:
                canonicalizer = _CandidateAutomorphismCanonicalizer(target, wbo_tol=iso_tolerance)
            locked_roles = defaultdict(list)
            for r, p in sorted(locked):
                if p in canonicalizer.atom_index:
                    locked_roles[p].append(('locked', r))
            coloring = canonicalizer._colored_vertices(candidate, group_domains=True,
                locked_roles={p: tuple(roles) for p, roles in locked_roles.items()})
            # Nauty observes the ordered partition, not our descriptive labels.
            coloring_key = tuple(vertices for _label, vertices in coloring)
            if coloring_key not in coloring_cache:
                coloring_cache[coloring_key] = canonicalizer.atom_generators(
                    candidate, colored_vertices=coloring)
            cache[key] = coloring_cache[coloring_key]
        symmetry = {**state, 'automorph_generators': [list(g) for g in cache[key]],
                    'automorph_group_source': 'conditioned_search_transition'}
        edges.append(replace(edge, match={**edge.match, 'symmetry': symmetry}))
    return replace(graph, transitions=tuple(edges)), {
        'completed_candidate_group_requests': requests,
        'completed_candidate_group_calculations': len(cache),
        'completed_candidate_group_cache_hits': requests - len(cache),
        'exact_coloring_group_calculations': len(coloring_cache),
    }
