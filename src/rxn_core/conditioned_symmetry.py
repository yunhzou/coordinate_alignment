"""Opt-in native conditioning over one immutable target symmetry graph.

The Python fragment and AAM graph abstractions are unchanged. The native object
reuses role colors and complete ordered partitions, returning the same ordered
generators as the reference finalizer, without permutation enumeration.
"""
from ._engine import ConditionedAutGraph
from .matcher.canonical import _CandidateAutomorphismCanonicalizer
from .search_symmetry import SymmetryWorkspace


class ConditionedSymmetryWorkspace(SymmetryWorkspace):
    def __init__(self, target, iso_tolerance, *, cache_bytes=64 * 1024 * 1024):
        super().__init__(target, iso_tolerance)
        self.canonicalizer = _CandidateAutomorphismCanonicalizer(target, wbo_tol=iso_tolerance)
        c = self.canonicalizer
        self.engine = ConditionedAutGraph(
            c.n_vertices,
            [(a, b) for a, neighbors in c.adjacency.items() for b in neighbors if a <= b],
            list(c.nodes), [c.atom_index[a] for a in c.nodes],
            [repr(c.atom_base_color[c.atom_index[a]]) for a in c.nodes],
            [(repr(('edge', bucket)), sorted(vertices)) for bucket, vertices in c.edge_color_classes],
            cache_bytes,
        )

    def conditioned_generators(self, candidate, locked):
        roles = self.canonicalizer._candidate_roles(candidate, group_domains=True)
        return self.engine.generators(locked, {p: repr(role) for p, role in roles.items()})

    def stats(self):
        return dict(self.engine.stats())
