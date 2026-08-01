"""Exact automorphism certificates for hierarchical partial mappings."""
from __future__ import annotations

from collections import defaultdict

from .orbits import _nauty_colored_wbo_graph
from .policy import (
    AttributeNodeMatchPolicy,
    ElementNodeMatchPolicy,
    as_node_match_policy,
)
from .state import _SymCand


class _CandidateAutomorphismCanonicalizer:
    """Canonicalize candidate roles under exact product automorphisms.

    The full WBO-colored product graph supplies the group.  Locked mappings
    individualize their product images.  Candidate mappings and symmetry pools
    are additional vertex roles.  Equal pynauty certificates therefore mean
    that one exact automorphism transports the complete hierarchical candidate
    state to the other while preserving every possible future product edge.
    """

    def __init__(self, g_P, p_orbits=None, locked_mapping=None,
                 node_policy=None, wbo_tol=None, base_cache=None):
        self.g_P = g_P
        self.node_policy = as_node_match_policy(node_policy)
        tolerance = float(
            wbo_tol if wbo_tol is not None else
            (getattr(p_orbits, 'wbo_tol', 0.2) or 0.2))
        # A cut worker reuses one immutable product graph and orbit object but
        # constructs many canonicalizers as the locked prefix changes.  The
        # subdivision graph is independent of that prefix, so retain it on the
        # graph-specific orbit object.  Only built-in immutable key policies
        # are cached; arbitrary user policies may carry mutable state.
        policy_key = None
        if isinstance(self.node_policy, ElementNodeMatchPolicy):
            policy_key = ('element',)
        elif isinstance(self.node_policy, AttributeNodeMatchPolicy):
            policy_key = ('attributes', self.node_policy.fields)
        cache = base_cache
        cache_key = None
        if (cache is None and p_orbits is not None and policy_key is not None
                and hasattr(p_orbits, '__dict__')):
            cache = getattr(p_orbits, '_candidate_canonical_bases', None)
            if cache is None:
                cache = {}
                p_orbits._candidate_canonical_bases = cache
        if cache is not None and policy_key is None:
            # A mutable/custom policy cannot safely share prepared colors.
            cache = None
        if cache is not None:
            cache_key = (g_P, tolerance, policy_key)
        base = cache.get(cache_key) if cache is not None else None
        if base is None:
            nodes, atom_index, _base_graph, pair_buckets, zero_bucket = (
                _nauty_colored_wbo_graph(
                    g_P, wbo_tol=tolerance, node_policy=self.node_policy))
            n_atoms = len(nodes)
            atom_base_color = {
                atom_index[p]: ('node', self.node_policy.key(g_P, p))
                for p in nodes
            }
            adjacency = defaultdict(set)
            edge_vertices_by_bucket = defaultdict(set)
            next_vertex = n_atoms
            for (a, b), bucket in sorted(pair_buckets.items()):
                if bucket == zero_bucket:
                    continue
                ai = atom_index[a]
                bi = atom_index[b]
                adjacency[ai].add(next_vertex)
                adjacency[bi].add(next_vertex)
                adjacency[next_vertex].update((ai, bi))
                edge_vertices_by_bucket[bucket].add(next_vertex)
                next_vertex += 1
            base = (
                tuple(nodes), atom_index, n_atoms, atom_base_color,
                next_vertex,
                {vertex: sorted(adjacency.get(vertex, ()))
                 for vertex in range(next_vertex)},
                tuple(sorted(edge_vertices_by_bucket.items(),
                             key=lambda item: item[0])),
            )
            if cache is not None:
                cache[cache_key] = base
        (self.nodes, self.atom_index, self.n_atoms, self.atom_base_color,
         self.n_vertices, self.adjacency, self.edge_color_classes) = base

        locked_roles = defaultdict(list)
        for r, p in sorted((int(r), int(p))
                           for r, p in dict(locked_mapping or {}).items()):
            if p in self.atom_index:
                locked_roles[p].append(('locked', r))
        self.locked_roles = {
            p: tuple(roles) for p, roles in locked_roles.items()
        }
    def _candidate_roles(self, cand, *, group_domains=False):
        if isinstance(cand, _SymCand):
            mapping = cand.mapping
            blocks = cand.blocks
        else:
            mapping = dict(cand)
            blocks = ()

        roles = defaultdict(list)
        block_r = {r for block in blocks for r in block.r_atoms}
        if group_domains and isinstance(cand, _SymCand):
            block_r.update(
                r for block in cand.automorph_blocks
                for r in block.r_atoms)
        for r, p in sorted(mapping.items()):
            if r not in block_r:
                roles[int(p)].append(('mapped', int(r)))
        for block in blocks:
            block_role = (
                'pool', tuple(int(r) for r in block.r_atoms),
                bool(block.extendable),
            )
            for p in block.p_atoms:
                roles[int(p)].append(block_role)
        if isinstance(cand, _SymCand):
            for block in cand.automorph_blocks:
                group_role = (
                    'automorph_domain', tuple(int(r) for r in block.r_atoms)
                )
                for p in block.p_atoms:
                    roles[int(p)].append(group_role)
        return {p: tuple(sorted(items, key=repr))
                for p, items in roles.items()}

    def _colored_vertices(self, cand, *, group_domains=False):
        candidate_roles = self._candidate_roles(
            cand, group_domains=group_domains)
        colors = defaultdict(set)
        for p in self.nodes:
            vertex = self.atom_index[p]
            role = (
                self.locked_roles.get(p, ()),
                candidate_roles.get(p, ()),
            )
            colors[('atom', self.atom_base_color[vertex], role)].add(vertex)
        for color_index, vertices in self.edge_color_classes:
            colors[('edge', color_index)].update(vertices)
        return tuple(
            (color, frozenset(vertices))
            for color, vertices in sorted(
                colors.items(), key=lambda item: repr(item[0])))

    def graph(self, cand, *, group_domains=False):
        import pynauty

        colored_vertices = self._colored_vertices(
            cand, group_domains=group_domains)
        return pynauty.Graph(
            self.n_vertices,
            directed=False,
            adjacency_dict=self.adjacency,
            vertex_coloring=[set(vertices)
                             for _, vertices in colored_vertices],
        )

    def certificate(self, cand):
        import pynauty
        colored_vertices = self._colored_vertices(cand)
        # pynauty canonicalizes a partition, whose cells are not themselves
        # named.  Preserve the semantic role attached to every cell as part of
        # the coarse certificate; the exact transporter below remains the
        # authoritative equivalence test.
        color_profile = tuple(
            (color, len(vertices)) for color, vertices in colored_vertices)
        return pynauty.certificate(self.graph(cand)), color_profile

    def atom_generators(self, cand):
        """Exact generators for a bounded completed candidate state."""
        import pynauty

        raw_generators = pynauty.autgrp(
            self.graph(cand, group_domains=True))[0]
        atom_by_index = {
            index: atom for atom, index in self.atom_index.items()}
        identity = tuple(range(self.n_atoms))
        generators = []
        seen = set()
        for raw in raw_generators:
            permutation = list(identity)
            for atom, index in self.atom_index.items():
                image_index = int(raw[index])
                if image_index not in atom_by_index:
                    raise RuntimeError(
                        "candidate automorphism mixed atom/edge vertices")
                permutation[int(atom)] = int(atom_by_index[image_index])
            permutation = tuple(permutation)
            if permutation != identity and permutation not in seen:
                seen.add(permutation)
                generators.append(permutation)
        return tuple(generators)
