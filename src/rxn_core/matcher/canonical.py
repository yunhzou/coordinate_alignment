"""Exact automorphism certificates for hierarchical partial mappings."""
from __future__ import annotations

from collections import defaultdict, deque

from .orbits import (
    _nauty_atom_generators,
    _nauty_colored_wbo_graph,
    _wbo_tolerance_bucket_lookup,
)
from .policy import (
    AttributeNodeMatchPolicy,
    ElementNodeMatchPolicy,
    as_node_match_policy,
)
from .state import _SymCand


class _PartialMappingCanonicalizer:
    """Exact joint certificate for a partial source-to-target relation.

    Source and target graphs retain distinct colors.  A subdivision vertex
    represents each mapped pair, so certificate equality proves that exact
    endpoint automorphisms transport one complete partial relation to the
    other.  Atom-orbit labels are never used as a substitute for that proof.
    """

    def __init__(self, g_R, g_P, *, wbo_tol=0.2, node_policy=None,
                 source_atom_tags=None, target_atom_tags=None,
                 base_cache=None):
        self.g_R = g_R
        self.g_P = g_P
        self.node_policy = as_node_match_policy(node_policy)
        self.wbo_tol = float(wbo_tol)
        source_atom_tags = dict(source_atom_tags or {})
        target_atom_tags = dict(target_atom_tags or {})
        policy_key = None
        if isinstance(self.node_policy, ElementNodeMatchPolicy):
            policy_key = ('element',)
        elif isinstance(self.node_policy, AttributeNodeMatchPolicy):
            policy_key = ('attributes', self.node_policy.fields)
        cache = base_cache
        if (cache is None and policy_key is not None
                and hasattr(g_R, 'graph')):
            cache = g_R.graph.setdefault(
                '_partial_mapping_canonical_bases', {})
        cache_key = None
        if cache is not None and policy_key is not None:
            cache_key = (
                g_P,
                self.wbo_tol,
                policy_key,
                tuple(sorted(source_atom_tags.items())),
                tuple(sorted(target_atom_tags.items())),
            )
        base = cache.get(cache_key) if cache_key is not None else None
        if base is not None:
            (
                self.r_nodes,
                self.p_nodes,
                self.r_index,
                self.p_index,
                self.atom_vertex_count,
                self.base_vertex_count,
                self.base_adjacency,
                self.base_colors,
            ) = base
            return

        self.r_nodes = tuple(sorted(g_R.nodes()))
        self.p_nodes = tuple(sorted(g_P.nodes()))
        self.r_index = {atom: index
                        for index, atom in enumerate(self.r_nodes)}
        p_offset = len(self.r_nodes)
        self.p_index = {atom: p_offset + index
                        for index, atom in enumerate(self.p_nodes)}
        self.atom_vertex_count = len(self.r_nodes) + len(self.p_nodes)

        adjacency = defaultdict(set)
        colors = defaultdict(set)
        for atom, vertex in self.r_index.items():
            colors[(
                'atom', 'source', self.node_policy.key(g_R, atom),
                source_atom_tags.get(atom),
            )].add(vertex)
        for atom, vertex in self.p_index.items():
            colors[(
                'atom', 'target', self.node_policy.key(g_P, atom),
                target_atom_tags.get(atom),
            )].add(vertex)

        next_vertex = self.atom_vertex_count
        for side, graph, atom_index in (
                ('source', g_R, self.r_index),
                ('target', g_P, self.p_index)):
            pair_buckets, zero_bucket = _wbo_tolerance_bucket_lookup(
                graph, self.wbo_tol)
            for (left, right), bucket in sorted(pair_buckets.items()):
                if bucket == zero_bucket:
                    continue
                vertex = next_vertex
                next_vertex += 1
                left_vertex = atom_index[left]
                right_vertex = atom_index[right]
                adjacency[vertex].update((left_vertex, right_vertex))
                adjacency[left_vertex].add(vertex)
                adjacency[right_vertex].add(vertex)
                colors[('bond', side, bucket)].add(vertex)

        self.base_vertex_count = next_vertex
        self.base_adjacency = {
            vertex: set(adjacency.get(vertex, ()))
            for vertex in range(next_vertex)
        }
        self.base_colors = {
            color: set(vertices) for color, vertices in colors.items()
        }
        if cache_key is not None:
            cache[cache_key] = (
                self.r_nodes,
                self.p_nodes,
                self.r_index,
                self.p_index,
                self.atom_vertex_count,
                self.base_vertex_count,
                self.base_adjacency,
                self.base_colors,
            )

    def certificate(self, mapping):
        """Return the exact endpoint-automorphism certificate of ``mapping``."""
        import pynauty

        pairs = tuple(sorted((int(source), int(target))
                             for source, target in dict(mapping).items()))
        if len({source for source, _target in pairs}) != len(pairs):
            raise ValueError("partial mapping repeats a source atom")
        if len({target for _source, target in pairs}) != len(pairs):
            raise ValueError("partial mapping repeats a target atom")
        adjacency = {
            vertex: set(neighbors)
            for vertex, neighbors in self.base_adjacency.items()
        }
        colors = {
            color: set(vertices)
            for color, vertices in self.base_colors.items()
        }
        next_vertex = self.base_vertex_count
        for source, target in pairs:
            if source not in self.r_index or target not in self.p_index:
                raise ValueError("partial mapping atom lies outside an endpoint")
            vertex = next_vertex
            next_vertex += 1
            source_vertex = self.r_index[source]
            target_vertex = self.p_index[target]
            adjacency[vertex] = {source_vertex, target_vertex}
            adjacency[source_vertex].add(vertex)
            adjacency[target_vertex].add(vertex)
            colors.setdefault(('partial_mapping',), set()).add(vertex)
        graph = pynauty.Graph(
            next_vertex,
            directed=False,
            adjacency_dict={
                vertex: sorted(adjacency.get(vertex, ()))
                for vertex in range(next_vertex)
            },
            vertex_coloring=[
                set(vertices)
                for _color, vertices in sorted(
                    colors.items(), key=lambda item: repr(item[0]))
            ],
        )
        color_profile = tuple(
            (repr(color), len(vertices))
            for color, vertices in sorted(
                colors.items(), key=lambda item: repr(item[0]))
        )
        return pynauty.certificate(graph), color_profile

    @staticmethod
    def _moved_generators(graph, tolerance, node_policy):
        return tuple(
            {
                int(atom): int(image)
                for atom, image in generator.items()
                if atom != image
            }
            for generator in _nauty_atom_generators(
                graph, wbo_tol=tolerance, node_policy=node_policy)
        )

    def _ensure_relation_generators(self):
        generators = getattr(self, "_relation_generators", None)
        if generators is None:
            generators = (
                self._moved_generators(
                    self.g_R, self.wbo_tol, self.node_policy),
                self._moved_generators(
                    self.g_P, self.wbo_tol, self.node_policy),
            )
            self._relation_generators = generators
        return generators

    def _cached_certificate(self, mapping):
        cache = getattr(self, "_equivalence_certificates", None)
        if cache is None:
            cache = {}
            self._equivalence_certificates = cache
        certificate = cache.get(mapping)
        if certificate is None:
            certificate = self.certificate(dict(mapping))
            cache[mapping] = certificate
        return certificate

    def equivalent(self, left, right):
        """Test exact equivalence under independent endpoint automorphisms.

        The nauty generators span the complete source and target
        automorphism groups.  Applying those generators to a partial
        relation and exhausting its orbit is therefore the exact transporter
        test, without rebuilding and canonicalizing both endpoint graphs for
        every relation.
        """
        left = tuple(sorted((int(source), int(target))
                            for source, target in dict(left).items()))
        right = tuple(sorted((int(source), int(target))
                             for source, target in dict(right).items()))
        if left == right:
            return True
        if len(left) != len(right):
            return False
        if (len({source for source, _target in left}) != len(left)
                or len({target for _source, target in left}) != len(left)
                or len({source for source, _target in right}) != len(right)
                or len({target for _source, target in right}) != len(right)):
            raise ValueError("partial mapping must be injective")

        if self.base_vertex_count <= len(left) * len(left):
            return (self._cached_certificate(left)
                    == self._cached_certificate(right))
        source_generators, target_generators = (
            self._ensure_relation_generators())

        queue = deque((left,))
        seen = {left}
        orbit_budget = max(
            1,
            sum(map(len, source_generators))
            + sum(map(len, target_generators)),
        )
        while queue:
            state = queue.popleft()
            source_atoms = {source for source, _target in state}
            target_atoms = {target for _source, target in state}
            for generator in source_generators:
                if source_atoms.isdisjoint(generator):
                    continue
                candidate = tuple(sorted(
                    (generator.get(source, source), target)
                    for source, target in state
                ))
                if candidate == right:
                    return True
                if candidate not in seen:
                    seen.add(candidate)
                    if len(seen) >= orbit_budget:
                        return (self._cached_certificate(left)
                                == self._cached_certificate(right))
                    queue.append(candidate)
            for generator in target_generators:
                if target_atoms.isdisjoint(generator):
                    continue
                candidate = tuple(sorted(
                    (source, generator.get(target, target))
                    for source, target in state
                ))
                if candidate == right:
                    return True
                if candidate not in seen:
                    seen.add(candidate)
                    if len(seen) >= orbit_budget:
                        return (self._cached_certificate(left)
                                == self._cached_certificate(right))
                    queue.append(candidate)
        return False


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
        import pynauty
        self.nauty_graph = pynauty.Graph(
            self.n_vertices,
            directed=False,
            adjacency_dict=self.adjacency,
        )

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

    def _colored_vertices(self, cand, *, group_domains=False, locked_roles=None):
        candidate_roles = self._candidate_roles(
            cand, group_domains=group_domains)
        colors = defaultdict(set)
        locked_roles = self.locked_roles if locked_roles is None else locked_roles
        for p in self.nodes:
            vertex = self.atom_index[p]
            role = (
                locked_roles.get(p, ()),
                candidate_roles.get(p, ()),
            )
            colors[('atom', self.atom_base_color[vertex], role)].add(vertex)
        for color_index, vertices in self.edge_color_classes:
            colors[('edge', color_index)].update(vertices)
        return tuple(
            (color, frozenset(vertices))
            for color, vertices in sorted(
                colors.items(), key=lambda item: repr(item[0])))

    def graph(self, cand, *, group_domains=False, colored_vertices=None):
        if colored_vertices is None:
            colored_vertices = self._colored_vertices(
                cand, group_domains=group_domains)
        self.nauty_graph.set_vertex_coloring([
            set(vertices) for _, vertices in colored_vertices])
        return self.nauty_graph

    def certificate(self, cand):
        colored_vertices, color_profile = self.coloring_profile(cand)
        return self.certificate_from_coloring(
            colored_vertices, color_profile)

    def coloring_profile(self, cand):
        """Return the exact vertex coloring and its cheap role profile."""
        colored_vertices = self._colored_vertices(cand)
        color_profile = tuple(
            (color, len(vertices)) for color, vertices in colored_vertices)
        return colored_vertices, color_profile

    def certificate_from_coloring(self, colored_vertices, color_profile):
        """Canonicalize one already prepared candidate coloring."""
        import pynauty
        return (
            pynauty.certificate(self.graph(
                None, colored_vertices=colored_vertices)),
            color_profile,
        )

    def atom_generators(self, cand, *, colored_vertices=None):
        """Exact generators for a bounded completed candidate state."""
        import pynauty
        from .._group_ops import project_generators

        raw_generators = pynauty.autgrp(
            self.graph(cand, group_domains=True, colored_vertices=colored_vertices))[0]
        return project_generators(raw_generators, tuple(self.atom_index), tuple(self.atom_index.values()))
