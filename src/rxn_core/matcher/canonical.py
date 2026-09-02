"""Exact automorphism certificates for hierarchical partial mappings."""
from __future__ import annotations

from collections import Counter, defaultdict

from .orbits import (
    _nauty_colored_wbo_graph,
    _wbo_tolerance_bucket_lookup,
)
from .policy import (
    AttributeNodeMatchPolicy,
    ElementNodeMatchPolicy,
    as_node_match_policy,
)
from .state import _VERIFY_ROLES, _SymCand, _cand_roles_from_scratch


class _PartialMappingCanonicalizer:
    """Exact joint certificate for a partial source-to-target relation.

    Source and target graphs retain distinct colors.  A subdivision vertex
    represents each mapped pair, so certificate equality proves that exact
    endpoint automorphisms transport one complete partial relation to the
    other.  Atom-orbit labels are never used as a substitute for that proof.
    """

    def __init__(self, g_R, g_P, *, wbo_tol=0.2, node_policy=None,
                 source_atom_tags=None, target_atom_tags=None):
        self.g_R = g_R
        self.g_P = g_P
        self.node_policy = as_node_match_policy(node_policy)
        self.wbo_tol = float(wbo_tol)
        source_atom_tags = dict(source_atom_tags or {})
        target_atom_tags = dict(target_atom_tags or {})
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
        self.wbo_tol = tolerance
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
        """Role dictionary of ``cand`` (see ``_cand_roles_from_scratch``).

        The dictionary depends on the candidate alone, and ``_SymCand`` is
        immutable, so the ``group_domains=False`` result is cached on the
        candidate; child constructors derive it from the parent's cached copy
        when one exists.  ``RXN_CORE_VERIFY_ROLES=1`` recomputes and asserts.
        """
        if group_domains or not isinstance(cand, _SymCand):
            return _cand_roles_from_scratch(cand, group_domains=group_domains)
        roles = cand._roles
        if roles is None:
            roles = _cand_roles_from_scratch(cand)
            cand._roles = roles
        elif _VERIFY_ROLES:
            expected = _cand_roles_from_scratch(cand)
            assert roles == expected, (
                "cached candidate roles differ from the from-scratch roles",
                roles, expected)
        return roles

    def _color_order_key(self, color):
        """repr of a colour key, memoised: cell order must be a function of
        the key alone and the same keys recur across candidates."""
        cache = self.__dict__.setdefault('_color_repr_cache', {})
        key = cache.get(color)
        if key is None:
            key = repr(color)
            cache[color] = key
        return key

    def _colored_vertices_from_roles(self, candidate_roles):
        colors = defaultdict(set)
        locked_roles = self.locked_roles
        atom_base_color = self.atom_base_color
        for p, vertex in self.atom_index.items():
            role = (
                locked_roles.get(p, ()),
                candidate_roles.get(p, ()),
            )
            colors[('atom', atom_base_color[vertex], role)].add(vertex)
        for color_index, vertices in self.edge_color_classes:
            colors[('edge', color_index)].update(vertices)
        order = self._color_order_key
        return tuple(
            (color, frozenset(vertices))
            for color, vertices in sorted(
                colors.items(), key=lambda item: order(item[0])))

    def _colored_vertices(self, cand, *, group_domains=False):
        return self._colored_vertices_from_roles(
            self._candidate_roles(cand, group_domains=group_domains))

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

    def _reusable_graph(self):
        """One pynauty graph over the fixed base adjacency, recoloured per
        certificate.  ``Graph.__init__`` validates the adjacency dictionary on
        every construction; ``set_vertex_coloring`` is the same call the
        constructor makes and replaces the partition completely, so a
        recoloured graph is indistinguishable from a freshly built one."""
        import pynauty
        graph = self.__dict__.get('_reusable_graph_object')
        if graph is None:
            graph = pynauty.Graph(
                self.n_vertices, directed=False,
                adjacency_dict=self.adjacency)
            self._reusable_graph_object = graph
        return graph

    def certificate_from_roles(self, candidate_roles):
        import pynauty
        colored_vertices = self._colored_vertices_from_roles(candidate_roles)
        # pynauty canonicalizes a partition, whose cells are not themselves
        # named.  Preserve the semantic role attached to every cell as part of
        # the coarse certificate; the exact transporter below remains the
        # authoritative equivalence test.
        color_profile = tuple(
            (color, len(vertices)) for color, vertices in colored_vertices)
        graph = self._reusable_graph()
        graph.set_vertex_coloring(
            [set(vertices) for _, vertices in colored_vertices])
        return pynauty.certificate(graph), color_profile

    def certificate(self, cand):
        return self.certificate_from_roles(self._candidate_roles(cand))

    def role_keys_applicable(self, orbits):
        """True when ``orbits`` is the exact orbit map of this target graph.

        The orbit-role key below is valid only for the automorphism orbits of
        the same WBO-coloured graph at the same tolerance and node policy that
        this canonicalizer colours.  The check compares the orbit map's own
        pair-bucket table with a fresh bucket lookup of ``g_P``; the verdict
        is cached on the orbit map per base graph, so the O(N^2) comparison
        runs once per worker.
        """
        from .orbits import _OrbitMap

        if not isinstance(orbits, _OrbitMap) or orbits.wbo_tol is None:
            return False
        if not isinstance(self.node_policy,
                          (ElementNodeMatchPolicy, AttributeNodeMatchPolicy)):
            return False
        verdicts = getattr(orbits, '_role_key_bases', None)
        if verdicts is None:
            verdicts = {}
            orbits._role_key_bases = verdicts
        base_id = id(self.adjacency)
        verdict = verdicts.get(base_id)
        if verdict is None:
            pair_buckets, zero_bucket = _wbo_tolerance_bucket_lookup(
                self.g_P, self.wbo_tol)
            verdict = bool(
                float(orbits.wbo_tol) == self.wbo_tol
                and orbits.zero_bucket == zero_bucket
                and set(orbits) == set(self.atom_index)
                and orbits.wbo_buckets == pair_buckets)
            verdicts[base_id] = verdict
        return verdict

    def role_key(self, cand, orbits):
        """Automorphism-invariant key of one candidate's role colouring.

        Two candidates merge only when a colour-preserving isomorphism maps
        one role colouring onto the other.  Such an isomorphism preserves the
        underlying WBO-coloured target graph, hence is one of its
        automorphisms and maps every atom into its own exact orbit while
        preserving locked and candidate roles.  Equal certificates therefore
        imply equal keys; the key is a necessary condition and never a merge
        rule by itself.  The returned flag is True when every role-carrying
        atom lies in a singleton orbit, in which case an equal key identifies
        the colouring atom by atom and the certificates are equal without
        running nauty.
        """
        return self.role_key_from_roles(self._candidate_roles(cand), orbits)

    def role_key_from_roles(self, candidate_roles, orbits):
        """See :meth:`role_key`; takes an already computed role dictionary."""
        sizes = getattr(orbits, '_orbit_sizes', None)
        if sizes is None:
            sizes = Counter(orbits.values())
            orbits._orbit_sizes = sizes
        locked_roles = self.locked_roles
        singleton = True
        items = []
        for p, role in candidate_roles.items():
            orbit = orbits[p]
            if sizes[orbit] > 1:
                singleton = False
            items.append((orbit, locked_roles.get(p, ()), role))
        # The key is only hashed and compared for equality, so represent the
        # multiset directly instead of sorting nested tuples by repr.
        return frozenset(Counter(items).items()), singleton

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
