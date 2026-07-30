"""Exact automorphism certificates for hierarchical partial mappings."""
from __future__ import annotations

from collections import defaultdict

from .orbits import _nauty_colored_wbo_graph
from .policy import as_node_match_policy
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
                 node_policy=None):
        self.g_P = g_P
        self.node_policy = as_node_match_policy(node_policy)
        self.nodes, self.atom_index, _base_graph, pair_buckets, zero_bucket = (
            _nauty_colored_wbo_graph(
                g_P,
                wbo_tol=float(getattr(p_orbits, 'wbo_tol', 0.2) or 0.2),
                node_policy=self.node_policy,
            )
        )
        self.n_atoms = len(self.nodes)
        self.atom_base_color = {
            self.atom_index[p]: ('node', self.node_policy.key(g_P, p))
            for p in self.nodes
        }
        adjacency = defaultdict(set)
        edge_vertices_by_bucket = defaultdict(set)
        next_vertex = self.n_atoms
        for (a, b), bucket in sorted(pair_buckets.items()):
            if bucket == zero_bucket:
                continue
            ai = self.atom_index[a]
            bi = self.atom_index[b]
            adjacency[ai].add(next_vertex)
            adjacency[bi].add(next_vertex)
            adjacency[next_vertex].update((ai, bi))
            edge_vertices_by_bucket[bucket].add(next_vertex)
            next_vertex += 1
        self.n_vertices = next_vertex
        self.adjacency = {
            vertex: sorted(adjacency.get(vertex, ()))
            for vertex in range(self.n_vertices)
        }
        self.edge_color_classes = tuple(sorted(
            edge_vertices_by_bucket.items(), key=lambda item: item[0]))

        locked_roles = defaultdict(list)
        for r, p in sorted((int(r), int(p))
                           for r, p in dict(locked_mapping or {}).items()):
            if p in self.atom_index:
                locked_roles[p].append(('locked', r))
        self.locked_roles = {
            p: tuple(roles) for p, roles in locked_roles.items()
        }
    def _candidate_roles(self, cand):
        if isinstance(cand, _SymCand):
            mapping = cand.mapping
            blocks = cand.blocks
        else:
            mapping = dict(cand)
            blocks = ()

        roles = defaultdict(list)
        block_r = {r for block in blocks for r in block.r_atoms}
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

    def _colored_vertices(self, cand):
        candidate_roles = self._candidate_roles(cand)
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

    def graph(self, cand):
        import pynauty

        colored_vertices = self._colored_vertices(cand)
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

    def transporter(self, source, target):
        """Return the exact product-label permutation ``source -> target``."""
        import pynauty

        if self.certificate(source) != self.certificate(target):
            raise ValueError("candidate states are not automorphically equivalent")

        def colors_by_vertex(cand):
            result = {}
            for color, vertices in self._colored_vertices(cand):
                for vertex in vertices:
                    result[vertex] = color
            return result

        source_colors = colors_by_vertex(source)
        target_colors = colors_by_vertex(target)

        # Canonical labels of separately constructed colored graphs may rename
        # whole partition cells.  Build one graph containing both candidates
        # instead.  Each semantic role cell is shared across the two halves,
        # so every automorphism preserves its actual meaning.  The paired hubs
        # make each half one connected component; consequently a generator
        # moving the source hub to the target hub transports the complete
        # candidate state in one exact operation.
        degree = self.n_vertices
        source_hub = 2 * degree
        target_hub = source_hub + 1
        adjacency = defaultdict(set)
        for vertex, neighbors in self.adjacency.items():
            for neighbor in neighbors:
                adjacency[vertex].add(neighbor)
                adjacency[vertex + degree].add(neighbor + degree)
            adjacency[vertex].add(source_hub)
            adjacency[source_hub].add(vertex)
            adjacency[vertex + degree].add(target_hub)
            adjacency[target_hub].add(vertex + degree)

        source_cells = defaultdict(set)
        target_cells = defaultdict(set)
        for vertex, color in source_colors.items():
            source_cells[color].add(vertex)
        for vertex, color in target_colors.items():
            target_cells[color].add(vertex + degree)
        if set(source_cells) != set(target_cells):
            raise ValueError("candidate semantic color domains differ")
        coloring = [
            set(source_cells[color]) | set(target_cells[color])
            for color in sorted(source_cells, key=repr)
        ]
        coloring.append({source_hub, target_hub})
        union_graph = pynauty.Graph(
            2 * degree + 2,
            directed=False,
            adjacency_dict={
                vertex: sorted(adjacency.get(vertex, ()))
                for vertex in range(2 * degree + 2)
            },
            vertex_coloring=coloring,
        )
        raw_generators = pynauty.autgrp(union_graph)[0]
        crossing = next(
            (tuple(map(int, generator)) for generator in raw_generators
             if int(generator[source_hub]) == target_hub),
            None,
        )
        if crossing is None:
            raise ValueError(
                "candidate states share a coarse pynauty certificate but are "
                "not exactly color-preserving equivalent")
        full = tuple(crossing[vertex] - degree for vertex in range(degree))
        if any(image < 0 or image >= degree for image in full):
            raise RuntimeError("candidate union transporter did not swap halves")
        atom_by_index = {index: atom for atom, index in self.atom_index.items()}
        atom_permutation = [0] * self.n_atoms
        for atom, index in self.atom_index.items():
            image_index = full[index]
            if image_index not in atom_by_index:
                raise RuntimeError("candidate transporter mixed atom/edge vertices")
            atom_permutation[int(atom)] = int(atom_by_index[image_index])
        return tuple(atom_permutation)
