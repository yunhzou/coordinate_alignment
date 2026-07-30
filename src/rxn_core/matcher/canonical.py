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

    def graph(self, cand):
        import pynauty

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
        coloring = [
            set(vertices)
            for _, vertices in sorted(colors.items(), key=lambda item: repr(item[0]))
        ]
        return pynauty.Graph(
            self.n_vertices,
            directed=False,
            adjacency_dict=self.adjacency,
            vertex_coloring=coloring,
        )

    def certificate(self, cand):
        import pynauty
        return pynauty.certificate(self.graph(cand))

    def transporter(self, source, target):
        """Return the exact product-label permutation ``source -> target``."""
        import pynauty

        graph_source = self.graph(source)
        graph_target = self.graph(target)
        if pynauty.certificate(graph_source) != pynauty.certificate(graph_target):
            raise ValueError("candidate states are not automorphically equivalent")
        source_label = tuple(map(int, pynauty.canon_label(graph_source)))
        target_label = tuple(map(int, pynauty.canon_label(graph_target)))

        def inverse(permutation):
            result = [0] * len(permutation)
            for atom, image in enumerate(permutation):
                result[image] = atom
            return tuple(result)

        source_inverse = inverse(source_label)
        target_inverse = inverse(target_label)
        candidates = (
            tuple(target_inverse[source_label[v]]
                  for v in range(self.n_vertices)),
            tuple(target_label[source_inverse[v]]
                  for v in range(self.n_vertices)),
        )

        def valid(permutation):
            for vertex in range(self.n_vertices):
                if {permutation[n] for n in self.adjacency[vertex]} != set(
                        self.adjacency[permutation[vertex]]):
                    return False
            source_colors = [set(cell) for cell in graph_source.vertex_coloring]
            target_colors = [set(cell) for cell in graph_target.vertex_coloring]
            return [
                {permutation[v] for v in cell} for cell in source_colors
            ] == target_colors

        full = next((permutation for permutation in candidates
                     if valid(permutation)), None)
        if full is None:
            raise RuntimeError("pynauty canonical labels yielded no transporter")
        atom_by_index = {index: atom for atom, index in self.atom_index.items()}
        atom_permutation = [0] * self.n_atoms
        for atom, index in self.atom_index.items():
            image_index = full[index]
            if image_index not in atom_by_index:
                raise RuntimeError("candidate transporter mixed atom/edge vertices")
            atom_permutation[int(atom)] = int(atom_by_index[image_index])
        return tuple(atom_permutation)
