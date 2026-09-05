"""Assembly over saved occupations; no matching, chemistry, or search budgets."""
from collections import defaultdict
from dataclasses import dataclass
from itertools import product
from itertools import count
from fractions import Fraction
import heapq

from ..search_graph import frozen_value
from .compressed_coverage import place_item
from .decision_graph import CoverageDecisionGraph
from .ranking import build_ranked_assembly, validate_atom_ownership, precursor_cost


def construction_pattern(items, target_labels, target_bonds):
    """Exact colored incidence-graph certificate, up to target symmetry.

    Atom, source-copy, fragment-unit and carried-bond roles remain distinct.
    This is an exact nauty certificate, not a WL hash or a coverage-only key.
    """
    import pynauty
    adjacency, colors = {}, defaultdict(set)
    def vertex(color):
        index = len(adjacency)
        adjacency[index] = set()
        colors[color].add(index)
        return index
    def edge(a, b):
        adjacency[a].add(b)
        adjacency[b].add(a)
    for label in target_labels:
        vertex(("atom", label))
    for a, b, order in target_bonds:
        bond = vertex(("target_bond", order))
        edge(a, bond)
        edge(b, bond)
    for item in items:
        parent = vertex(("source_copy",))
        for part in item["target_fragment_atoms"]:
            fragment = vertex(("matched_fragment",))
            edge(parent, fragment)
            for atom in part:
                edge(fragment, atom)
        for a, b in item["preserved_target_bonds"]:
            bond = vertex(("carried_bond",))
            edge(parent, bond)
            edge(a, bond)
            edge(b, bond)
    graph = pynauty.Graph(len(adjacency), adjacency_dict={a: sorted(b) for a, b in adjacency.items()},
                          vertex_coloring=[colors[c] for c in sorted(colors, key=repr)])
    profile = tuple((repr(c), len(colors[c])) for c in sorted(colors, key=repr))
    return repr(profile) + ":" + pynauty.certificate(graph).hex()


@dataclass(frozen=True)
class AssemblyProblem:
    pools: dict
    decisions: CoverageDecisionGraph
    target_edges: tuple
    target_labels: tuple
    target_bonds: tuple

    @classmethod
    def from_index(cls, index, target):
        pools = defaultdict(dict)
        for items in index.groups.values():
            for item in items:
                for occupation in item["target_occupations"]:
                    placed = place_item(item, occupation)
                    mask = sum(1 << atom for atom in placed["covered_target_atoms"])
                    shape = (mask, frozen_value(sorted(placed["target_fragment_atoms"])),
                             frozen_value(placed["preserved_target_bonds"]),
                             frozen_value(placed["attachment_atoms_target"]))
                    # Retain different source fragments, target partitions and
                    # attachment relations even when their coverage is equal.
                    key = (item["structure_key"], item["precursor_id"],
                           frozen_value(placed["retained_fragments"]),
                           frozen_value(placed["target_fragment_atoms"]),
                           frozen_value(placed["attachment_atoms_target"]),
                           frozen_value(placed["preserved_target_bonds"]),
                           frozen_value(placed["mapping"]))
                    pools[shape].setdefault(key, placed)
        shapes = sorted(pools, key=lambda s: (-s[0].bit_count(), s))
        decisions = CoverageDecisionGraph.build((s[0] for s in shapes), target.GetNumAtoms())
        return cls({index: tuple(pools[shape][key] for key in sorted(pools[shape]))
                    for index, shape in enumerate(shapes)}, decisions,
                   tuple((b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in target.GetBonds()),
                   tuple((a.GetAtomicNum(), a.GetFormalCharge(), a.GetIsotope()) for a in target.GetAtoms()),
                   tuple((b.GetBeginAtomIdx(), b.GetEndAtomIdx(), b.GetBondTypeAsDouble()) for b in target.GetBonds()))

    def assemblies(self):
        """Lazily evaluate every supplier combination for every cover."""
        for slots in self.decisions.paths():
            for items in product(*(self.pools[slot] for slot in slots)):
                connections = validate_atom_ownership(items, self.target_edges)
                if connections is None:
                    raise AssertionError("coverage decision produced an incomplete target")
                assembly = build_ranked_assembly(items, connections)
                assembly["pattern_key"] = construction_pattern(items, self.target_labels, self.target_bonds)
                yield assembly

    def ranked_assemblies(self):
        """Exact best-first traversal with admissible suffix cost bounds.

        Bounds ignore future distinct-source penalties, so they can only be
        optimistic. No state or candidate is removed based on a guessed rank.
        A popped terminal therefore certifies the next globally ranked result.
        """
        graph = self.decisions
        lower = {-2: (Fraction(), 0)}
        for index, node in enumerate(graph.nodes):
            options = []
            if node.skip != -1:
                options.append(lower[node.skip])
            if node.take != -1:
                costs = tuple(precursor_cost(item) for item in self.pools[node.region])
                options.append((min(c[0] for c in costs) + lower[node.take][0],
                                min(c[1] for c in costs) + lower[node.take][1]))
            lower[index] = (min(c[0] for c in options), min(c[1] for c in options))
        queue, serial = [], count()
        atom_count = graph.full.bit_count()
        def push(node_id, items, structures, adjusted, total):
            if node_id == -1:
                return
            bound = lower[node_id]
            ids = tuple(sorted(item["precursor_id"] for item in items)) if node_id == -2 else ()
            rank = (len(structures), -Fraction(atom_count) / (adjusted + bound[0]),
                    -Fraction(atom_count, total + bound[1]), ids)
            heapq.heappush(queue, (rank, next(serial), node_id, items, structures, adjusted, total))
        if graph.root == -1:
            return
        push(graph.root, (), frozenset(), Fraction(), 0)
        while queue:
            _rank, _serial, node_id, items, structures, adjusted, total = heapq.heappop(queue)
            if node_id == -2:
                connections = validate_atom_ownership(items, self.target_edges)
                assembly = build_ranked_assembly(items, connections)
                assembly["pattern_key"] = construction_pattern(items, self.target_labels, self.target_bonds)
                yield assembly
                continue
            node = graph.nodes[node_id]
            push(node.skip, items, structures, adjusted, total)
            for item in self.pools[node.region]:
                cost, atoms = precursor_cost(item)
                push(node.take, items + (item,), structures | {item["structure_key"]},
                     adjusted + cost, total + atoms)
