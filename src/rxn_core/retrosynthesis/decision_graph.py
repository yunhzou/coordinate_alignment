"""Exact, shared decision graph for fragment-region covers.

Nodes are (next region, covered atoms); equal suffix problems share a node.
No beam, depth, branch, or result budget participates in construction.
The graph represents all sets of supplied occupation slots, including covers
with overlap and redundant regions. Different slots may cover identical atoms:
their internal fragment relations can differ. Repeated species are unrestricted.
Internal atom permutations never enter this search.
"""
from dataclasses import dataclass


@dataclass(frozen=True)
class CoverageDecision:
    region: int
    skip: int
    take: int


@dataclass(frozen=True)
class CoverageDecisionGraph:
    masks: tuple[int, ...]
    nodes: tuple[CoverageDecision, ...]
    root: int
    full: int

    @classmethod
    def build(cls, masks, atom_count, *, maximum_regions=None):
        full = (1 << atom_count) - 1
        masks = tuple(masks)
        if any(mask <= 0 or mask & ~full for mask in masks):
            raise ValueError("occupation contains atoms outside target")
        suffix = [0] * (len(masks) + 1)
        capacities = [()] * (len(masks) + 1)
        for index in range(len(masks) - 1, -1, -1):
            suffix[index] = suffix[index + 1] | masks[index]
            if maximum_regions is not None:
                capacities[index] = tuple(sorted(
                    (*capacities[index + 1], masks[index].bit_count()), reverse=True
                )[:maximum_regions])
        # -1 = impossible; -2 = covered. Build iteratively for large banks.
        memo, nodes = {}, []
        stack = [(0, 0, 0, False)]
        while stack:
            index, covered, used, expanded = stack.pop()
            key = index, covered, used
            if key in memo:
                continue
            if index == len(masks) or used == maximum_regions:
                memo[key] = -2 if covered == full else -1
                continue
            if (covered | suffix[index]) != full:
                memo[key] = -1
                continue
            # Even disjoint largest remaining slots cannot fill this many
            # atoms within an explicit caller copy limit: a proof, not a beam.
            if (maximum_regions is not None and
                    (full ^ covered).bit_count() > sum(capacities[index][:maximum_regions - used])):
                memo[key] = -1
                continue
            skip = index + 1, covered, used
            take = index + 1, covered | masks[index], (used + 1 if maximum_regions is not None else 0)
            if not expanded:
                stack.append((index, covered, used, True))
                stack.append((*skip, False))
                if take != skip:
                    stack.append((*take, False))
                continue
            left = memo[skip]
            right = memo[take]
            if right == -1:
                memo[key] = left
            else:
                memo[key] = len(nodes)
                nodes.append(CoverageDecision(index, left, right))
        return cls(masks, tuple(nodes), memo[0, 0, 0], full)

    def paths(self):
        stack = [(self.root, ())]
        while stack:
            node_id, selected = stack.pop()
            if node_id == -1:
                continue
            if node_id == -2:
                yield selected
                continue
            node = self.nodes[node_id]
            stack.append((node.skip, selected))
            stack.append((node.take, selected + (node.region,)))

    def covers(self):
        for path in self.paths():
            yield tuple(self.masks[index] for index in path)

    @staticmethod
    def _union(masks):
        union = 0
        for mask in masks:
            union |= mask
        return union

    def to_record(self):
        return {"schema": "rxn_core.coverage_decisions/v1",
                "masks": self.masks, "full": self.full, "root": self.root,
                "nodes": [(n.region, n.skip, n.take) for n in self.nodes]}
