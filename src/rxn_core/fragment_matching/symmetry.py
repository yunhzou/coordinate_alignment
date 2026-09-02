"""Exact target-group action on compressed fragment candidates."""
from __future__ import annotations

from dataclasses import replace

from ..matcher import _nauty_atom_generators
from ..subgraph import _coerce_graph


class FragmentOrbitLimitExceeded(RuntimeError):
    """The exact target-coverage orbit exceeded its configured limit."""

    def __init__(self, count, limit):
        self.count = int(count)
        self.limit = int(limit)
        super().__init__(
            f"fragment target-coverage orbit cap hit: {count}>{limit}")


def materialize_target_coverage_orbit(
        candidate, target, *, iso_tolerance=0.5, limit=100_000,
        generators=None):
    """Materialize one witness per exact coverage/attachment group orbit.

    Internal atom permutations that leave both target ownership and attachment
    sets unchanged remain compressed.  This is the assembly-relevant quotient;
    it does not enumerate concrete atom bijections.
    """
    if limit < 1:
        raise ValueError("fragment orbit limit must be positive")
    graph = _coerce_graph(target, 0.2)
    generators = tuple(generators or _nauty_atom_generators(
        graph, wbo_tol=iso_tolerance))

    def key(item):
        return (
            tuple(sorted(item.covered_target_atoms)),
            tuple(sorted(item.attachment_atoms_target)),
        )

    seen = {key(candidate): candidate}
    queue = [candidate]
    cursor = 0
    while cursor < len(queue):
        current = queue[cursor]
        cursor += 1
        for generator in generators:
            mapping = tuple(sorted(
                (int(source), int(generator.get(target_atom, target_atom)))
                for source, target_atom in current.mapping
            ))
            transformed = replace(
                current,
                mapping=mapping,
                covered_target_atoms=tuple(sorted(
                    int(generator.get(atom, atom))
                    for atom in current.covered_target_atoms)),
                attachment_atoms_target=tuple(sorted(
                    int(generator.get(atom, atom))
                    for atom in current.attachment_atoms_target)),
            )
            state = key(transformed)
            if state in seen:
                continue
            if len(seen) >= limit:
                raise FragmentOrbitLimitExceeded(len(seen) + 1, limit)
            seen[state] = transformed
            queue.append(transformed)
    return tuple(seen[state] for state in sorted(seen))
