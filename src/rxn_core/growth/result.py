"""Result containers for fragment growth."""
from __future__ import annotations


class _IsoResult(dict):
    __slots__ = ('deferred_edges', 'fragment', 'symmetry')

    def __init__(self, mapping=None, deferred_edges=(), fragment=(),
                 symmetry=None):
        super().__init__(mapping or {})
        self.deferred_edges = frozenset(tuple(sorted(e)) for e in deferred_edges)
        self.fragment = frozenset(fragment)
        self.symmetry = symmetry or {}


class IslandBranchLimitExceeded(RuntimeError):
    """One fragment-growth seed exceeded its live canonical branch limit."""

    def __init__(self, count, limit, *, seed=None):
        self.count = int(count)
        self.limit = int(limit)
        self.seed = seed
        message = f"fragment subtree branch cap hit: {count}>{limit}"
        if seed is not None:
            message += f" seed={seed}"
        super().__init__(message)
