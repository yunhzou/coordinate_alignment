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
