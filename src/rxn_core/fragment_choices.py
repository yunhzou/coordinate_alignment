"""Conditional greedy fragments and lazy early closures, without bijections."""
from dataclasses import dataclass
import time

from .fragment import FragmentPlacement, FragmentMatchConfig, FragmentMatchContext, FragmentMatchResult
from .growth import native
from .growth.result import _IsoResult
from .matcher import _nauty_orbits


@dataclass(frozen=True)
class FragmentClosure:
    index: int
    size: int
    next_atom: int
    candidates: int
    atoms: frozenset


class FragmentChoices:
    """Retain one native conditional growth; finalize earlier frontiers lazily.

    Contexts are fixed for the lifetime of the session. No source bonds are
    removed. A closure records its frontier bonds as deferred, not broken.
    """
    def __init__(self, source, target, *, seed, context=None, config=None):
        if not native.available():
            raise ValueError('FragmentChoices requires the native engine')
        self.source = source
        self.config = config or FragmentMatchConfig()
        context = context or FragmentMatchContext()
        if self.config.allow_mapped_seed or self.config.node_policy is not None or not self.config.orbit_dedup:
            raise ValueError('Experimental choices require an unmapped seed and native element/orbit policy')
        if context.growth_replay is not None:
            raise ValueError('FragmentChoices does not support external growth replay')
        orbits = context.target_orbits
        if orbits is None:
            orbits = _nauty_orbits(target, wbo_tol=self.config.iso_tolerance)
        self.r = native.source_graph(source)
        self.p = native.target_graph(target, orbits)
        if self.r is None or self.p is None:
            raise ValueError('FragmentChoices requires prepared WBO graphs')
        mapping = [-1] * len(self.r.nodes)
        for r, p in context.locked_mapping.items():
            mapping[self.r.index[r]] = self.p.index[p]
        if mapping[self.r.index[seed]] >= 0:
            raise ValueError('FragmentChoices seed is already mapped')
        started = time.perf_counter()
        self.engine = native._engine.FragmentChoices(self.r.graph, self.p.graph,
            self.r.index[seed], mapping, self.config.graph_floor, self.config.iso_tolerance,
            self.config.branch_limit, [(self.r.index[r], i) for r, i in (context.islands or {}).items()],
            [(self.r.index[a], self.r.index[b]) for a, b in context.deferred_edges])
        self.normal = self._decode(self.engine.normal(), time.perf_counter() - started)
        self.closures = tuple(FragmentClosure(i, size, self.r.nodes[n], count,
                                              frozenset(self.r.nodes[a] for a in atoms))
                              for i, size, n, count, atoms in self.engine.options()
                              if size >= self.config.minimum_size)
        self.frontier_atoms = frozenset().union(*(m.fragment for m in self.normal.matches),
                                                *(c.atoms for c in self.closures))

    def _decode(self, raw, elapsed):
        placements = []
        for iso in raw['isos']:
            match = _IsoResult(
                {self.r.nodes[r]: self.p.nodes[p] for r, p in iso['mapping'].items()},
                deferred_edges=[(self.r.nodes[a], self.r.nodes[b]) for a, b in iso['deferred_edges']],
                fragment=[self.r.nodes[r] for r in iso['fragment']],
                symmetry=native._translate_symmetry(iso['symmetry'], self.r.nodes, self.p.nodes))
            if len(match.fragment) >= self.config.minimum_size:
                placements.append(FragmentPlacement.from_match(match, self.source))
        return FragmentMatchResult(tuple(placements), raw['capped'],
            raw['cap_count'] if raw['capped'] else len(placements), self.config.branch_limit, elapsed)

    def close(self, closure):
        started = time.perf_counter()
        raw = self.engine.close(closure.index)
        return self._decode(raw, time.perf_counter() - started)
