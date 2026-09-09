"""Experimental dependency-checked cut replay, separate from search policy.

One session belongs to a fixed uncut source and target. It memoizes complete
conditional growth calls and resumes compressed growth at the first topology
read affected by a changed cut. It never selects seeds, cuts, or branches.
Cache eviction affects time only. Production ``search_aam`` does not enable it.
"""
from dataclasses import dataclass

from .growth import native


@dataclass(frozen=True)
class CutReplayView:
    source: object
    target: object
    source_view: object
    target_view: object
    engine: object
    cuts: tuple


class _CutSession:
    def __init__(self, source, target, target_orbits):
        self.source = source
        self.target = target
        self.source_view = native.source_graph(source)
        self.target_view = native.target_graph(target, target_orbits)
        if self.source_view is None or self.target_view is None:
            raise ValueError("cut replay requires native WBO graphs and exact target orbits")

    def for_cut(self, cut):
        index = self.source_view.index
        source = self.source.copy()
        source.remove_edges_from(cut)
        return CutReplayView(source, self.target, self.source_view, self.target_view, self.engine,
            tuple((index[a], index[b]) for a, b in cut))

    def stats(self):
        return dict(self.engine.stats())


class CutReplay(_CutSession):
    def __init__(self, source, target, target_orbits, *, cache_bytes=64*1024*1024,
                 checkpoints=True, reference_only=False):
        super().__init__(source, target, target_orbits)
        self.engine = native._engine.GrowthReplay(
            self.source_view.graph, self.target_view.graph, cache_bytes, checkpoints, reference_only)


class FragmentRepair(_CutSession):
    """Reuse a later fragment only when its full conditional result is unchanged.

    Target occupancy and observed source/boundary inputs are dependencies.
    Changed earlier mapping histories do not themselves invalidate an otherwise
    independent fragment. Search ordering and admission remain the caller's.
    """
    def __init__(self, source, target, target_orbits, *, cache_bytes=64*1024*1024):
        super().__init__(source, target, target_orbits)
        self.engine = native._engine.FragmentRepair(
            self.source_view.graph, self.target_view.graph, cache_bytes)
