"""Reusable conditional fragment matching over the native-enabled AAM kernel."""
from __future__ import annotations

from dataclasses import dataclass, field
from collections.abc import Mapping
import time

from .growth import IslandBranchLimitExceeded, grow_island
from .matcher import _nauty_orbits


@dataclass(frozen=True)
class FragmentPlacement(Mapping):
    """Matched R–P fragment: witness plus correlated compressed alternatives."""
    mapping: dict
    fragment: frozenset
    deferred_edges: frozenset
    symmetry: dict
    preserved_bonds: tuple

    def __iter__(self):
        return iter(self.mapping)

    def __len__(self):
        return len(self.mapping)

    def __getitem__(self, atom):
        return self.mapping[atom]

    def items(self):
        return self.mapping.items()


@dataclass(frozen=True)
class FragmentMatchContext:
    locked_mapping: dict = field(default_factory=dict)
    islands: dict | None = None
    deferred_edges: tuple = ()
    source_orbits: object = None
    target_orbits: object = None


@dataclass(frozen=True)
class FragmentMatchConfig:
    graph_floor: float = 0.2
    iso_tolerance: float = 0.5
    minimum_size: int = 1
    branch_limit: int = 100
    node_policy: object = None
    allow_mapped_seed: bool = False
    orbit_dedup: bool = True


@dataclass(frozen=True)
class FragmentMatchResult:
    matches: tuple
    capped: bool
    branch_count: int
    branch_limit: int
    elapsed_seconds: float


def match_fragment(source, target, *, seed, context=None, config=None,
                   events=None, profile=None, profile_context=None):
    """Grow one fragment, preserving compressed placements and cap evidence.

    The caller chooses seeds and commits results. This function neither groups
    mechanisms nor requires balanced endpoint composition. Inputs are prepared
    weighted NetworkX graphs, as used by the underlying growth kernel.
    """
    context = context or FragmentMatchContext()
    config = config or FragmentMatchConfig()
    started = time.perf_counter()
    from .subgraph import _coerce_graph
    source = _coerce_graph(source, config.graph_floor)
    target = _coerce_graph(target, config.graph_floor)
    source_orbits, target_orbits = context.source_orbits, context.target_orbits
    if config.orbit_dedup:
        if source_orbits is None:
            source_orbits = _nauty_orbits(source, wbo_tol=config.iso_tolerance,
                                          node_policy=config.node_policy)
        if target_orbits is None:
            target_orbits = _nauty_orbits(target, wbo_tol=config.iso_tolerance,
                                          node_policy=config.node_policy)
    try:
        matches = grow_island(
            source, target, seed, context.locked_mapping,
            graph_floor=config.graph_floor, iso_tol=config.iso_tolerance,
            min_lock_size=config.minimum_size, max_branches=config.branch_limit,
            node_policy=config.node_policy, islands_R=context.islands,
            p_orbits=target_orbits, r_orbits=source_orbits,
            prior_deferred_edges=context.deferred_edges,
            allow_mapped_seed=config.allow_mapped_seed,
            events=events, profile=profile, profile_context=profile_context)
    except IslandBranchLimitExceeded as exc:
        return FragmentMatchResult((), True, exc.count, exc.limit,
                                   time.perf_counter() - started)
    placements = tuple(FragmentPlacement(dict(match), match.fragment,
        match.deferred_edges, match.symmetry,
        tuple(sorted(tuple(sorted((a, b))) for a, b in source.edges()
                     if a in match.fragment and b in match.fragment
                     and tuple(sorted((a, b))) not in match.deferred_edges)))
        for match in matches)
    return FragmentMatchResult(placements, False, len(matches),
                               config.branch_limit, time.perf_counter() - started)
