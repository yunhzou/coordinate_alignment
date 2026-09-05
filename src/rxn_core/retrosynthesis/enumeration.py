"""Exact region-cover enumeration backed by the shared decision graph."""
from dataclasses import dataclass
from .decision_graph import CoverageDecisionGraph


@dataclass(frozen=True)
class CoverageEnumerationConfig:
    maximum_precursors: int | None = None


@dataclass(frozen=True)
class CoverageEnumerationResult:
    patterns: tuple[tuple[int, ...], ...]
    complete: bool
    truncated: bool


def enumerate_coverage_patterns(masks, atom_count, rank_pattern, *, config=None):
    config = config or CoverageEnumerationConfig()
    graph = CoverageDecisionGraph.build(masks, atom_count,
                                        maximum_regions=config.maximum_precursors)
    patterns = tuple(sorted(graph.covers(), key=lambda p: (rank_pattern(p, atom_count), p)))
    return CoverageEnumerationResult(patterns, True, False)
