"""Opt-in, largest-block-first recommendation; not the exact full workflow.

Bank queries grow connected fragments only. Selected sources receive full
augmented AAM before new suppliers are introduced. Each branch contains whole
correlated per-copy occupations, never independently unioned atom witnesses.
"""
from dataclasses import dataclass
from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED
from fractions import Fraction
from heapq import heappop, heappush
from itertools import count
from time import perf_counter

import numpy as np

from ..frag import WeightedGraph
from ..fragment_matching import (
    FragmentCandidate, FragmentDetectionConfig, detect_fragments,
    materialize_target_coverage_orbit, prepare_fragment_target,
)
from ..fragment_matching.connected import find_connected_fragments
from ..subgraph import _coerce_graph
from .config import DEFAULT_ISO_TOLERANCE


def _connected_query(source_id, source, context, config):
    started = perf_counter()
    result = find_connected_fragments(source, context, source_id=source_id, config=config)
    return source_id, result, perf_counter() - started


def _query_results(sources, context, config, workers):
    if workers == 1:
        for source_id, source in sources:
            yield _connected_query(source_id, source, context, config)
        return
    # Scheduling bound only: every source is visited, with one in-flight task
    # per worker so a slow first result cannot accumulate the entire bank.
    sources = iter(sources)
    with ProcessPoolExecutor(workers) as pool:
        pending = set()
        for _ in range(workers):
            item = next(sources, None)
            if item is not None:
                pending.add(pool.submit(_connected_query, *item, context, config))
        while pending:
            done, pending = wait(pending, return_when=FIRST_COMPLETED)
            for future in done:
                yield future.result()
                item = next(sources, None)
                if item is not None:
                    pending.add(pool.submit(_connected_query, *item, context, config))


@dataclass(frozen=True)
class BetaPlacement:
    candidate: FragmentCandidate
    # Query-local P index -> original P index. R indices never change.
    target_atoms: tuple[int, ...]
    refined: bool = False

    @property
    def mapping(self):
        return tuple((a, self.target_atoms[p]) for a, p in self.candidate.mapping)

    @property
    def covered_atoms(self):
        return frozenset(p for _, p in self.mapping)

    @property
    def key(self):
        return (self.candidate.source_id, self.mapping,
                self.candidate.retained_fragments, self.refined)


@dataclass(frozen=True)
class BetaRecommendation:
    placements: tuple[BetaPlacement, ...]
    uncovered_target_atoms: tuple[int, ...]


@dataclass(frozen=True)
class BetaResult:
    recommendations: tuple[BetaRecommendation, ...]
    best_partial: BetaRecommendation
    capped_searches: int
    elapsed_seconds: float
    events: tuple[dict, ...]
    # Recommendation order is heuristic even when every target atom is covered.
    exhaustive: bool = False


class FragmentQueryBank:
    """Cache connected gap queries separately from selected-R augmentation.

    The caller supplies explicit-H graphs. `checkpoint` receives typed evidence
    immediately after each source computation, before occupation expansion.
    """

    def __init__(self, sources, target, *, config=None, checkpoint=None, workers=1):
        self.config = config or FragmentDetectionConfig(
            iso_tolerance=DEFAULT_ISO_TOLERANCE, branch_limit=100)
        sources = tuple(sources)
        self.sources = dict(sources)
        if len(self.sources) != len(sources):
            raise ValueError('source IDs must be unique')
        if workers < 1:
            raise ValueError('workers must be positive')
        self.workers = workers
        self.target = _coerce_graph(target, self.config.graph_floor)
        self.context = prepare_fragment_target(self.target, config=self.config)
        self.checkpoint = checkpoint
        self.connected = {}
        self.augmented = {}
        self.capped_searches = 0
        self.events = []

    def _save(self, stage, source_id, region, result, seconds):
        self.capped_searches += bool(result.capped_seed_count)
        event = dict(stage=stage, source_id=source_id, target_atoms=region,
                     seconds=seconds,
                     candidates=len(result.candidates),
                     capped=bool(result.capped_seed_count), complete=result.complete)
        self.events.append(event)
        if self.checkpoint is not None:
            self.checkpoint(event, result)

    def query(self, region):
        """Search only the induced missing target, including singleton H gaps."""
        region = tuple(sorted(region))
        if region in self.connected:
            return self.connected[region]
        graph = WeightedGraph([dict(self.target.nodes[a]) for a in region],
            np.asarray(self.target.graph['wbo_matrix'])[np.ix_(region, region)])
        context = prepare_fragment_target(graph, config=self.config)
        elements = {self.target.nodes[a]['element'] for a in region}
        blocks = {}
        def eligible():
            for source_id, source in self.sources.items():
                source_graph = _coerce_graph(source, self.config.graph_floor)
                # Safe impossibility test, not a minimum-overlap threshold.
                if any(data['element'] in elements for _, data in source_graph.nodes(data=True)):
                    yield source_id, source
        for source_id, result, seconds in _query_results(
                eligible(), context, self.config, self.workers):
            self._save('connected_query', source_id, region, result, seconds)
            for candidate in result.candidates:
                for variant in materialize_target_coverage_orbit(candidate, context.graph,
                        iso_tolerance=self.config.iso_tolerance,
                        generators=context.automorphism_generators):
                    block = BetaPlacement(variant, region)
                    blocks.setdefault(block.key, block)
        self.connected[region] = tuple(sorted(blocks.values(), key=lambda p: p.key))
        return self.connected[region]

    def refine(self, block):
        """Replace one provisional R copy with a whole augmented occupation.

        Preserve its provisional target coverage, but permit a different valid
        internal assignment. Never splice two incompatible R mappings together.
        The full selected-R result is cached across gaps and branch alternatives.
        """
        source_id = block.candidate.source_id
        if source_id not in self.augmented:
            started = perf_counter()
            result = detect_fragments(self.sources[source_id], self.context,
                source_id=source_id, config=self.config)
            region = tuple(range(len(self.target)))
            self._save('selected_augmentation', source_id, region, result, perf_counter() - started)
            blocks = {}
            for candidate in result.candidates:
                for variant in materialize_target_coverage_orbit(candidate, self.target,
                        iso_tolerance=self.config.iso_tolerance,
                        generators=self.context.automorphism_generators):
                    placed = BetaPlacement(variant, region, True)
                    blocks.setdefault(placed.key, placed)
            self.augmented[source_id] = tuple(blocks.values())
        return tuple(candidate for candidate in self.augmented[source_id]
                     if block.covered_atoms <= candidate.covered_atoms)


def recommend_big_blocks(bank, *, recommendations=1):
    """Best-first beta search, with lazy refinement and retained alternatives.

    No beam or hidden reactant-count cap. `recommendations` is an explicit
    output stopping rule, not a claim of exact global ranking. Repeated sources
    and overlapping target coverage are allowed. A stalled branch is retained
    as a partial result while other block choices are explored.
    """
    if recommendations < 1:
        raise ValueError('recommendations must be positive')
    started = perf_counter()
    target_atoms = frozenset(range(len(bank.target)))
    queue, seen, serial = [], set(), count()
    answers = {}
    best = BetaRecommendation((), tuple(sorted(target_atoms)))

    def push(placements):
        placements = tuple(sorted(placements, key=lambda p: p.key))
        key = tuple(p.key for p in placements)
        if key in seen:
            return
        seen.add(key)
        covered = frozenset(a for p in placements for a in p.covered_atoms)
        fragments = sum(len(p.candidate.retained_fragments) for p in placements)
        total_atoms = sum(p.candidate.retained_size + sum(map(len,
            p.candidate.leftover_fragments)) for p in placements)
        rank = (len(target_atoms - covered), fragments,
                len({p.candidate.source_id for p in placements}),
                -Fraction(len(covered), total_atoms or 1))
        heappush(queue, (rank, next(serial), placements, covered))

    push(())
    while queue and len(answers) < recommendations:
        _, _, placements, covered = heappop(queue)
        missing = target_atoms - covered
        pending = next((i for i, p in enumerate(placements) if not p.refined), None)
        if pending is not None:
            for replacement in bank.refine(placements[pending]):
                push(placements[:pending] + (replacement,) + placements[pending + 1:])
            continue
        if len(missing) < len(best.uncovered_target_atoms):
            best = BetaRecommendation(placements, tuple(sorted(missing)))
        if not missing:
            answer = BetaRecommendation(placements, ())
            answers.setdefault(tuple(p.key for p in placements), answer)
            best = answer
            continue
        for block in bank.query(missing):
            if block.covered_atoms & missing:
                push(placements + (block,))

    return BetaResult(tuple(answers.values()), best, bank.capped_searches,
                      perf_counter() - started, tuple(bank.events))
