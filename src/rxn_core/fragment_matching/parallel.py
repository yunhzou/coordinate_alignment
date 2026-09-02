"""Parallel execution policy for exact R–P fragment detection."""
from __future__ import annotations

from dataclasses import dataclass
import multiprocessing as mp

from ..matcher import _nauty_orbits
from .detection import (
    _InitialFamilyAccumulator,
    _detect_fragments_from_initial,
    _grow_initial_seed,
    _initial_seed_order,
    _prepare_fragment_detection,
    detect_fragments,
)
from .models import FragmentDetectionConfig


@dataclass(frozen=True)
class FragmentDetectionExecution:
    seed_workers: int = 1

    def __post_init__(self):
        if self.seed_workers < 1:
            raise ValueError("seed workers must be positive")


_SEED_STATE = None


def _initialize_seed_worker(
        source, target, config, source_orbits, target_orbits):
    global _SEED_STATE
    _SEED_STATE = (
        source, target, config, source_orbits, target_orbits)


def _run_seed(seed):
    source, target, config, source_orbits, target_orbits = _SEED_STATE
    return _grow_initial_seed(
        source, target, seed, config, source_orbits, target_orbits)


def _parallel_initial_fragment_placements(
        source, target, config, target_orbits, target_region_atoms,
        seed_workers):
    if config.seed_mode != "all":
        raise ValueError(
            "parallel seed execution requires the exact all-seed mode")
    source_orbits = _nauty_orbits(
        source, wbo_tol=config.iso_tolerance)
    seed_order, seed_limited = _initial_seed_order(source, config)
    accumulator = _InitialFamilyAccumulator(
        source, target, config, target_region_atoms)
    capped_seed_count = 0
    maximum_branch_count = 0
    candidate_capped = False
    seed_attempt_count = 0
    worker_count = min(int(seed_workers), len(seed_order))

    context = mp.get_context("fork")
    pool = context.Pool(
        worker_count,
        initializer=_initialize_seed_worker,
        initargs=(
            source, target, config, source_orbits, target_orbits),
    )
    try:
        results = pool.imap(_run_seed, seed_order, chunksize=1)
        for placements, capped, branch_count in results:
            seed_attempt_count += 1
            maximum_branch_count = max(
                maximum_branch_count, branch_count)
            if capped:
                capped_seed_count += 1
                continue
            candidate_capped = accumulator.add(placements)
            if candidate_capped:
                pool.terminate()
                break
        else:
            pool.close()
    finally:
        pool.join()

    return (
        tuple(accumulator.families.values()),
        capped_seed_count,
        maximum_branch_count,
        candidate_capped,
        seed_limited,
        seed_attempt_count,
        len(seed_order) - seed_attempt_count,
        False,
    )


def detect_fragments_parallel(
        source, target, *, source_id="",
        config: FragmentDetectionConfig | None = None,
        execution: FragmentDetectionExecution | None = None,
        target_region_atoms=None):
    """Run independent exact R–P seed experiments in worker processes."""
    config = config or FragmentDetectionConfig()
    execution = execution or FragmentDetectionExecution()
    if execution.seed_workers == 1:
        return detect_fragments(
            source,
            target,
            source_id=source_id,
            config=config,
            target_region_atoms=target_region_atoms,
        )
    source_graph, target_context, region = _prepare_fragment_detection(
        source, target, config, target_region_atoms)
    initial_search = _parallel_initial_fragment_placements(
        source_graph,
        target_context.graph,
        config,
        target_context.atom_orbits,
        region,
        execution.seed_workers,
    )
    return _detect_fragments_from_initial(
        source_graph,
        target_context,
        initial_search,
        source_id=str(source_id),
        config=config,
        region=region,
    )
