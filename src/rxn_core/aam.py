"""Typed public AAM search API."""
from __future__ import annotations

import time

from .alignment.post_aam import AAMBranch, AAMHierarchy, AtomBijection
from .alignment.sweep import (
    _candidate_from_symmetry_state,
    _freeze_analytical,
    attach_completed_candidate_groups,
    cut_sweep as _execute_cut_sweep,
)
from .domain import (
    AAMMechanism,
    AAMProblem,
    AAMResult,
    AAMSearchConfig,
    AAMSearchMetrics,
)
from .frag import build_graph


_ATTACH_PARALLEL_MIN_BRANCHES = 64


def _attach_chunk(args):
    """Worker: finalize one contiguous chunk of branches."""
    branches, elements, wbo, graph_floor, wbo_tol = args
    graph_product = build_graph(elements, wbo, bond_cut=graph_floor)
    return attach_completed_candidate_groups(
        branches, graph_product, wbo_tol=wbo_tol)


def _serial_attach_metrics(raw_branches):
    """Request/calculation/hit counts exactly as the serial pass reports them.

    The serial finalization keys its generator cache by (locked prefix,
    frozen fragment state).  Walking the same keys in the same order without
    running nauty reproduces its accounting, so a chunked parallel run can
    report identical metrics even though every chunk keeps its own cache.
    """
    metrics = {
        'completed_candidate_group_requests': 0,
        'completed_candidate_group_calculations': 0,
        'completed_candidate_group_cache_hits': 0,
    }
    seen = set()
    for raw_branch in raw_branches:
        locked = {}
        hierarchy = raw_branch.get('hierarchy') or {}
        for fragment in hierarchy.get('fragments') or ():
            metrics['completed_candidate_group_requests'] += 1
            state = fragment.get('symmetry') or {}
            candidate = _candidate_from_symmetry_state(state)
            key = (tuple(sorted(locked.items())), _freeze_analytical(state))
            if key in seen:
                metrics['completed_candidate_group_cache_hits'] += 1
            else:
                metrics['completed_candidate_group_calculations'] += 1
                seen.add(key)
            for r, p in candidate.mapping.items():
                locked[int(r)] = int(p)
    return metrics


def _attach_exact_fragment_groups(problem, config, pool, metrics, workers=1):
    """Finalize candidate-carried groups once, as part of AAM output.

    Each branch's generators depend only on that branch (its fragments are
    processed in order with the locked prefix they define), so branches are
    finalized in contiguous chunks across ``workers`` processes when there are
    enough of them; concatenating the chunks reproduces the serial order and
    contents, and the cache metrics are computed by replaying the serial
    cache walk.
    """
    locations, raw_branches = [], []
    for entry in pool.values():
        branches = list(entry.get("branches") or ())
        if not branches:
            raise ValueError("AAM mechanism lacks completed branch records")
        locations.append((entry, len(branches)))
        raw_branches.extend(branches)
    workers = max(1, int(workers))
    if workers > 1 and len(raw_branches) >= _ATTACH_PARALLEL_MIN_BRANCHES:
        import multiprocessing as mp
        chunk_count = min(workers, max(1, len(raw_branches) // 32))
        size, extra = divmod(len(raw_branches), chunk_count)
        chunks, start = [], 0
        for index in range(chunk_count):
            stop = start + size + (1 if index < extra else 0)
            chunks.append(raw_branches[start:stop])
            start = stop
        payload = [
            (chunk, tuple(problem.product.elements), problem.product.wbo,
             config.graph_floor, config.iso_tolerance)
            for chunk in chunks if chunk]
        with mp.Pool(len(payload)) as attach_pool:
            results = attach_pool.map(_attach_chunk, payload, chunksize=1)
        completed = [branch for chunk in results for branch in chunk]
        group_metrics = _serial_attach_metrics(raw_branches)
    else:
        graph_product = build_graph(
            problem.product.elements, problem.product.wbo,
            bond_cut=config.graph_floor)
        completed, group_metrics = attach_completed_candidate_groups(
            raw_branches, graph_product, wbo_tol=config.iso_tolerance,
            return_metrics=True)
    offset = 0
    for entry, count in locations:
        entry["branches"] = completed[offset:offset + count]
        offset += count
    metrics = dict(metrics or {})
    metrics.update(group_metrics)
    return pool, metrics


def _branch_from_record(raw):
    mapping_family = dict(raw.get("mapping_family") or {})
    # The branch representative is the concrete AAM source mapping.  A
    # canonical representative of a later compiled coset belongs to the
    # analytical family object and must not overwrite branch provenance.
    representative = raw.get("mapping")
    if representative is None:
        raise ValueError("AAM branch record lacks its representative mapping")
    hierarchy = raw.get("hierarchy")
    if not hierarchy:
        raise ValueError("AAM branch record lacks its fragment hierarchy")
    raw_generators = raw.get("target_group_generators")
    if raw_generators is None:
        raw_generators = mapping_family.get("target_generators")
    target_group = None
    if raw_generators is not None:
        from .alignment.post_aam import PermutationGroup
        target_group = PermutationGroup.from_generator_mappings(
            len(representative), raw_generators)
    return AAMBranch(
        representative=AtomBijection.from_mapping(representative),
        hierarchy=AAMHierarchy.from_record(hierarchy),
        encounter_count=int(raw.get("encounter_count", 1)),
        cuts=tuple(tuple(map(int, cut)) for cut in raw.get("cuts") or ()),
        covered_path_count=int(raw.get("covered_path_count", 1)),
        mapping_family=mapping_family,
        path_provenance=tuple(
            dict(item) for item in raw.get("path_provenance") or ()),
        target_group=target_group,
    )


def _result_from_pool(problem, config, pool, metrics, elapsed_seconds):
    mechanisms = []
    for key, entry in pool.items():
        representative = AtomBijection.from_mapping(entry["mapping"])
        raw_branches = tuple(entry.get("branches") or ())
        if not raw_branches:
            raise ValueError("AAM mechanism lacks completed branches")
        branches = tuple(_branch_from_record(raw) for raw in raw_branches)
        mechanisms.append(AAMMechanism(
            key=tuple(key),
            representative=representative,
            branches=branches,
            cuts=tuple(entry.get("cuts") or ()),
            includes_uncut_search=bool(entry.get("has_no_cut", False)),
            encounter_count=int(entry.get("dedup_count", 1)),
        ))
    metrics = dict(metrics or {})
    metrics["retained_branch_count"] = sum(
        len(mechanism.branches) for mechanism in mechanisms)
    return AAMResult(
        problem=problem,
        config=config,
        mechanisms=tuple(mechanisms),
        metrics=AAMSearchMetrics.from_record(metrics, elapsed_seconds),
    )


def search_aam(problem: AAMProblem, config: AAMSearchConfig | None = None,
               *, workers: int = 1) -> AAMResult:
    """Search all configured no-cut/one-cut AAM branches.

    This function performs graph search and mechanism classification only.
    It does not apply index chirality, geometry ranking, RMSD selection, TS
    scoring, serialization, or viewer logic.
    """
    if not isinstance(problem, AAMProblem):
        raise TypeError("search_aam requires an AAMProblem")
    config = config or AAMSearchConfig()
    started = time.perf_counter()
    pool, metrics = _execute_cut_sweep(
        problem.reactant.elements,
        problem.reactant.wbo,
        problem.product.elements,
        problem.product.wbo,
        n_workers=max(1, int(workers)),
        cut_floor=config.cut_floor,
        graph_floor=config.graph_floor,
        iso_tol=config.iso_tolerance,
        dwbo_threshold=config.event_threshold,
        metal_dwbo_threshold=config.metal_event_threshold,
        symmetry_wbo_tol=config.iso_tolerance,
        n_seeds=config.seed_count,
        max_branches=config.branch_limit,
        chunksize=config.task_chunksize,
        symmetry_repair=config.symmetry_repair,
        symmetry_repair_min_changes=config.symmetry_repair_min_changes,
        symmetry_repair_max_evals=(
            config.symmetry_repair_max_evaluations),
        anchor_map=dict(config.anchors),
        return_metrics=True,
    )
    pool, metrics = _attach_exact_fragment_groups(
        problem, config, pool, metrics, workers=max(1, int(workers)))
    return _result_from_pool(
        problem, config, pool, metrics,
        elapsed_seconds=time.perf_counter() - started)
