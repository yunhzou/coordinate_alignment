"""Typed public AAM search API."""
from __future__ import annotations

import time

from .alignment.post_aam import AAMBranch, AAMHierarchy, AtomBijection
from .alignment.sweep import (
    _candidate_from_symmetry_state,
    _completed_branch_with_groups,
    _freeze_analytical,
    attach_completed_candidate_groups,
    completed_candidate_group_generators,
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


_ATTACH_STATE = None


def _attach_init(branches, elements, wbo, graph_floor, wbo_tol):
    """Worker initializer: the branch list and product graph for this pool.

    Under the fork start method the branches are inherited by memory instead
    of being pickled per task; under spawn they are pickled once per worker.
    """
    global _ATTACH_STATE
    graph_product = build_graph(elements, wbo, bond_cut=graph_floor)
    _ATTACH_STATE = (branches, graph_product, wbo_tol)


def _attach_range(bounds):
    """Worker: exact groups for one contiguous range of branches.

    Returns the per-fragment generator sequences of every branch in the range
    and the fragments' cache keys in processing order.  The parent applies the
    generators to its own copies of the branches, so no branch record travels
    back.
    """
    start, stop = bounds
    branches, graph_product, wbo_tol = _ATTACH_STATE
    keys = []
    groups = list(completed_candidate_group_generators(
        branches[start:stop], graph_product, wbo_tol=wbo_tol, key_sink=keys))
    return groups, keys


def _attach_metrics_from_keys(keys):
    """Serial cache accounting replayed over the fragments' cache keys.

    The serial pass keys its generator cache by (locked prefix, frozen
    fragment state) and processes fragments in branch order; walking the same
    keys in the same order reproduces its request/calculation/hit counts.
    """
    metrics = {
        'completed_candidate_group_requests': 0,
        'completed_candidate_group_calculations': 0,
        'completed_candidate_group_cache_hits': 0,
    }
    seen = set()
    for key in keys:
        metrics['completed_candidate_group_requests'] += 1
        if key in seen:
            metrics['completed_candidate_group_cache_hits'] += 1
        else:
            metrics['completed_candidate_group_calculations'] += 1
            seen.add(key)
    return metrics


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
        bounds, start = [], 0
        for chunk in chunks:
            if chunk:
                bounds.append((start, start + len(chunk)))
            start += len(chunk)
        with mp.Pool(len(bounds), initializer=_attach_init,
                     initargs=(raw_branches, tuple(problem.product.elements),
                               problem.product.wbo, config.graph_floor,
                               config.iso_tolerance)) as attach_pool:
            results = attach_pool.map(_attach_range, bounds, chunksize=1)
        groups = [group for chunk_groups, _ in results for group in chunk_groups]
        completed = [
            _completed_branch_with_groups(raw_branch, fragment_generators)
            for raw_branch, fragment_generators in zip(raw_branches, groups)]
        group_metrics = _attach_metrics_from_keys(
            [key for _, keys in results for key in keys])
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


def _branch_from_record(raw, memo=None):
    """Typed branch from a pool record.

    ``memo`` (optional dict) shares the immutable permutation groups,
    permutations and domains that recur across the records of one result;
    equal raw inputs give equal frozen objects, so this changes nothing but
    the construction cost.
    """
    memo = {} if memo is None else memo
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
        group_key = None
        if all(isinstance(generator, (list, tuple))
               for generator in raw_generators):
            group_key = ("grp", len(representative),
                         tuple(tuple(generator) for generator in raw_generators))
            target_group = memo.get(group_key)
        if target_group is None:
            target_group = PermutationGroup.from_generator_mappings(
                len(representative), raw_generators)
            if group_key is not None:
                memo[group_key] = target_group
    return AAMBranch(
        representative=AtomBijection.from_mapping(representative),
        hierarchy=AAMHierarchy.from_record(hierarchy, _memo=memo),
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
    memo = {}
    for key, entry in pool.items():
        representative = AtomBijection.from_mapping(entry["mapping"])
        raw_branches = tuple(entry.get("branches") or ())
        if not raw_branches:
            raise ValueError("AAM mechanism lacks completed branches")
        branches = tuple(_branch_from_record(raw, memo) for raw in raw_branches)
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
