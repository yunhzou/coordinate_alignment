"""Typed public AAM search API."""
from __future__ import annotations

import time

from .alignment.post_aam import AAMBranch, AAMHierarchy, AtomBijection
from .alignment.sweep import cut_sweep as _execute_cut_sweep
from .domain import (
    AAMMechanism,
    AAMProblem,
    AAMResult,
    AAMSearchConfig,
    AAMSearchMetrics,
)


def _branch_from_record(raw, fallback_hierarchy):
    mapping_family = dict(raw.get("mapping_family") or {})
    representative = (
        mapping_family.get("representative_mapping")
        or raw.get("mapping"))
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
        hierarchy=AAMHierarchy.from_record(
            raw.get("hierarchy") or fallback_hierarchy),
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
        fallback_hierarchy = entry.get("branch_symmetry") or {}
        branches = tuple(
            _branch_from_record(raw, fallback_hierarchy)
            for raw in entry.get("branches") or ({
                "mapping": entry["mapping"],
                "hierarchy": fallback_hierarchy,
                "encounter_count": entry.get("dedup_count", 1),
            },))
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
    return _result_from_pool(
        problem, config, pool, metrics,
        elapsed_seconds=time.perf_counter() - started)

