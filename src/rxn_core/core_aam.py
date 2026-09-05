"""Exact partial AAM for mechanism-local endpoint/TS core assignments."""
from __future__ import annotations

import time
from collections import Counter

from .alignment.branch import BranchLimitExceeded, _generate_seed_orders, find_islands
from .alignment.post_aam import AAMHierarchy
from .alignment.sweep import _branch_symmetry_record, _core_mapping_variants
from .domain import (
    AAMSearchConfig,
    AtomAssignment,
    CoreAAMBranch,
    CoreAAMResult,
    MolecularEndpoint,
)
from .frag import build_graph
from .matcher import _nauty_orbits


def search_core_assignments(
        source: MolecularEndpoint, target: MolecularEndpoint, core_atoms, *,
        config: AAMSearchConfig | None = None,
        assignment_limit: int = 20_000) -> CoreAAMResult:
    """Map a source core into a target without losing tuple correlations.

    Search stops once the requested source core has been assigned.  Every
    final branch is expanded only by its candidate-carried exact target
    automorphism action.  Deduplication is then performed on complete core
    tuples, never on independent vertex orbits.
    """
    config = config or AAMSearchConfig()
    core = tuple(sorted(set(map(int, core_atoms))))
    if not core:
        raise ValueError("partial AAM requires a non-empty source core")
    if source.atom_count != target.atom_count:
        raise ValueError("partial AAM currently requires equal atom counts")
    if Counter(source.elements) != Counter(target.elements):
        raise ValueError("partial AAM requires equal endpoint compositions")
    started = time.perf_counter()
    graph_source = build_graph(
        source.elements, source.wbo, bond_cut=config.graph_floor)
    graph_target = build_graph(
        target.elements, target.wbo, bond_cut=config.graph_floor)
    source_orbits = _nauty_orbits(
        graph_source, wbo_tol=config.iso_tolerance)
    target_orbits = _nauty_orbits(
        graph_target, wbo_tol=config.iso_tolerance)

    typed_branches = []
    unique = {}
    capped = 0
    for seed_index, order in enumerate(_generate_seed_orders(
            graph_source, n_trials=config.seed_count)):
        graph = find_islands(
            graph_source, graph_target, list(order),
            iso_tol=config.iso_tolerance,
            max_branches=config.branch_limit,
            dwbo_threshold=config.event_threshold,
            metal_dwbo_threshold=config.metal_event_threshold,
            symmetry_wbo_tol=config.iso_tolerance,
            core_R=core,
            stop_when_core_mapped=True,
            p_orbits=target_orbits,
            r_orbits=source_orbits,
        )
        capped += int(graph.capped)
        for raw in graph.paths():
            if not all(atom in raw.mapping for atom in core):
                continue
            variants = _core_mapping_variants(
                raw, core, assignment_limit,
                g_P=graph_target, p_orbits=target_orbits)
            exact = tuple(AtomAssignment.from_mapping(mapping)
                          for mapping in variants)
            representative = AtomAssignment.from_mapping({
                atom: raw.mapping[atom] for atom in core})
            hierarchy = AAMHierarchy.from_record(
                _branch_symmetry_record(raw))
            typed_branches.append(CoreAAMBranch(
                representative=representative,
                hierarchy=hierarchy,
                exact_assignments=exact,
                seed_index=seed_index,
                search_path=raw,
            ))
            for assignment in exact:
                unique.setdefault(assignment.pairs, assignment)
                if len(unique) > int(assignment_limit):
                    raise BranchLimitExceeded(
                        assignment_limit, branch_count=len(unique),
                        stage="exact_core_assignment_union")

    if capped and not typed_branches:
        raise BranchLimitExceeded(
            config.branch_limit, branch_count=config.branch_limit + 1,
            stage="all_core_seed_paths_capped")
    return CoreAAMResult(
        source=source,
        target=target,
        core_atoms=core,
        branches=tuple(typed_branches),
        assignments=tuple(unique[key] for key in sorted(unique)),
        elapsed_seconds=time.perf_counter() - started,
        capped_seed_count=capped,
    )
