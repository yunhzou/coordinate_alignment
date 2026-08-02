"""Typed R/P post-processing built exclusively from analytical AAM data."""
from __future__ import annotations

import time

from .alignment.index_chirality import (
    analyze_group_chirality_branch,
    analytical_family_static_context,
    fixed_mapping_aligned_rmsd,
    select_index_chirality_assignment,
)
from .analytical import _hierarchy_record
from .domain import AnalyticalAAMResult, RPMechanism, RPResult
from .frag import classify_bonds
from .alignment.post_aam import AtomBijection


def _select_branch_mapping(result, mechanism):
    problem = result.aam.problem
    config = result.aam.config
    reactant, product = problem.reactant, problem.product
    static_context = analytical_family_static_context(
        reactant.elements, reactant.wbo,
        product.elements, product.wbo,
        graph_floor=config.graph_floor,
        dwbo_threshold=config.event_threshold,
        metal_dwbo_threshold=config.metal_event_threshold)
    family_source_mappings = [
        branch.aam_branch.representative.as_dict()
        for branch in mechanism.branches]
    successes = []
    failures = []
    for branch_index, branch in enumerate(mechanism.branches):
        family = branch.family
        hierarchy = _hierarchy_record(branch.aam_branch.hierarchy)
        event_cosets = family.exact_event_coset_representatives(
            reactant.wbo, product.wbo, reactant.elements,
            dwbo_threshold=config.event_threshold,
            metal_dwbo_threshold=config.metal_event_threshold)
        for coset_index, source_mapping in enumerate(event_cosets):
            try:
                coset_family = (
                    family.with_coset_representative(source_mapping)
                    if dict(family.source_mapping) != dict(source_mapping)
                    else family)
                group_chirality = analyze_group_chirality_branch(
                    source_mapping,
                    reactant.elements, reactant.coordinates, reactant.wbo,
                    product.elements, product.coordinates, product.wbo,
                    graph_floor=config.graph_floor)
                selection = select_index_chirality_assignment(
                    source_mapping, hierarchy,
                    reactant.elements, reactant.coordinates, reactant.wbo,
                    product.elements, product.coordinates, product.wbo,
                    graph_floor=config.graph_floor,
                    symmetry_wbo_tol=config.iso_tolerance,
                    dwbo_threshold=config.event_threshold,
                    metal_dwbo_threshold=config.metal_event_threshold,
                    anchor_map=dict(config.anchors),
                    group_chirality_frames=group_chirality.defined_frames,
                    static_context=static_context,
                    branch_family_mappings=family_source_mappings,
                    aam_family_generators=coset_family.target_generators,
                    compiled_aam_family=coset_family)
                selected = selection.selected_mapping
                rmsd = fixed_mapping_aligned_rmsd(
                    selected, reactant.coordinates, product.coordinates)
                successes.append((
                    round(float(rmsd), 12),
                    tuple(selected[index]
                          for index in range(problem.atom_count)),
                    int(branch_index), int(coset_index),
                    dict(selected), dict(selection.metadata)))
            except Exception as exc:
                from .alignment.index_chirality import IndexChiralityConflict
                if not isinstance(exc, IndexChiralityConflict):
                    raise
                failures.append({
                    "branch_index": int(branch_index),
                    "event_coset_index": int(coset_index),
                    "reason": str(exc),
                    "diagnostics": getattr(exc, "diagnostics", None),
                })
    if not successes:
        from .alignment.index_chirality import IndexChiralityConflict
        raise IndexChiralityConflict(
            "no exact AAM event coset satisfies index chirality",
            diagnostics={"failures": failures})
    successes.sort(key=lambda item: item[:4])
    _rounded, _mapping_key, branch_index, coset_index, mapping, metadata = (
        successes[0])
    metadata["selected_analytical_branch_index"] = int(branch_index)
    metadata["selected_event_coset_index"] = int(coset_index)
    metadata["analytical_branch_count"] = len(mechanism.branches)
    metadata["event_coset_failure_count"] = len(failures)
    return branch_index, mapping, metadata


def select_rp_mappings(
        analytical: AnalyticalAAMResult) -> RPResult:
    """Apply index chirality and fixed-mapping RMSD to exact AAM families."""
    if not isinstance(analytical, AnalyticalAAMResult):
        raise TypeError("select_rp_mappings requires an AnalyticalAAMResult")
    started = time.perf_counter()
    problem = analytical.aam.problem
    config = analytical.aam.config
    selected_mechanisms = []
    for mechanism in analytical.mechanisms:
        branch_index, mapping, chirality = _select_branch_mapping(
            analytical, mechanism)
        broken, formed, _core_r, _core_p = classify_bonds(
            mapping,
            problem.reactant.wbo,
            problem.product.wbo,
            dwbo_threshold=config.event_threshold,
            elements_R=problem.reactant.elements,
            elements_P=problem.product.elements,
            metal_dwbo_threshold=config.metal_event_threshold)
        inverse = {product_atom: reactant_atom
                   for reactant_atom, product_atom in mapping.items()}
        broken_bonds = tuple(
            (int(left), int(right)) for left, right, _wr, _wp in broken)
        formed_bonds = tuple(
            (int(inverse[left]), int(inverse[right]))
            for left, right, _wr, _wp in formed)
        core_atoms = tuple(sorted({
            atom for bond in (*broken_bonds, *formed_bonds) for atom in bond
        }))
        rmsd = fixed_mapping_aligned_rmsd(
            mapping,
            problem.reactant.coordinates,
            problem.product.coordinates)
        selected_mechanisms.append(RPMechanism(
            analytical=mechanism,
            mapping=AtomBijection.from_mapping(mapping),
            broken_bonds=broken_bonds,
            formed_bonds=formed_bonds,
            core_atoms=core_atoms,
            fixed_mapping_rmsd=rmsd,
            chirality=chirality,
            selected_branch_index=branch_index,
        ))
    return RPResult(
        analytical=analytical,
        mechanisms=tuple(selected_mechanisms),
        elapsed_seconds=time.perf_counter() - started)


def align_reaction(problem, *, search_config=None, workers=1,
                   post_workers=None):
    """Convenience composition of typed AAM, family, and R/P stages."""
    from .aam import search_aam
    from .analytical import compile_mapping_families

    aam = search_aam(problem, search_config, workers=workers)
    families = compile_mapping_families(
        aam, workers=post_workers or workers, minimum_events_only=True)
    return select_rp_mappings(families)

