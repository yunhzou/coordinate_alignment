"""Compile completed AAM branches into exact analytical mapping families."""
from __future__ import annotations

from dataclasses import dataclass
from types import SimpleNamespace
import time

from .aam import _branch_from_record
from .alignment.post_aam import AAMHierarchy
from .alignment.sweep import attach_completed_candidate_groups
from .domain import (
    AAMResult,
    AnalyticalAAMResult,
    AnalyticalBranch,
    AnalyticalMechanism,
)
from .frag import build_graph


def _domain_record(domain):
    return {
        "r_atoms": list(domain.r_atoms),
        "p_atoms": list(domain.p_atoms),
        "source": domain.source,
        "extendable": bool(domain.extendable),
    }


def _hierarchy_record(hierarchy: AAMHierarchy):
    fragments = []
    for fragment in hierarchy.fragments:
        symmetry = {
            "witness": dict(fragment.representative_assignments),
            "blocks": [
                _domain_record(domain)
                for domain in fragment.symmetry_domains
            ],
            "exact_fixed": list(fragment.exact_fixed),
            "multiplicity": int(fragment.multiplicity),
            "automorph_blocks": [
                _domain_record(domain)
                for domain in fragment.automorph_domains
            ],
        }
        if fragment.target_generators is not None:
            symmetry["automorph_generators"] = [
                list(generator.images)
                for generator in fragment.target_generators
            ]
        fragments.append({
            "fragment_index": int(fragment.fragment_index),
            "island_idx": int(fragment.island_index),
            "fragment": list(fragment.r_atoms),
            "deferred_edges": [
                list(edge) for edge in fragment.deferred_edges],
            "symmetry": symmetry,
        })
    return {
        "rule": "typed_aam_hierarchy",
        "fragments": fragments,
        "blocks": [],
    }


def _branch_record(branch):
    return {
        "mapping": branch.representative.as_dict(),
        "cuts": [list(cut) for cut in branch.cuts],
        "encounter_count": int(branch.encounter_count),
        "covered_path_count": int(branch.covered_path_count),
        "hierarchy": _hierarchy_record(branch.hierarchy),
        "path_provenance": [dict(item) for item in branch.path_provenance],
    }


def _pipeline_input_adapter(aam: AAMResult):
    problem = aam.problem
    return SimpleNamespace(
        elR=list(problem.reactant.elements),
        xyzR=problem.reactant.coordinates,
        wboR=problem.reactant.wbo,
        elP=list(problem.product.elements),
        xyzP=problem.product.coordinates,
        wboP=problem.product.wbo,
    )


def _pipeline_config(aam: AAMResult, workers):
    config = aam.config
    return {
        "graph_floor": config.graph_floor,
        "iso_tol": config.iso_tolerance,
        "dwbo_threshold": config.event_threshold,
        "metal_dwbo_threshold": config.metal_event_threshold,
        "anchor_map": dict(config.anchors),
        "post_aam_workers": max(1, int(workers)),
    }


def compile_mapping_families(
        aam: AAMResult, *, workers: int = 1,
        minimum_events_only: bool = False) -> AnalyticalAAMResult:
    """Compile and maximally deduplicate exact cosets for AAM mechanisms.

    Geometry, chirality, and RMSD do not participate in this transformation.
    """
    if not isinstance(aam, AAMResult):
        raise TypeError("compile_mapping_families requires an AAMResult")
    # Imported here while the old orchestration module is being dismantled;
    # these exact compilers will move intact into this module in the next
    # refactor checkpoint.
    from .pipeline import (
        _dedupe_analytical_mapping_families,
        analytical_family_static_context,
    )

    started = time.perf_counter()
    inputs = _pipeline_input_adapter(aam)
    config = _pipeline_config(aam, workers)
    graph_product = build_graph(
        inputs.elP, inputs.wboP, bond_cut=aam.config.graph_floor)
    static_context = analytical_family_static_context(
        inputs.elR, inputs.wboR, inputs.elP, inputs.wboP,
        graph_floor=aam.config.graph_floor,
        dwbo_threshold=aam.config.event_threshold,
        metal_dwbo_threshold=aam.config.metal_event_threshold)
    sources = (
        aam.minimum_event_mechanisms()
        if minimum_events_only else aam.mechanisms)
    mechanisms = []
    for mechanism in sources:
        raw = [_branch_record(branch) for branch in mechanism.branches]
        completed = attach_completed_candidate_groups(
            raw, graph_product,
            wbo_tol=aam.config.iso_tolerance)
        maximal = _dedupe_analytical_mapping_families(
            inputs, completed, config, static_context=static_context)
        branches = []
        for record in maximal:
            family = record.get("_mapping_family_object")
            if family is None:
                raise RuntimeError(
                    "analytical compiler omitted its exact family object")
            typed_branch = _branch_from_record(
                record, record.get("hierarchy") or {})
            branches.append(AnalyticalBranch(typed_branch, family))
        mechanisms.append(AnalyticalMechanism(
            source=mechanism, branches=tuple(branches)))
    return AnalyticalAAMResult(
        aam=aam,
        mechanisms=tuple(mechanisms),
        elapsed_seconds=time.perf_counter() - started)

