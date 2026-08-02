"""Compile completed AAM branches into exact analytical mapping families."""
from __future__ import annotations

import copy
from collections import defaultdict
import multiprocessing as mp
import time

import numpy as np

from .aam import _branch_from_record
from .alignment.post_aam import AAMHierarchy
from .alignment.sweep import attach_completed_candidate_groups
from .alignment.index_chirality import (
    analytical_family_static_context,
    compile_analytical_mapping_family,
)
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


def _payload_key(record):
    mapping = tuple(sorted(
        (int(source), int(target))
        for source, target in record["mapping"].items()))
    fragments = tuple(sorted(
        tuple(sorted(map(int, fragment.get("fragment") or ())))
        for fragment in (record.get("hierarchy") or {}).get("fragments") or ({
            "fragment": [source for source, _target in mapping],
        },)))
    return mapping, fragments


def _merge_fragment_generators(kept, incoming):
    kept_fragments = {
        tuple(sorted(map(int, item.get("fragment") or ()))): item
        for item in (kept.get("hierarchy") or {}).get("fragments") or ()}
    for item in (incoming.get("hierarchy") or {}).get("fragments") or ():
        key = tuple(sorted(map(int, item.get("fragment") or ())))
        target = kept_fragments.get(key)
        if target is None:
            continue
        left = target.setdefault("symmetry", {})
        right = item.get("symmetry") or {}
        if ("automorph_generators" not in left
                and "automorph_generators" not in right):
            continue
        generators = {
            tuple(map(int, generator))
            for generator in left.get("automorph_generators") or ()}
        generators.update(tuple(map(int, generator)) for generator in
                          right.get("automorph_generators") or ())
        left["automorph_generators"] = [
            list(generator) for generator in sorted(generators)]


_COMPILER_CONTEXT = None


def _init_compiler(context):
    global _COMPILER_CONTEXT
    _COMPILER_CONTEXT = context


def _compile_payload(payload):
    mapping, hierarchy = payload
    context = _COMPILER_CONTEXT
    return compile_analytical_mapping_family(
        mapping, hierarchy,
        context["elements_R"], context["wbo_R"],
        context["elements_P"], context["wbo_P"],
        graph_floor=context["graph_floor"],
        symmetry_wbo_tol=context["symmetry_wbo_tol"],
        dwbo_threshold=context["dwbo_threshold"],
        metal_dwbo_threshold=context["metal_dwbo_threshold"],
        anchor_map=context["anchor_map"],
        static_context=context["static_context"],
    )


def _compile_unique(records, aam, workers, static_context):
    problem, config = aam.problem, aam.config
    context = {
        "elements_R": problem.reactant.elements,
        "wbo_R": problem.reactant.wbo,
        "elements_P": problem.product.elements,
        "wbo_P": problem.product.wbo,
        "graph_floor": config.graph_floor,
        "symmetry_wbo_tol": config.iso_tolerance,
        "dwbo_threshold": config.event_threshold,
        "metal_dwbo_threshold": config.metal_event_threshold,
        "anchor_map": dict(config.anchors),
        "static_context": static_context,
    }
    payloads = [(record["mapping"], record["hierarchy"])
                for record in records]
    count = min(len(payloads), max(1, int(workers)),
                48 if len(payloads) >= 128 else 8)
    if count <= 1 or len(payloads) < 16 or mp.current_process().daemon:
        _init_compiler(context)
        return [_compile_payload(payload) for payload in payloads]
    with mp.get_context("fork").Pool(
            count, initializer=_init_compiler, initargs=(context,)) as pool:
        return pool.map(_compile_payload, payloads)


def _maximal_families(records, aam, workers, static_context):
    """Compile once per exact relation and retain maximal exact cosets."""
    groups = {}
    for source_index, original in enumerate(records):
        record = copy.deepcopy(original)
        key = _payload_key(record)
        provenance = {
            "source_branch_index": source_index,
            "cuts": copy.deepcopy(record.get("cuts") or ()),
            "encounter_count": int(record.get("encounter_count", 1)),
        }
        group = groups.get(key)
        if group is None:
            groups[key] = {"record": record, "provenance": [provenance]}
            continue
        _merge_fragment_generators(group["record"], record)
        group["provenance"].append(provenance)
        group["record"]["encounter_count"] = (
            int(group["record"].get("encounter_count", 1))
            + int(record.get("encounter_count", 1)))
        cuts = {
            tuple(map(int, cut))
            for cut in group["record"].get("cuts") or ()}
        cuts.update(tuple(map(int, cut))
                    for cut in record.get("cuts") or ())
        group["record"]["cuts"] = [list(cut) for cut in sorted(cuts)]
    grouped = list(groups.values())
    compiled = _compile_unique(
        [group["record"] for group in grouped], aam, workers, static_context)
    entries = []
    for group, family in zip(grouped, compiled):
        record = group["record"]
        record["path_provenance"] = group["provenance"]
        record["covered_path_count"] = len(group["provenance"])
        entries.append([record, family])

    def merge(kept, removed):
        kept[0]["encounter_count"] = (
            int(kept[0].get("encounter_count", 1))
            + int(removed[0].get("encounter_count", 1)))
        kept[0].setdefault("path_provenance", []).extend(
            removed[0].get("path_provenance") or ())
        kept[0]["covered_path_count"] = len(kept[0]["path_provenance"])
        cuts = {tuple(map(int, cut))
                for cut in kept[0].get("cuts") or ()}
        cuts.update(tuple(map(int, cut))
                    for cut in removed[0].get("cuts") or ())
        kept[0]["cuts"] = [list(cut) for cut in sorted(cuts)]

    buckets = defaultdict(list)
    unique = []
    for entry in entries:
        equivalent = next((candidate for candidate in
                           buckets[entry[1].equivalence_bucket]
                           if entry[1].equivalent(candidate[1])), None)
        if equivalent is not None:
            merge(equivalent, entry)
        else:
            buckets[entry[1].equivalence_bucket].append(entry)
            unique.append(entry)
    unique.sort(key=lambda entry: -(
        np.log10(entry[1].group_order[0]) + entry[1].group_order[1]))
    maximal = []
    for entry in unique:
        covering = next((candidate for candidate in maximal
                         if entry[1].is_subset_of(candidate[1])), None)
        if covering is not None:
            merge(covering, entry)
        else:
            maximal.append(entry)
    return maximal


def compile_mapping_families(
        aam: AAMResult, *, workers: int = 1,
        minimum_events_only: bool = False) -> AnalyticalAAMResult:
    """Compile and maximally deduplicate exact cosets for AAM mechanisms.

    Geometry, chirality, and RMSD do not participate in this transformation.
    """
    if not isinstance(aam, AAMResult):
        raise TypeError("compile_mapping_families requires an AAMResult")
    started = time.perf_counter()
    problem = aam.problem
    graph_product = build_graph(
        problem.product.elements, problem.product.wbo,
        bond_cut=aam.config.graph_floor)
    static_context = analytical_family_static_context(
        problem.reactant.elements, problem.reactant.wbo,
        problem.product.elements, problem.product.wbo,
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
        maximal = _maximal_families(
            completed, aam, workers, static_context)
        branches = []
        for record, family in maximal:
            typed_branch = _branch_from_record(
                record, record.get("hierarchy") or {})
            branches.append(AnalyticalBranch(typed_branch, family))
        mechanisms.append(AnalyticalMechanism(
            source=mechanism, branches=tuple(branches)))
    return AnalyticalAAMResult(
        aam=aam,
        mechanisms=tuple(mechanisms),
        elapsed_seconds=time.perf_counter() - started)
