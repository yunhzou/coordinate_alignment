"""Compile completed AAM branches into exact analytical mapping families."""
from __future__ import annotations

import copy
from collections import defaultdict
import multiprocessing as mp
import time

import numpy as np

from .mechanisms import _branch_from_record
from .alignment.post_aam import AAMHierarchy
from .search_graph import frozen_value
from .alignment.index_chirality import (
    analytical_family_static_context,
    compile_analytical_mapping_family,
)
from .domain import (
    AAMResult,
    MechanismResult,
    MappingFamilyResult,
    AnalyticalAAMResult,
    AnalyticalBranch,
    AnalyticalMechanism,
)


def _hierarchy_record(hierarchy: AAMHierarchy):
    return hierarchy.to_record()


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
    fragments = []
    for fragment in record["hierarchy"]["fragments"]:
        symmetry = {key: value for key, value in fragment['symmetry'].items()
                    if key not in {'multiplicity', 'automorph_group_source'}}
        fragments.append((fragment['fragment'], fragment['deferred_edges'], symmetry))
    return mapping, frozen_value(fragments)


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
        include_event_relations=context["include_event_relations"],
    )


def _compile_unique(records, aam, workers, static_context, include_event_relations):
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
        "include_event_relations": include_event_relations,
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


def _maximal_families(records, aam, workers, static_context, *, include_event_relations=True):
    """Compile once per exact relation and retain maximal exact cosets."""
    groups = {}
    for source_index, original in enumerate(records):
        record = copy.deepcopy(original)
        key = _payload_key(record)
        provenance = {
            "source_branch_index": source_index,
            "cuts": copy.deepcopy(record.get("cuts") or ()),
            "encounter_count": int(record.get("encounter_count", 1)),
            "search_paths": copy.deepcopy(record.get("path_provenance", ())),
        }
        group = groups.get(key)
        if group is None:
            groups[key] = {"record": record, "provenance": [provenance]}
            continue
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
        [group["record"] for group in grouped], aam, workers, static_context,
        include_event_relations)
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


def compile_mechanism_families(
        grouped: MechanismResult, *, workers: int = 1,
        minimum_events_only: bool = False) -> AnalyticalAAMResult:
    """Compile and maximally deduplicate exact cosets for AAM mechanisms.

    Geometry, chirality, and RMSD do not participate in this transformation.
    """
    if not isinstance(grouped, MechanismResult):
        raise TypeError("compile_mechanism_families requires group_mechanisms(aam)")
    aam = grouped.aam
    started = time.perf_counter()
    problem = aam.problem
    static_context = analytical_family_static_context(
        problem.reactant.elements, problem.reactant.wbo,
        problem.product.elements, problem.product.wbo,
        graph_floor=aam.config.graph_floor,
        dwbo_threshold=aam.config.event_threshold,
        metal_dwbo_threshold=aam.config.metal_event_threshold)
    sources = (
        grouped.minimum_event_mechanisms()
        if minimum_events_only else grouped.mechanisms)
    mechanisms = []
    for mechanism in sources:
        completed = [_branch_record(branch) for branch in mechanism.branches]
        maximal = _maximal_families(
            completed, aam, workers, static_context)
        branches = []
        for record, family in maximal:
            typed_branch = _branch_from_record(record)
            branches.append(AnalyticalBranch(typed_branch, family))
        mechanisms.append(AnalyticalMechanism(
            source=mechanism, branches=tuple(branches)))
    return AnalyticalAAMResult(
        aam=aam,
        mechanisms=tuple(mechanisms),
        elapsed_seconds=time.perf_counter() - started)


def compile_mapping_families(aam: AAMResult, *, workers: int = 1):
    """Compile complete raw branch relations without mechanism grouping."""
    if not isinstance(aam, AAMResult):
        raise TypeError('compile_mapping_families requires an AAMResult')
    started = time.perf_counter()
    problem, config = aam.problem, aam.config
    records = []
    for branch in aam.branches:
        mapping = branch.representative.as_dict()
        if len(mapping) != problem.atom_count:
            raise ValueError('full mapping-family compilation requires complete branches')
        records.append({'mapping': mapping, 'hierarchy': branch.hierarchy.to_record(),
            'cuts': sorted({cut for path in branch.paths for cut in path.context.cuts}),
            'encounter_count': len(branch.paths),
            'path_provenance': [{'terminal': path.terminal, 'transitions': path.transitions}
                                for path in branch.paths]})
    context = analytical_family_static_context(problem.reactant.elements, problem.reactant.wbo,
        problem.product.elements, problem.product.wbo, graph_floor=config.graph_floor,
        dwbo_threshold=config.event_threshold, metal_dwbo_threshold=config.metal_event_threshold)
    compiled = _maximal_families(records, aam, workers, context,
                                 include_event_relations=False)
    return MappingFamilyResult(aam,
        tuple(AnalyticalBranch(_branch_from_record(raw), family) for raw, family in compiled),
        time.perf_counter() - started)
