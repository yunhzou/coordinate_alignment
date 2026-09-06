"""Optional balanced-reaction event grouping of recorded AAM paths."""
from __future__ import annotations

import time

from .alignment.post_aam import AAMBranch, AAMHierarchy, AtomBijection, PermutationGroup
from .alignment.sweep import (_branch_symmetry_record, _cut_sweep_cfg,
                              _score_branch_mapping, _pool_add,
                              _MechanismEventCanonicalizer)
from .domain import AAMMechanism, MechanismResult
from .frag import build_graph, expand_mapping
from .matcher import _nauty_orbits


def _branch_from_record(raw):
    mapping_family = dict(raw.get('mapping_family') or {})
    generators = raw.get('target_group_generators', mapping_family.get('target_generators'))
    return AAMBranch(
        representative=AtomBijection.from_mapping(raw['mapping']),
        hierarchy=AAMHierarchy.from_record(raw['hierarchy']),
        encounter_count=int(raw.get('encounter_count', 1)),
        cuts=tuple(tuple(cut) for cut in raw.get('cuts', ())),
        covered_path_count=int(raw.get('covered_path_count', 1)),
        mapping_family=mapping_family,
        path_provenance=tuple(raw.get('path_provenance', ())),
        target_group=None if generators is None else PermutationGroup.from_generator_mappings(
            len(raw['mapping']), generators))


def group_mechanisms(aam):
    """Reproduce balanced event scoring without mutating raw search evidence."""
    started = time.perf_counter()
    problem, config = aam.problem, aam.config
    if not problem.balanced:
        raise ValueError('Balanced mechanism grouping requires balanced endpoints; '
                         'partial AAM mappings remain available in aam.graph')
    source = build_graph(problem.reactant.elements, problem.reactant.wbo,
                         bond_cut=config.graph_floor)
    target = build_graph(problem.product.elements, problem.product.wbo,
                         bond_cut=config.graph_floor)
    ro = _nauty_orbits(source, wbo_tol=config.iso_tolerance)
    po = _nauty_orbits(target, wbo_tol=config.iso_tolerance)
    canonicalizer = _MechanismEventCanonicalizer(source, wbo_tol=config.iso_tolerance)
    cfg = _cut_sweep_cfg(graph_floor=config.graph_floor, iso_tol=config.iso_tolerance,
        dwbo_threshold=config.event_threshold, metal_dwbo_threshold=config.metal_event_threshold,
        symmetry_repair=config.symmetry_repair,
        symmetry_repair_min_changes=config.symmetry_repair_min_changes,
        symmetry_repair_max_evals=config.symmetry_repair_max_evaluations,
        n_atoms=problem.atom_count, anchor_map=dict(config.anchors))
    pool, scored_mappings = {}, {}
    for path in aam.graph.paths():
        mapping = expand_mapping(path.mapping, source, target)
        mapping_key = tuple(sorted(mapping.items()))
        if mapping_key not in scored_mappings:
            scored_mappings[mapping_key] = _score_branch_mapping(mapping, source, target,
                problem.reactant.wbo, problem.product.wbo, source, po, ro, (), cfg,
                problem.reactant.elements, problem.product.elements,
                event_canonicalizer=canonicalizer)
        scored = scored_mappings[mapping_key]
        if scored is None:
            continue
        signature, selected = scored
        hierarchy = _branch_symmetry_record(path, symmetry_fragments=path.fragments)
        _pool_add(pool, signature, selected, path.context.cuts, hierarchy,
                  provenance={'terminal': path.terminal, 'transitions': path.transitions,
                    'selection_action': tuple(sorted((image, selected[atom])
                        for atom, image in mapping.items() if image != selected[atom]))})
    mechanisms = tuple(AAMMechanism(key=tuple(key),
        representative=AtomBijection.from_mapping(entry['mapping']),
        branches=tuple(_branch_from_record(raw) for raw in entry['branches']),
        cuts=tuple(entry['cuts']), includes_uncut_search=entry['has_no_cut'],
        encounter_count=entry['dedup_count']) for key, entry in pool.items())
    return MechanismResult(aam, mechanisms, time.perf_counter()-started)
