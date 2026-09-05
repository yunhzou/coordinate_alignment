"""Strict serialization for fragment-detection records."""
from __future__ import annotations

from ..alignment.post_aam import AAMHierarchy
from ..search_graph import AAMSearchGraph, SearchPath
from .models import FragmentCandidate, FragmentDetectionResult, FragmentDerivation


FRAGMENT_DETECTION_SCHEMA = "rxn_core.fragment_detection/v3"


class _GraphArchive:
    def __init__(self, graphs=()):
        self.graphs = []
        self.ids = {}
        for graph in graphs:
            self.add(graph)

    def add(self, graph):
        if id(graph) not in self.ids:
            self.ids[id(graph)] = len(self.graphs)
            self.graphs.append(graph)
        return self.ids[id(graph)]

    def reference(self, path):
        return path.to_reference(self.add(path.graph))

    def records(self):
        return [graph.to_record() for graph in self.graphs]


def fragment_candidate_to_record(candidate: FragmentCandidate, *, archive=None):
    standalone = archive is None
    archive = _GraphArchive() if standalone else archive
    record = {
        "mapping": [list(item) for item in candidate.mapping],
        "retained_atoms": list(candidate.retained_atoms),
        "covered_target_atoms": list(candidate.covered_target_atoms),
        "leftover_fragments": [
            list(item) for item in candidate.leftover_fragments
        ],
        "boundary_bonds": [list(item) for item in candidate.boundary_bonds],
        "attachment_atoms_source": list(candidate.attachment_atoms_source),
        "attachment_atoms_target": list(candidate.attachment_atoms_target),
        "copied_residual_placements": [
            list(item) for item in candidate.copied_residual_placements
        ],
        "augmented_target_atom_count": candidate.augmented_target_atom_count,
        "retained_fragments": [
            list(item) for item in candidate.retained_fragments
        ],
        "aam_hierarchy": candidate.aam_hierarchy.to_record(),
        "derivations": [{
            "initial_paths": [archive.reference(p) for p in d.initial_paths],
            "residual_paths": [archive.reference(p) for p in d.residual_paths],
            "target_action": [list(pair) for pair in d.target_action],
        } for d in candidate.derivations],
    }
    if standalone:
        record["search_graphs"] = archive.records()
    return record


def fragment_candidate_from_record(record, *, search_graphs=None):
    graphs = (tuple(AAMSearchGraph.from_record(g) for g in record.get("search_graphs", ()))
              if search_graphs is None else search_graphs)
    return FragmentCandidate(
        source_id=str(record.get("source_id", "")),
        mapping=tuple(tuple(map(int, item))
                      for item in record.get("mapping") or ()),
        retained_atoms=tuple(map(int, record.get("retained_atoms") or ())),
        covered_target_atoms=tuple(map(
            int, record.get("covered_target_atoms") or ())),
        leftover_fragments=tuple(
            tuple(map(int, item))
            for item in record.get("leftover_fragments") or ()),
        boundary_bonds=tuple(
            tuple(map(int, item))
            for item in record.get("boundary_bonds") or ()),
        attachment_atoms_source=tuple(map(
            int, record.get("attachment_atoms_source") or ())),
        attachment_atoms_target=tuple(map(
            int, record.get("attachment_atoms_target") or ())),
        copied_residual_placements=tuple(
            tuple(map(int, item))
            for item in record.get("copied_residual_placements") or ()),
        augmented_target_atom_count=int(
            record.get("augmented_target_atom_count", 0)),
        retained_fragments=tuple(
            tuple(map(int, item))
            for item in record.get("retained_fragments") or ()),
        aam_hierarchy=AAMHierarchy.from_record(
            record.get("aam_hierarchy") or {}),
        derivations=tuple(FragmentDerivation(
            tuple(SearchPath.from_reference(p, graphs) for p in d["initial_paths"]),
            tuple(SearchPath.from_reference(p, graphs) for p in d["residual_paths"]),
            tuple(tuple(pair) for pair in d["target_action"]))
            for d in record.get("derivations", ())),
    )


def fragment_detection_to_record(
        result: FragmentDetectionResult, *, row_index, representation,
        candidates=None):
    selected = result.candidates if candidates is None else tuple(candidates)
    archive = _GraphArchive(result.search_graphs)
    candidates = [fragment_candidate_to_record(candidate, archive=archive)
                  for candidate in selected]
    return {
        "schema": FRAGMENT_DETECTION_SCHEMA,
        "row_index": int(row_index),
        "source_id": result.source_id,
        "representation": representation,
        "status": result.status,
        "complete": result.complete,
        "branch_limit": result.branch_limit,
        "maximum_branch_count": result.maximum_branch_count,
        "capped_seed_count": result.capped_seed_count,
        "best_fragment_size": result.best_fragment_size,
        "initial_placement_encounters": result.initial_placement_encounters,
        "initial_family_count": result.initial_family_count,
        "best_initial_family_count": result.best_initial_family_count,
        "seed_attempt_count": result.seed_attempt_count,
        "seed_pruned_count": result.seed_pruned_count,
        "rough_stop_hit": result.rough_stop_hit,
        "candidates": candidates,
        "search_graphs": archive.records(),
    }
