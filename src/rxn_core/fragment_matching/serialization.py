"""Strict serialization for fragment-detection records."""
from __future__ import annotations
import copy
from dataclasses import asdict, fields, replace

from ..alignment.post_aam import AAMHierarchy, AAMHierarchyChain
from ..search_graph import AAMSearchGraph, SearchPath
from .models import FragmentCandidate, FragmentDetectionResult, FragmentDerivation


FRAGMENT_DETECTION_SCHEMA = "rxn_core.fragment_detection/v7"


class _GraphArchive:
    def __init__(self, graphs=()):
        self.graphs = []
        self.ids = {}
        self.fragments = []
        self.fragment_ids = {}
        self.hierarchies = {}
        self.generators = []
        self.generator_ids = {}
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
        # Archive encoding shares groups; the ordinary AAM graph API is unchanged.
        return [{"schema": "rxn_core.aam_search_graph_refs/v1",
                 "contexts": [asdict(c) for c in graph.contexts], "roots": graph.roots,
                 "states": [asdict(s) for s in graph.states],
                 "transitions": [{**{f.name: getattr(edge, f.name) for f in fields(edge)
                                      if f.name != "match"},
                                  "match": None if edge.match is None else self.encode_fragment(edge.match)}
                                 for edge in graph.transitions],
                 "stops": [asdict(s) for s in graph.stops]} for graph in self.graphs]

    def encode_fragment(self, fragment):
        symmetry = fragment["symmetry"]
        encoded = copy.deepcopy({k: v for k, v in symmetry.items() if k != "automorph_generators"})
        if "automorph_generators" in symmetry:
            ids = []
            for raw in symmetry["automorph_generators"]:
                images = tuple(raw)
                if images not in self.generator_ids:
                    self.generator_ids[images] = len(self.generators)
                    self.generators.append(images)
                ids.append(self.generator_ids[images])
            encoded["automorph_generator_ids"] = ids
        return {**copy.deepcopy({k: v for k, v in fragment.items() if k != "symmetry"}),
                "symmetry": encoded}

    def hierarchy_reference(self, hierarchy):
        return {"segments": [self.segment_reference(base, action)
                             for base, action in hierarchy.segments]}

    def segment_reference(self, base, action):
        # Hold the immutable base, not just its id, for this archive's lifetime.
        key = id(base)
        if key not in self.hierarchies:
            ids = []
            for fragment in base.fragments:
                if fragment not in self.fragment_ids:
                    self.fragment_ids[fragment] = len(self.fragments)
                    self.fragments.append(fragment)
                ids.append(self.fragment_ids[fragment])
            self.hierarchies[key] = base, ids
        return {"fragment_ids": self.hierarchies[key][1],
                "target_action": action}

    def fragment_records(self):
        return [self.encode_fragment(f) for f in
                AAMHierarchy(tuple(self.fragments)).to_record()["fragments"]]


def fragment_archive_from_record(record):
    """Resolve a shared archive once; all candidates reuse its immutable objects."""
    from ..alignment.post_aam import AtomPermutation
    generators = tuple(AtomPermutation(tuple(g)) for g in record["generators"])
    def decode_fragment(raw):
        symmetry = raw["symmetry"]
        decoded = {k: v for k, v in symmetry.items() if k != "automorph_generator_ids"}
        if "automorph_generator_ids" in symmetry:
            decoded["automorph_generators"] = tuple(generators[i].images
                for i in symmetry["automorph_generator_ids"])
        return {**raw, "symmetry": decoded}
    graphs = []
    for graph in record.get("search_graphs", ()):
        if graph["schema"] != "rxn_core.aam_search_graph_refs/v1":
            raise ValueError("unsupported fragment search archive")
        decoded = {**graph, "schema": "rxn_core.aam_search_graph/v1",
                   "transitions": [{**edge, "match": None if edge["match"] is None else
                                    decode_fragment(edge["match"])} for edge in graph["transitions"]]}
        graphs.append(AAMSearchGraph.from_record(decoded))
    fragments = []
    for raw in record["hierarchy_fragments"]:
        symmetry = raw["symmetry"]
        # Reuse validated permutations instead of reconstructing them per fragment.
        fragment = AAMHierarchy.from_record({"fragments": ({**raw, "symmetry": {
            k: v for k, v in symmetry.items() if k != "automorph_generator_ids"}},)}).fragments[0]
        if "automorph_generator_ids" in symmetry:
            fragment = replace(fragment, target_generators=tuple(generators[i]
                for i in symmetry["automorph_generator_ids"]))
        fragments.append(fragment)
    return tuple(graphs), tuple(fragments)


def repack_fragment_detection_v4(record):
    """Explicit, lossless archive migration. Does not perform any matching."""
    if record["schema"] != "rxn_core.fragment_detection/v4":
        raise ValueError("repacking requires a v4 augmented-AAM record")
    archive = _GraphArchive(tuple(AAMSearchGraph.from_record(g) for g in record["search_graphs"]))
    candidates = [dict(c, aam_hierarchy=archive.hierarchy_reference(
        AAMHierarchy.from_record(c["aam_hierarchy"]))) for c in record["candidates"]]
    return {**record, "schema": FRAGMENT_DETECTION_SCHEMA, "candidates": candidates,
            "search_graphs": archive.records(), "hierarchy_fragments": archive.fragment_records(),
            "generators": archive.generators}


def repack_fragment_detection_v6(record):
    """Explicit lossless migration to segmented hierarchy references; no AAM."""
    if record["schema"] != "rxn_core.fragment_detection/v6":
        raise ValueError("repacking requires a v6 augmented-AAM record")
    return {**record, "schema": FRAGMENT_DETECTION_SCHEMA,
            "candidates": [dict(c, aam_hierarchy={"segments": [c["aam_hierarchy"]]})
                           for c in record["candidates"]]}


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
        "fragment_classes": list(candidate.fragment_classes),
        "preserved_source_bonds": [list(edge) for edge in candidate.preserved_source_bonds],
        "aam_hierarchy": archive.hierarchy_reference(candidate.aam_hierarchy),
        "derivations": [{
            "initial_paths": [archive.reference(p) for p in d.initial_paths],
            "residual_paths": [archive.reference(p) for p in d.residual_paths],
            "target_action": [list(pair) for pair in d.target_action],
            "occupation_projected": d.occupation_projected,
        } for d in candidate.derivations],
    }
    if standalone:
        record["search_graphs"] = archive.records()
        record["hierarchy_fragments"] = archive.fragment_records()
        record["generators"] = archive.generators
    return record


def fragment_candidate_from_record(record, *, search_graphs=None, hierarchy_fragments=None):
    if search_graphs is None:
        graphs, fragments = fragment_archive_from_record(record)
    else:
        graphs, fragments = search_graphs, hierarchy_fragments
    reference = record["aam_hierarchy"]
    parts = tuple(AAMHierarchy(tuple(fragments[i] for i in segment["fragment_ids"])).relabel_target(
        segment["target_action"]) for segment in reference["segments"])
    hierarchy = parts[0] if len(parts) == 1 else AAMHierarchyChain(parts)
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
        aam_hierarchy=hierarchy,
        fragment_classes=tuple(record.get("fragment_classes", ())),
        preserved_source_bonds=tuple(tuple(edge) for edge in record.get("preserved_source_bonds", ())),
        derivations=tuple(FragmentDerivation(
            tuple(SearchPath.from_reference(p, graphs) for p in d["initial_paths"]),
            tuple(SearchPath.from_reference(p, graphs) for p in d["residual_paths"]),
            tuple(tuple(pair) for pair in d["target_action"]),
            bool(d.get("occupation_projected", False)))
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
        "hierarchy_fragments": archive.fragment_records(),
        "generators": archive.generators,
    }
