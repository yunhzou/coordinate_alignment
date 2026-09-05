"""Public typed assembly API using the common exact coverage decision graph."""
from collections import defaultdict
from fractions import Fraction
from itertools import product

from ..fragment_matching import materialize_target_coverage_orbit
from ..subgraph import _coerce_graph
from .decision_graph import CoverageDecisionGraph
from .models import RetroAssembly, RetroAssemblySearchResult
from .ranking import validate_atom_ownership


def assemble_fragment_cover(target, candidates, *, maximum_precursors=None,
        assembly_limit=None, require_attachment_bonds=False,
        allow_repeated_precursors=True, orbit_limit=None, iso_tolerance=0.5):
    """All full covers; optional explicit caller restrictions and output slicing.

    No cap limits traversal. Repeated copies at different target occupations
    remain distinct slots. The raw API uses source IDs as precursor identities;
    the catalog API additionally has canonical structure identities.
    """
    graph = _coerce_graph(target, 0.2)
    pools = defaultdict(dict)
    for candidate in candidates:
        for variant in materialize_target_coverage_orbit(candidate, graph,
                iso_tolerance=iso_tolerance, limit=orbit_limit):
            mask = sum(1 << atom for atom in variant.covered_target_atoms)
            mapping = dict(variant.mapping)
            shape = (mask, tuple(sorted(tuple(sorted(mapping[a] for a in part))
                                         for part in variant.retained_fragments)),
                     tuple(sorted(tuple(sorted((mapping[a], mapping[b])))
                                  for a, b in variant.preserved_source_bonds)),
                     variant.attachment_atoms_target)
            key = (variant.source_id, variant.mapping, variant.retained_fragments,
                   variant.boundary_bonds, variant.attachment_atoms_target)
            pools[shape].setdefault(key, variant)
    shapes = sorted(pools, key=lambda s: (-s[0].bit_count(), s))
    decisions = CoverageDecisionGraph.build((s[0] for s in shapes), len(graph),
                                             maximum_regions=maximum_precursors)
    assemblies = []
    for slots in decisions.paths():
        for selected in product(*(tuple(pools[shapes[slot]].values()) for slot in slots)):
            if not allow_repeated_precursors and len({c.source_id for c in selected}) != len(selected):
                continue
            records = []
            for candidate in selected:
                mapping = dict(candidate.mapping)
                carried = {tuple(sorted((mapping[a], mapping[b])))
                           for a, b in candidate.preserved_source_bonds}
                records.append({
                    "covered_target_atoms": candidate.covered_target_atoms,
                    "attachment_atoms_target": candidate.attachment_atoms_target,
                    "preserved_target_bonds": carried,
                })
            formed = validate_atom_ownership(records, graph.edges(), require_attachment_bonds)
            if formed is None:
                continue
            assemblies.append(RetroAssembly(tuple(selected), tuple(map(tuple, formed)),
                tuple(sorted((c.source_id, a, b) for c in selected for a, b in c.boundary_bonds))))
    def rank(assembly):
        total = sum(c.retained_size + sum(map(len, c.leftover_fragments)) for c in assembly.candidates)
        return (len(set(assembly.precursor_ids)), -Fraction(len(graph), total),
                tuple(sorted(assembly.precursor_ids)),
                tuple(c.mapping for c in assembly.candidates))
    assemblies.sort(key=rank)
    shown = assemblies if assembly_limit is None else assemblies[:assembly_limit]
    return RetroAssemblySearchResult(tuple(shown),
        "matched" if assemblies else "no_cover", True, assembly_limit)
