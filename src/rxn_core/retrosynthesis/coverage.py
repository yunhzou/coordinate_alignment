"""Target-coverage assembly from precursor fragment candidates."""
from __future__ import annotations

from ..fragment_matching import materialize_target_coverage_orbit
from ..matcher import _nauty_atom_generators
from ..subgraph import _coerce_graph
from .models import (
    RetroAssembly,
    RetroAssemblySearchResult,
)


def assemble_fragment_cover(
        target, candidates, *, maximum_precursors=2,
        assembly_limit=1_000, require_attachment_bonds=False,
        allow_repeated_precursors=True, orbit_limit=100_000,
        iso_tolerance=0.5):
    """Combine candidates into non-overlapping complete target covers."""
    if maximum_precursors < 1 or assembly_limit < 1:
        raise ValueError("precursor and assembly limits must be positive")
    graph_P = _coerce_graph(target, 0.2)
    target_atoms = frozenset(map(int, graph_P.nodes()))
    generators = _nauty_atom_generators(
        graph_P, wbo_tol=iso_tolerance)
    expanded = []
    for candidate in candidates:
        variants = materialize_target_coverage_orbit(
            candidate,
            graph_P,
            iso_tolerance=iso_tolerance,
            limit=orbit_limit,
            generators=generators,
        )
        expanded.extend(variants)
    ordered = sorted(expanded, key=lambda candidate: (
        -candidate.retained_size,
        candidate.source_id,
        candidate.covered_target_atoms,
        candidate.mapping,
    ))
    assemblies = []
    capped = False

    def emit(selected):
        owner = {}
        for index, candidate in enumerate(selected):
            for atom in candidate.covered_target_atoms:
                owner[atom] = index
        formed = []
        for left, right in graph_P.edges():
            if owner[left] == owner[right]:
                continue
            if require_attachment_bonds:
                left_candidate = selected[owner[left]]
                right_candidate = selected[owner[right]]
                if (left not in left_candidate.attachment_atoms_target
                        or right not in right_candidate.attachment_atoms_target):
                    return
            formed.append(tuple(sorted((int(left), int(right)))))
        broken = tuple(sorted(
            (candidate.source_id, int(left), int(right))
            for candidate in selected
            for left, right in candidate.boundary_bonds
        ))
        assemblies.append(RetroAssembly(
            candidates=tuple(selected),
            formed_bonds=tuple(sorted(formed)),
            broken_bonds=broken,
        ))

    def visit(start, selected, covered, used_precursors):
        nonlocal capped
        if capped:
            return
        if covered == target_atoms:
            emit(selected)
            if len(assemblies) >= assembly_limit:
                capped = True
            return
        if len(selected) >= maximum_precursors:
            return
        for index in range(start, len(ordered)):
            candidate = ordered[index]
            coverage = frozenset(candidate.covered_target_atoms)
            if coverage & covered:
                continue
            if (not allow_repeated_precursors
                    and candidate.source_id in used_precursors):
                continue
            visit(
                index + 1,
                selected + [candidate],
                covered | coverage,
                used_precursors | {candidate.source_id},
            )

    visit(0, [], frozenset(), set())
    assemblies.sort(key=lambda assembly: (
        len(assembly.candidates),
        len(assembly.formed_bonds),
        len(assembly.broken_bonds),
        assembly.precursor_ids,
    ))
    return RetroAssemblySearchResult(
        assemblies=tuple(assemblies),
        status=("capped" if capped
                else ("matched" if assemblies else "no_cover")),
        complete=not capped,
        assembly_limit=int(assembly_limit),
    )
