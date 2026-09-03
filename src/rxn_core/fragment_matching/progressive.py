"""Progressive multi-precursor matching against uncovered target atoms."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..frag import WeightedGraph
from .detection import detect_fragments
from .models import FragmentDetectionConfig


@dataclass(frozen=True)
class ProgressiveFragmentPlacement:
    source_id: str
    mapping: tuple[tuple[int, int], ...]
    retained_fragments: tuple[tuple[int, ...], ...]


@dataclass(frozen=True)
class ProgressiveFragmentMatchingResult:
    placements: tuple[ProgressiveFragmentPlacement, ...]
    uncovered_target_atoms: tuple[int, ...]


def _induced(graph, atoms):
    atoms = tuple(sorted(atoms))
    network = graph.to_networkx()
    return (
        WeightedGraph(
            [dict(network.nodes[atom]) for atom in atoms],
            np.asarray(graph.weights)[np.ix_(atoms, atoms)],
        ),
        atoms,
    )


def progressive_fragment_matching(
        sources, target, *, config: FragmentDetectionConfig | None = None):
    """Repeatedly apply AAM to residual R and currently uncovered P.

    Every candidate returned by fragment detection participates in the
    deterministic marginal-coverage choice. Source and target atoms are
    consumed only by an actual AAM mapping, so the combined result never
    fabricates or reuses an assignment.
    """
    config = config or FragmentDetectionConfig(seed_mode="all")
    sources = tuple((str(source_id), graph) for source_id, graph in sources)
    remaining_sources = [set(range(len(graph.nodes))) for _id, graph in sources]
    remaining_target = set(range(len(target.nodes)))
    mappings = [[] for _source in sources]
    fragments = [[] for _source in sources]

    while remaining_target:
        options = []
        for source_index, (source_id, source) in enumerate(sources):
            if not remaining_sources[source_index]:
                continue
            source_subgraph, source_atoms = _induced(
                source, remaining_sources[source_index])
            target_subgraph, target_atoms = _induced(target, remaining_target)
            detection = detect_fragments(
                source_subgraph,
                target_subgraph,
                source_id=source_id,
                config=config,
            )
            for candidate in detection.candidates:
                mapping = tuple(sorted(
                    (source_atoms[source_atom], target_atoms[target_atom])
                    for source_atom, target_atom in candidate.mapping
                ))
                if not mapping:
                    continue
                retained_fragments = tuple(
                    tuple(sorted(source_atoms[atom] for atom in fragment))
                    for fragment in candidate.retained_fragments
                    if fragment
                ) or (tuple(source_atom for source_atom, _ in mapping),)
                heavy_atoms = sum(
                    source.nodes[source_atom]["element"] != "H"
                    for source_atom, _target_atom in mapping
                )
                options.append((
                    (-len(mapping), -heavy_atoms,
                     len(candidate.boundary_bonds), source_id, mapping),
                    source_index,
                    mapping,
                    retained_fragments,
                ))
        if not options:
            break
        _rank, source_index, mapping, retained_fragments = min(options)
        mappings[source_index].extend(mapping)
        fragments[source_index].extend(retained_fragments)
        remaining_sources[source_index].difference_update(
            source_atom for source_atom, _target_atom in mapping)
        remaining_target.difference_update(
            target_atom for _source_atom, target_atom in mapping)

    return ProgressiveFragmentMatchingResult(
        placements=tuple(
            ProgressiveFragmentPlacement(
                source_id=source_id,
                mapping=tuple(sorted(mapping)),
                retained_fragments=tuple(fragment),
            )
            for (source_id, _graph), mapping, fragment in zip(
                sources, mappings, fragments, strict=True)
        ),
        uncovered_target_atoms=tuple(sorted(remaining_target)),
    )
