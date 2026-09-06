"""Connected-fragment discovery without residual augmentation.

This is the same initial search used by full detection, exposed independently.
Search paths and compressed symmetry remain attached to every candidate.
"""
from dataclasses import dataclass

from .detection import _initial_fragment_placements, _prepare_fragment_detection
from .graph_ops import partition_at_retained_fragment
from .models import FragmentCandidate, FragmentDerivation, FragmentDetectionConfig


@dataclass(frozen=True)
class ConnectedFragmentResult:
    source_id: str
    candidates: tuple[FragmentCandidate, ...]
    search_graphs: tuple
    capped_seed_count: int
    maximum_branch_count: int
    seed_attempt_count: int
    complete: bool


def find_connected_fragments(source, target, *, source_id="", config=None):
    """Return discovered connected families, largest first; do not augment.

    'Largest' means largest discovered by the configured AAM search, not a
    certificate of a globally maximum common subgraph.
    """
    config = config or FragmentDetectionConfig()
    source, context, _ = _prepare_fragment_detection(source, target, config)
    (families, capped, branches, candidate_capped, seed_limited, attempts,
     pruned, rough_stop, graphs) = _initial_fragment_placements(
        source, context.graph, config, target_orbits=context.atom_orbits)
    candidates = []
    for family in families:
        retained = family.retained_atoms
        mapping = dict(family.representative_mapping)
        _, boundary, leftovers = partition_at_retained_fragment(source, retained)
        boundary = tuple(sorted(set(boundary) | set(family.search_paths[0].deferred_edges)))
        attachments = tuple(sorted({a for edge in boundary for a in edge if a in mapping}))
        candidates.append(FragmentCandidate(
            source_id=str(source_id), mapping=family.representative_mapping,
            retained_atoms=retained, covered_target_atoms=tuple(sorted(mapping.values())),
            leftover_fragments=leftovers, boundary_bonds=boundary,
            attachment_atoms_source=attachments,
            attachment_atoms_target=tuple(sorted(mapping[a] for a in attachments)),
            copied_residual_placements=(), augmented_target_atom_count=len(context.graph),
            retained_fragments=(retained,), fragment_classes=(0,),
            preserved_source_bonds=tuple(sorted(tuple(sorted((a, b)))
                for a, b in source.edges() if a in mapping and b in mapping
                and tuple(sorted((a, b))) not in boundary)),
            aam_hierarchy=family.search_paths[0].hierarchy,
            derivations=(FragmentDerivation(family.search_paths),),
        ))
    candidates.sort(key=lambda c: (-c.retained_size, c.covered_target_atoms, c.mapping))
    approximate = config.seed_mode != "all" and bool(pruned or rough_stop)
    return ConnectedFragmentResult(str(source_id), tuple(candidates), graphs,
        capped + int(candidate_capped), branches, attempts,
        not (capped or candidate_capped or seed_limited or approximate))
