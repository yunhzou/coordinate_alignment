"""Immutable records for augmented fragment detection."""
from __future__ import annotations

from dataclasses import dataclass
from dataclasses import field
from typing import Any

from ..alignment.post_aam import AAMHierarchy
from ..search_graph import AAMSearchGraph, SearchPath


@dataclass(frozen=True)
class FragmentDerivation:
    """Recorded searches behind an occupation, in their original atom frame.

    initial_paths[0] is the selected initial witness; remaining initial paths
    are equivalent discoveries, not independent fragments to concatenate.
    target_action transports the selected result without rewriting evidence.
    """
    initial_paths: tuple[SearchPath, ...]
    residual_paths: tuple[SearchPath, ...] = ()
    target_action: tuple[tuple[int, int], ...] = ()


@dataclass(frozen=True)
class FragmentDetectionConfig:
    graph_floor: float = 0.2
    iso_tolerance: float = 0.5
    minimum_fragment_size: int = 1
    branch_limit: int = 100
    candidate_limit: int = 512
    seed_limit: int | None = None
    seed_mode: str = "all"
    rough_retention_threshold: float = 0.5
    maximum_boundary_bonds: int | None = None
    maximum_leftover_fragments: int | None = None

    def __post_init__(self):
        if self.graph_floor <= 0 or self.iso_tolerance <= 0:
            raise ValueError(
                "graph floor and isomorphism tolerance must be positive")
        if self.minimum_fragment_size < 1:
            raise ValueError("minimum fragment size must be positive")
        if self.branch_limit < 1 or self.candidate_limit < 1:
            raise ValueError("branch and candidate limits must be positive")
        if self.seed_limit is not None and self.seed_limit < 1:
            raise ValueError("seed limit must be positive")
        if self.seed_mode not in {
                "all", "fragment_cover", "orbit_representatives"}:
            raise ValueError(
                "seed mode must be 'all', 'fragment_cover', or "
                "'orbit_representatives'")
        if not 0 < self.rough_retention_threshold <= 1:
            raise ValueError("rough retention threshold must be in (0, 1]")


@dataclass(frozen=True)
class FragmentTargetContext:
    """Reusable target graph and symmetry data for bank-scale detection."""

    graph: Any
    atom_orbits: Any
    automorphism_generators: tuple
    graph_floor: float
    iso_tolerance: float


@dataclass(frozen=True)
class FragmentCandidate:
    """Target-owned fragments detected from one source graph."""

    source_id: str
    mapping: tuple[tuple[int, int], ...]
    retained_atoms: tuple[int, ...]
    covered_target_atoms: tuple[int, ...]
    leftover_fragments: tuple[tuple[int, ...], ...]
    boundary_bonds: tuple[tuple[int, int], ...]
    attachment_atoms_source: tuple[int, ...]
    attachment_atoms_target: tuple[int, ...]
    copied_residual_placements: tuple[tuple[int, int], ...]
    augmented_target_atom_count: int
    retained_fragments: tuple[tuple[int, ...], ...] = ()
    aam_hierarchy: AAMHierarchy = field(
        default_factory=lambda: AAMHierarchy(()))
    derivations: tuple[FragmentDerivation, ...] = ()

    @property
    def atom_mapping(self):
        return dict(self.mapping)

    @property
    def retained_size(self):
        return len(self.retained_atoms)


@dataclass(frozen=True)
class FragmentDetectionResult:
    source_id: str
    candidates: tuple[FragmentCandidate, ...]
    status: str
    complete: bool
    branch_limit: int
    maximum_branch_count: int
    capped_seed_count: int
    best_fragment_size: int
    initial_placement_encounters: int
    initial_family_count: int
    best_initial_family_count: int
    seed_attempt_count: int
    seed_pruned_count: int
    rough_stop_hit: bool
    search_graphs: tuple[AAMSearchGraph, ...] = ()
