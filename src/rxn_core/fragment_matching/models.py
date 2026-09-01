"""Immutable records for augmented fragment detection."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any


@dataclass(frozen=True)
class FragmentDetectionConfig:
    graph_floor: float = 0.2
    iso_tolerance: float = 0.5
    minimum_fragment_size: int = 1
    branch_limit: int = 100
    candidate_limit: int = 512
    seed_limit: int | None = None
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
