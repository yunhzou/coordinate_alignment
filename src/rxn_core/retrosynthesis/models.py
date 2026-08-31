"""Immutable records for precursor assembly."""
from __future__ import annotations

from dataclasses import dataclass

from ..fragment_matching.models import FragmentCandidate


@dataclass(frozen=True)
class RetroAssembly:
    candidates: tuple[FragmentCandidate, ...]
    formed_bonds: tuple[tuple[int, int], ...]
    broken_bonds: tuple[tuple[str, int, int], ...]

    @property
    def precursor_ids(self):
        return tuple(candidate.source_id for candidate in self.candidates)


@dataclass(frozen=True)
class RetroAssemblySearchResult:
    assemblies: tuple[RetroAssembly, ...]
    status: str
    complete: bool
    assembly_limit: int
