"""Public immutable domain objects for atom mapping and reaction analysis.

Computational stages exchange these objects.  JSON dictionaries belong only
at serialization/view boundaries and are deliberately absent from this
module's public contracts.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field
from types import MappingProxyType
from typing import Any, Mapping

import numpy as np

from .alignment.post_aam import AAMBranch, AtomBijection


def _readonly_array(value, *, shape=None, ndim=None):
    array = np.array(value, dtype=float, copy=True)
    if ndim is not None and array.ndim != int(ndim):
        raise ValueError(f"expected a {ndim}-dimensional array")
    if shape is not None and array.shape != tuple(shape):
        raise ValueError(f"expected array shape {tuple(shape)}, got {array.shape}")
    array.setflags(write=False)
    return array


def _frozen_mapping(value: Mapping[str, Any] | None):
    return MappingProxyType(dict(value or {}))


@dataclass(frozen=True)
class MolecularEndpoint:
    """One fully materialized molecular endpoint.

    Atom indices are local to the endpoint.  Coordinates and the complete WBO
    matrix are owned immutable copies, so later stages cannot silently mutate
    an AAM problem after search.
    """

    elements: tuple[str, ...]
    coordinates: np.ndarray
    wbo: np.ndarray
    label: str = ""
    energy: float | None = None
    metadata: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self):
        elements = tuple(str(element) for element in self.elements)
        if not elements:
            raise ValueError("a molecular endpoint must contain atoms")
        coordinates = _readonly_array(
            self.coordinates, shape=(len(elements), 3))
        wbo = _readonly_array(
            self.wbo, shape=(len(elements), len(elements)))
        if not np.allclose(wbo, wbo.T, atol=1e-12, rtol=0.0):
            raise ValueError("WBO matrix must be symmetric")
        if not np.all(np.isfinite(coordinates)) or not np.all(np.isfinite(wbo)):
            raise ValueError("endpoint arrays must contain finite values")
        object.__setattr__(self, "elements", elements)
        object.__setattr__(self, "coordinates", coordinates)
        object.__setattr__(self, "wbo", wbo)
        object.__setattr__(self, "label", str(self.label))
        object.__setattr__(
            self, "energy", None if self.energy is None else float(self.energy))
        object.__setattr__(self, "metadata", _frozen_mapping(self.metadata))

    @property
    def atom_count(self):
        return len(self.elements)


@dataclass(frozen=True)
class AAMProblem:
    """Reactant/product endpoint pair supplied to AAM search."""

    reactant: MolecularEndpoint
    product: MolecularEndpoint
    name: str = ""

    def __post_init__(self):
        if self.reactant.atom_count != self.product.atom_count:
            raise ValueError("AAM endpoints must have equal atom counts")
        if Counter(self.reactant.elements) != Counter(self.product.elements):
            raise ValueError("AAM endpoints must have equal compositions")
        object.__setattr__(self, "name", str(self.name))

    @property
    def atom_count(self):
        return self.reactant.atom_count


@dataclass(frozen=True)
class AAMSearchConfig:
    """All hypotheses that affect AAM search and mechanism classification."""

    cut_floor: float = 0.2
    graph_floor: float = 0.2
    iso_tolerance: float = 1.0
    event_threshold: float = 0.5
    metal_event_threshold: float | None = 0.3
    seed_count: int = 3
    branch_limit: int = 100
    task_chunksize: int = 1
    symmetry_repair: bool = True
    symmetry_repair_min_changes: int = 1
    symmetry_repair_max_evaluations: int = 20_000
    anchors: tuple[tuple[int, int], ...] = ()

    def __post_init__(self):
        if self.cut_floor <= 0 or self.graph_floor <= 0:
            raise ValueError("graph and cut floors must be positive")
        if self.iso_tolerance <= 0:
            raise ValueError("isomorphism tolerance must be positive")
        if self.event_threshold <= 0:
            raise ValueError("event threshold must be positive")
        if self.seed_count < 1 or self.branch_limit < 1:
            raise ValueError("seed count and branch limit must be positive")
        if self.task_chunksize < 1:
            raise ValueError("task chunksize must be positive")
        anchors = tuple(sorted(
            (int(source), int(target)) for source, target in self.anchors))
        if len({source for source, _target in anchors}) != len(anchors):
            raise ValueError("an AAM source atom cannot have multiple anchors")
        if len({target for _source, target in anchors}) != len(anchors):
            raise ValueError("an AAM target atom cannot have multiple anchors")
        object.__setattr__(self, "anchors", anchors)


@dataclass(frozen=True)
class AAMSearchMetrics:
    elapsed_seconds: float
    cut_count: int
    raw_result_count: int
    retained_branch_count: int
    max_live_branches: int
    max_growth_candidates: int
    subtree_branch_cap_count: int
    worker_returned_branch_count: int
    parent_merge_seconds: float

    @classmethod
    def from_record(cls, record, elapsed_seconds):
        record = dict(record or {})
        return cls(
            elapsed_seconds=float(elapsed_seconds),
            cut_count=int(record.get("cuts", 0)),
            raw_result_count=int(record.get("raw_result_count", 0)),
            retained_branch_count=int(record.get("retained_branch_count", 0)),
            max_live_branches=int(record.get("max_live_branches", 0)),
            max_growth_candidates=int(record.get("max_growth_candidates", 0)),
            subtree_branch_cap_count=int(
                record.get("subtree_branch_cap_count", 0)),
            worker_returned_branch_count=int(
                record.get("worker_returned_branch_count", 0)),
            parent_merge_seconds=float(record.get("parent_merge_seconds", 0.0)),
        )


@dataclass(frozen=True)
class AAMMechanism:
    """One exact event class and every completed AAM branch that realizes it."""

    key: tuple
    representative: AtomBijection
    branches: tuple[AAMBranch, ...]
    cuts: tuple[tuple[int, int], ...]
    includes_uncut_search: bool
    encounter_count: int

    def __post_init__(self):
        if not self.branches:
            raise ValueError("an AAM mechanism must retain at least one branch")
        if any(branch.representative.degree != self.representative.degree
               for branch in self.branches):
            raise ValueError("AAM mechanism branches disagree on atom count")
        object.__setattr__(self, "key", tuple(self.key))
        object.__setattr__(self, "branches", tuple(self.branches))
        object.__setattr__(self, "cuts", tuple(sorted(
            tuple(sorted(map(int, cut))) for cut in self.cuts)))
        object.__setattr__(self, "encounter_count", int(self.encounter_count))

    @property
    def event_count(self):
        return len(self.key[0]) + len(self.key[1])


@dataclass(frozen=True)
class AAMResult:
    """Complete output of AAM search, before geometry post-processing."""

    problem: AAMProblem
    config: AAMSearchConfig
    mechanisms: tuple[AAMMechanism, ...]
    metrics: AAMSearchMetrics

    def __post_init__(self):
        mechanisms = tuple(self.mechanisms)
        if any(mechanism.representative.degree != self.problem.atom_count
               for mechanism in mechanisms):
            raise ValueError("AAM result contains a wrong-degree mapping")
        if len({mechanism.key for mechanism in mechanisms}) != len(mechanisms):
            raise ValueError("AAM mechanism keys must be unique")
        object.__setattr__(self, "mechanisms", mechanisms)

    def minimum_event_mechanisms(self):
        if not self.mechanisms:
            return ()
        minimum = min(mechanism.event_count for mechanism in self.mechanisms)
        return tuple(
            mechanism for mechanism in self.mechanisms
            if mechanism.event_count == minimum)

