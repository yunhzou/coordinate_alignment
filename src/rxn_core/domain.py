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

from .alignment.post_aam import AAMBranch, AAMHierarchy, AtomBijection
from .search_graph import AAMSearchGraph


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
        object.__setattr__(self, "name", str(self.name))

    @property
    def balanced(self):
        return Counter(self.reactant.elements) == Counter(self.product.elements)

    @property
    def source_atom_count(self):
        return self.reactant.atom_count

    @property
    def target_atom_count(self):
        return self.product.atom_count

    @property
    def atom_count(self):
        """Source endpoint size; use target_atom_count for the other side."""
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
    completed_group_requests: int = 0
    completed_group_calculations: int = 0
    completed_group_cache_hits: int = 0
    worker_search_seconds: float | None = None
    checkpoint_seconds: float | None = None
    symmetry_finalization_seconds: float | None = None
    checkpoint_restore_and_finalize_seconds: float | None = None

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
            completed_group_requests=int(record.get(
                "completed_candidate_group_requests", 0)),
            completed_group_calculations=int(record.get(
                "completed_candidate_group_calculations", 0)),
            completed_group_cache_hits=int(record.get(
                "completed_candidate_group_cache_hits", 0)),
            worker_search_seconds=record.get('worker_search_seconds'),
            checkpoint_seconds=record.get('checkpoint_seconds'),
            symmetry_finalization_seconds=record.get('symmetry_finalization_seconds'),
            checkpoint_restore_and_finalize_seconds=record.get('checkpoint_restore_and_finalize_seconds'),
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
    """Raw fragment-decision graph, independent of mechanism classification."""

    problem: AAMProblem
    config: AAMSearchConfig
    graph: AAMSearchGraph
    metrics: AAMSearchMetrics

    @property
    def branches(self):
        return self.graph.branches()


@dataclass(frozen=True)
class MechanismResult:
    """Optional event grouping with a reference to the unmodified search."""

    aam: AAMResult
    mechanisms: tuple[AAMMechanism, ...]
    elapsed_seconds: float

    def minimum_event_mechanisms(self):
        if not self.mechanisms:
            return ()
        minimum = min(mechanism.event_count for mechanism in self.mechanisms)
        return tuple(
            mechanism for mechanism in self.mechanisms
            if mechanism.event_count == minimum)


@dataclass(frozen=True)
class MappingFamilyResult:
    """Compiled branch relations without mechanism grouping."""

    aam: AAMResult
    branches: tuple[AnalyticalBranch, ...]
    elapsed_seconds: float


@dataclass(frozen=True)
class AnalyticalBranch:
    """One maximal exact mapping coset compiled from completed AAM paths."""

    aam_branch: AAMBranch
    family: Any

    def __post_init__(self):
        if self.family is None:
            raise ValueError("an analytical branch requires an exact family")
        if int(self.family.degree) != self.aam_branch.representative.degree:
            raise ValueError("analytical family degree differs from AAM branch")

    @property
    def representative(self):
        return AtomBijection.from_mapping(self.family.representative_mapping)


@dataclass(frozen=True)
class AnalyticalMechanism:
    source: AAMMechanism
    branches: tuple[AnalyticalBranch, ...]

    def __post_init__(self):
        branches = tuple(self.branches)
        if not branches:
            raise ValueError("an analytical mechanism must contain families")
        object.__setattr__(self, "branches", branches)

    @property
    def key(self):
        return self.source.key


@dataclass(frozen=True)
class AnalyticalAAMResult:
    """AAM mechanisms after exact maximal-coset compilation."""

    aam: AAMResult
    mechanisms: tuple[AnalyticalMechanism, ...]
    elapsed_seconds: float

    def __post_init__(self):
        mechanisms = tuple(self.mechanisms)
        if len({mechanism.key for mechanism in mechanisms}) != len(mechanisms):
            raise ValueError("analytical mechanism keys must be unique")
        object.__setattr__(self, "mechanisms", mechanisms)
        object.__setattr__(self, "elapsed_seconds", float(self.elapsed_seconds))


@dataclass(frozen=True)
class RPMechanism:
    """One geometry-selected mapping for an analytical mechanism family."""

    analytical: AnalyticalMechanism
    mapping: AtomBijection
    broken_bonds: tuple[tuple[int, int], ...]
    formed_bonds: tuple[tuple[int, int], ...]
    core_atoms: tuple[int, ...]
    fixed_mapping_rmsd: float
    chirality: Mapping[str, Any]
    selected_branch_index: int

    def __post_init__(self):
        if self.mapping.degree != self.analytical.source.representative.degree:
            raise ValueError("selected R/P mapping has the wrong degree")
        object.__setattr__(self, "broken_bonds", tuple(sorted(
            tuple(sorted(map(int, bond))) for bond in self.broken_bonds)))
        object.__setattr__(self, "formed_bonds", tuple(sorted(
            tuple(sorted(map(int, bond))) for bond in self.formed_bonds)))
        object.__setattr__(self, "core_atoms", tuple(sorted(
            map(int, self.core_atoms))))
        object.__setattr__(
            self, "fixed_mapping_rmsd", float(self.fixed_mapping_rmsd))
        object.__setattr__(self, "chirality", _frozen_mapping(self.chirality))
        object.__setattr__(
            self, "selected_branch_index", int(self.selected_branch_index))


@dataclass(frozen=True)
class RPResult:
    """Final R/P mappings after analytical chirality and RMSD processing."""

    analytical: AnalyticalAAMResult
    mechanisms: tuple[RPMechanism, ...]
    elapsed_seconds: float

    def __post_init__(self):
        object.__setattr__(self, "mechanisms", tuple(self.mechanisms))
        object.__setattr__(self, "elapsed_seconds", float(self.elapsed_seconds))


@dataclass(frozen=True)
class ResolvedMechanism:
    """Geometry-selected R/P mechanism independent of AAM provenance.

    This is the exact interface required by downstream TS processing.  It can
    be constructed from an ``RPResult`` or a persisted, already validated R/P
    selection without repeating the R/P search.
    """

    mapping: AtomBijection
    broken_bonds: tuple[tuple[int, int], ...]
    formed_bonds: tuple[tuple[int, int], ...]
    core_atoms: tuple[int, ...]

    def __post_init__(self):
        object.__setattr__(self, "broken_bonds", tuple(sorted(
            tuple(sorted(map(int, bond))) for bond in self.broken_bonds)))
        object.__setattr__(self, "formed_bonds", tuple(sorted(
            tuple(sorted(map(int, bond))) for bond in self.formed_bonds)))
        expected_core = tuple(sorted({
            atom for bond in (*self.broken_bonds, *self.formed_bonds)
            for atom in bond
        }))
        core = tuple(sorted(map(int, self.core_atoms)))
        if core != expected_core:
            raise ValueError("resolved mechanism core differs from its events")
        object.__setattr__(self, "core_atoms", core)


@dataclass(frozen=True)
class ReactionContext:
    """Materialized endpoint pair and selected mechanisms consumed by TS."""

    problem: AAMProblem
    config: AAMSearchConfig
    mechanisms: tuple[ResolvedMechanism, ...]

    def __post_init__(self):
        mechanisms = tuple(self.mechanisms)
        if any(item.mapping.degree != self.problem.atom_count
               for item in mechanisms):
            raise ValueError("resolved mechanism has the wrong mapping degree")
        object.__setattr__(self, "mechanisms", mechanisms)


@dataclass(frozen=True)
class VibrationalModes:
    """Frequencies and Cartesian normal modes in one target atom order."""

    frequencies: np.ndarray
    displacements: np.ndarray

    def __post_init__(self):
        frequencies = _readonly_array(self.frequencies, ndim=1)
        displacements = _readonly_array(self.displacements, ndim=3)
        if displacements.shape[0] != frequencies.shape[0]:
            raise ValueError("mode and frequency counts differ")
        if displacements.shape[2] != 3:
            raise ValueError("normal-mode displacements must be Cartesian")
        object.__setattr__(self, "frequencies", frequencies)
        object.__setattr__(self, "displacements", displacements)


@dataclass(frozen=True)
class TransitionStateTarget:
    """A molecular target and its vibrational analysis."""

    molecule: MolecularEndpoint
    vibrations: VibrationalModes
    kind: str = "candidate"

    def __post_init__(self):
        if self.vibrations.displacements.shape[1] != self.molecule.atom_count:
            raise ValueError("normal modes and target atom counts differ")
        object.__setattr__(self, "kind", str(self.kind))


@dataclass(frozen=True)
class AtomAssignment:
    """An injective partial atom assignment with an explicit source domain."""

    pairs: tuple[tuple[int, int], ...]

    def __post_init__(self):
        pairs = tuple(sorted((int(a), int(b)) for a, b in self.pairs))
        if len({a for a, _ in pairs}) != len(pairs):
            raise ValueError("partial assignment repeats a source atom")
        if len({b for _, b in pairs}) != len(pairs):
            raise ValueError("partial assignment repeats a target atom")
        object.__setattr__(self, "pairs", pairs)

    @classmethod
    def from_mapping(cls, mapping):
        return cls(tuple(dict(mapping).items()))

    def as_dict(self):
        return dict(self.pairs)

    @property
    def source_atoms(self):
        return tuple(source for source, _target in self.pairs)


@dataclass(frozen=True)
class CoreAAMBranch:
    """One partial-AAM growth branch and its exact correlated core orbit."""

    representative: AtomAssignment
    hierarchy: AAMHierarchy
    exact_assignments: tuple[AtomAssignment, ...]
    seed_index: int
    search_path: object = None

    def __post_init__(self):
        exact = tuple(self.exact_assignments)
        if not exact:
            raise ValueError("a core branch must realize an exact assignment")
        domain = self.representative.source_atoms
        if any(item.source_atoms != domain for item in exact):
            raise ValueError("core branch assignments disagree on their domain")
        object.__setattr__(self, "exact_assignments", exact)
        object.__setattr__(self, "seed_index", int(self.seed_index))


@dataclass(frozen=True)
class CoreAAMResult:
    """Exact correlated core assignments from endpoint to target."""

    source: MolecularEndpoint
    target: MolecularEndpoint
    core_atoms: tuple[int, ...]
    branches: tuple[CoreAAMBranch, ...]
    assignments: tuple[AtomAssignment, ...]
    elapsed_seconds: float
    capped_seed_count: int = 0

    def __post_init__(self):
        core = tuple(sorted(map(int, self.core_atoms)))
        if any(atom < 0 or atom >= self.source.atom_count for atom in core):
            raise ValueError("core atom lies outside the source molecule")
        assignments = tuple(self.assignments)
        if any(item.source_atoms != core for item in assignments):
            raise ValueError("core assignment has the wrong source domain")
        if any(target < 0 or target >= self.target.atom_count
               for item in assignments for _source, target in item.pairs):
            raise ValueError("core assignment lies outside the target molecule")
        object.__setattr__(self, "core_atoms", core)
        object.__setattr__(self, "branches", tuple(self.branches))
        object.__setattr__(self, "assignments", assignments)
        object.__setattr__(self, "elapsed_seconds", float(self.elapsed_seconds))
        object.__setattr__(self, "capped_seed_count", int(self.capped_seed_count))


@dataclass(frozen=True)
class TSScore:
    """Mode score for one exact mechanism-local core assignment."""

    assignment: AtomAssignment
    sources: frozenset[str]
    score: float
    overlap: float
    wbo_progress: float
    mode_index: int
    frequency: float
    event_terms: tuple[Mapping[str, Any], ...]


@dataclass(frozen=True)
class TSScoringConfig:
    event_weight_power: float = 1.0
    wbo_progress_power: float = 1.0
    prefer_endpoint_consensus: bool = True
    core_assignment_limit: int = 20_000

    def __post_init__(self):
        if self.event_weight_power < 0 or self.wbo_progress_power < 0:
            raise ValueError("TS scoring powers cannot be negative")
        if self.core_assignment_limit < 1:
            raise ValueError("core assignment limit must be positive")


@dataclass(frozen=True)
class TSMechanismResult:
    mechanism: ResolvedMechanism
    target: TransitionStateTarget
    reactant_core_aam: CoreAAMResult | None
    product_core_aam: CoreAAMResult | None
    candidates: tuple[TSScore, ...]
    selected: TSScore | None
    status: str = "scored"
    reason: str | None = None


@dataclass(frozen=True)
class TSResult:
    reaction: ReactionContext
    mechanisms: tuple[TSMechanismResult, ...]
    elapsed_seconds: float

    def __post_init__(self):
        object.__setattr__(self, "mechanisms", tuple(self.mechanisms))
        object.__setattr__(self, "elapsed_seconds", float(self.elapsed_seconds))
