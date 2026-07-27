"""AAM-constrained index-orientation postprocessing for NEB endpoints.

The atom-to-atom mapper remains authoritative.  Starting from one selected
mechanism witness, this module materializes only candidates explicitly encoded
by that witness:

* its correlated, concrete fragment ``alternates``; and
* its closed, complete factorial shuffle blocks.

Open growth states and lossy orbit summaries never authorize a candidate.
Every candidate must preserve the selected mechanism's R-frame bond events.
Among those candidates, all defined persistent local frames must retain their
R-index orientation.  A proper (non-reflecting) Kabsch fit then selects the
lowest-displacement encoded candidate; mapping changes are only a final tie.
"""
from __future__ import annotations

import hashlib
import itertools
import json
import math
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Iterable, Mapping, Sequence

import numpy as np


FORMAT_VERSION = "rxn_core.neb_orientation/v3"


class NebOrientationError(ValueError):
    """Base error for invalid or inconsistent postprocessing input."""


class OrientationConflict(NebOrientationError):
    """No AAM-listed shuffle assignment satisfies the orientation rules."""


@dataclass(frozen=True)
class ShuffleBlock:
    """One complete shuffle pool belonging to the selected AAM witness."""

    block_id: str
    r_atoms: tuple[int, ...]
    p_atoms: tuple[int, ...]
    witness_index: int
    fragment_index: int | None
    block_index: str

    def to_dict(self) -> dict:
        return {
            "id": self.block_id,
            "r_atoms": list(self.r_atoms),
            "p_atoms": list(self.p_atoms),
            "witness_index": self.witness_index,
            "fragment_index": self.fragment_index,
            "block_index": self.block_index,
        }


@dataclass(frozen=True)
class EncodedCandidate:
    """One full mapping materialized from the selected compressed witness."""

    mapping: dict[int, int]
    provenance_paths: tuple[dict, ...]


@dataclass(frozen=True)
class OrientationFrame:
    """One R-index-ordered, center-based local orientation constraint."""

    center: int
    neighbors: tuple[int, ...]
    reactant_orientation: float
    block_ids: tuple[str, ...]

    @property
    def atoms(self) -> tuple[int, int, int, int]:
        return (self.center, *self.neighbors)

    def to_dict(self) -> dict:
        return {
            "center": self.center,
            "neighbors_R_index_order": list(self.neighbors),
            "orientation_model": (
                "affine_four_neighbor_tetrahedron"
                if len(self.neighbors) == 4
                else "center_to_three_neighbor_vectors"
            ),
            "reactant_normalized_orientation": self.reactant_orientation,
            "shuffle_blocks": list(self.block_ids),
        }


@dataclass(frozen=True)
class CandidateFamily:
    """The allowed local mapping family extracted from one mechanism."""

    source_mapping: dict[int, int]
    witness_index: int
    candidates: tuple[EncodedCandidate, ...]
    blocks: tuple[ShuffleBlock, ...]
    fixed_r_atoms: tuple[int, ...]
    discarded_blocks: tuple[dict, ...]


@dataclass(frozen=True)
class RigidTransform:
    """One proper row-vector rigid transform: ``x @ rotation + translation``."""

    rotation: np.ndarray
    translation: np.ndarray
    anchor_r_atoms: tuple[int, ...]
    anchor_rank: int
    rmsd: float

    @property
    def determinant(self) -> float:
        return float(np.linalg.det(self.rotation))

    def apply(self, coords: np.ndarray) -> np.ndarray:
        return np.asarray(coords, dtype=float) @ self.rotation + self.translation


@dataclass(frozen=True)
class NebOrientationResult:
    """Selected NEB mapping and its complete audit information."""

    source_mapping: dict[int, int]
    selected_mapping: dict[int, int]
    selected_candidate: EncodedCandidate
    family: CandidateFamily
    frames: tuple[OrientationFrame, ...]
    transform: RigidTransform
    aligned_product_native_order: np.ndarray
    aligned_product_r_order: np.ndarray
    block_choices: tuple[dict, ...]
    undefined_frame_count: int
    source_violation_count: int
    final_violation_count: int
    max_mutable_displacement: float
    mutable_rmsd: float
    event_signature_unchanged: bool
    broken_bonds_R: tuple[tuple[int, int], ...]
    formed_bonds_R: tuple[tuple[int, int], ...]
    formed_bonds_P: tuple[tuple[int, int], ...]

    def to_dict(self) -> dict:
        changed = [
            {
                "r_atom": r,
                "source_p_atom": self.source_mapping[r],
                "selected_p_atom": self.selected_mapping[r],
            }
            for r in sorted(self.source_mapping)
            if self.source_mapping[r] != self.selected_mapping[r]
        ]
        return {
            "schema_version": FORMAT_VERSION,
            "status": (
                "endpoint_orientation_conflict"
                if self.final_violation_count
                else (
                    "endpoint_orientation_consistent"
                    if self.undefined_frame_count == 0
                    else (
                        "endpoint_orientation_consistent_with_"
                        "undefined_frames"
                    )
                )
            ),
            "source_mapping_RP": {
                str(r): self.source_mapping[r]
                for r in sorted(self.source_mapping)
            },
            "selected_neb_mapping_RP": {
                str(r): self.selected_mapping[r]
                for r in sorted(self.selected_mapping)
            },
            "source_mapping_sha256": mapping_sha256(self.source_mapping),
            "selected_mapping_sha256": mapping_sha256(self.selected_mapping),
            "selected_witness_index": self.family.witness_index,
            "selected_candidate_provenance": list(
                self.selected_candidate.provenance_paths),
            "encoded_candidate_count": len(self.family.candidates),
            "fixed_r_atoms": list(self.family.fixed_r_atoms),
            "mutable_r_atoms": sorted(
                set(self.source_mapping) - set(self.family.fixed_r_atoms)),
            "allowed_shuffle_blocks": [
                block.to_dict() for block in self.family.blocks
            ],
            "discarded_block_records": list(self.family.discarded_blocks),
            "orientation_frames": [frame.to_dict() for frame in self.frames],
            "undefined_frame_count": self.undefined_frame_count,
            "source_orientation_violation_count": (
                self.source_violation_count),
            "final_orientation_violation_count": self.final_violation_count,
            "block_choices": list(self.block_choices),
            "mapping_changes": changed,
            "mapping_change_count": len(changed),
            "selected_broken_bonds_R": [
                list(pair) for pair in self.broken_bonds_R],
            "selected_formed_bonds_R": [
                list(pair) for pair in self.formed_bonds_R],
            "selected_formed_bonds_P": [
                list(pair) for pair in self.formed_bonds_P],
            "geometry_tiebreak": {
                "rule": (
                    "minimum mutable-atom RMSD after one proper "
                    "whole-product fit on fixed atoms; mapping changes "
                    "only break remaining ties"
                ),
                "proper_rotation_determinant": self.transform.determinant,
                "rotation_row_vector": self.transform.rotation.tolist(),
                "translation_angstrom": self.transform.translation.tolist(),
                "anchor_r_atoms": list(self.transform.anchor_r_atoms),
                "anchor_coordinate_rank": self.transform.anchor_rank,
                "anchor_rmsd_angstrom": self.transform.rmsd,
                "maximum_mutable_displacement_angstrom": (
                    self.max_mutable_displacement),
                "mutable_rmsd_angstrom": self.mutable_rmsd,
            },
            "invariants": {
                "complete_bijection": True,
                "elements_preserved": True,
                "fixed_pairs_unchanged": True,
                "shuffle_pool_domains_preserved": True,
                "event_signature_unchanged": self.event_signature_unchanged,
                "selected_mapping_is_encoded_by_selected_aam_witness": (
                    self.family.witness_index >= 0),
                "selected_mapping_is_native_core_assignment": (
                    self.family.witness_index < 0),
                "posthoc_orbit_permutation_used": False,
                "reflection_used": False,
            },
        }


def _int_mapping(mapping: Mapping[int, int]) -> dict[int, int]:
    return {int(r): int(p) for r, p in dict(mapping).items()}


def mapping_sha256(mapping: Mapping[int, int]) -> str:
    normalized = _int_mapping(mapping)
    payload = json.dumps(
        [normalized[r] for r in sorted(normalized)],
        separators=(",", ":"),
    ).encode("ascii")
    return hashlib.sha256(payload).hexdigest()


def validate_mapping(
    mapping: Mapping[int, int],
    elements_R: Sequence[str],
    elements_P: Sequence[str],
) -> dict[int, int]:
    """Return a complete, bijective, element-preserving integer mapping."""
    normalized = _int_mapping(mapping)
    atom_count = len(elements_R)
    if len(elements_P) != atom_count:
        raise NebOrientationError("reactant and product atom counts differ")
    if set(normalized) != set(range(atom_count)):
        raise NebOrientationError(
            "mapping must cover every reactant atom exactly once")
    if set(normalized.values()) != set(range(atom_count)):
        raise NebOrientationError("mapping is not bijective")
    for r, p in normalized.items():
        if str(elements_R[r]) != str(elements_P[p]):
            raise NebOrientationError(
                f"element mismatch at R{r}->P{p}: "
                f"{elements_R[r]} != {elements_P[p]}")
    return normalized


def _is_complete_factorial(block: Mapping, size: int) -> bool:
    value = str(block.get("assignments", "")).replace(" ", "")
    return value == f"{size}!" or value == str(math.factorial(size))


def _exact_fixed_atoms(local_symmetry: Mapping) -> set[int]:
    fixed: set[int] = set()
    for fragment in local_symmetry.get("fragments") or ():
        symmetry = fragment.get("symmetry") or {}
        fixed.update(int(r) for r in symmetry.get("exact_fixed") or ())
    return fixed


def _block_record(block: Mapping, reason: str) -> dict:
    return {
        "reason": reason,
        "source": block.get("source"),
        "fragment_index": block.get("fragment_index"),
        "block_index": block.get("block_index"),
        "r_atoms": [int(v) for v in block.get("r_atoms") or ()],
        "p_atoms": [int(v) for v in block.get("p_atoms") or ()],
        "assignments": block.get("assignments"),
        "extendable": block.get("extendable"),
    }


def _fragment_choice_sets(
    local_symmetry: Mapping,
    source_mapping: Mapping[int, int],
) -> tuple[tuple[dict, ...], ...]:
    """Return concrete primary/alternate choices for each mapped fragment.

    A fragment alternate is a correlated mapping, not a bag of independently
    swappable atoms.  Its R domain and P image set must therefore be identical
    to the fragment primary before it can replace that primary.
    """
    choice_sets: list[tuple[dict, ...]] = []
    for fragment_position, fragment in enumerate(
            local_symmetry.get("fragments") or ()):
        fragment_index = int(fragment.get(
            "fragment_index", fragment_position))
        symmetry = fragment.get("symmetry") or {}
        primary = _int_mapping(symmetry.get("witness") or {})
        if not primary:
            continue
        for r, p in primary.items():
            if source_mapping.get(r) != p:
                raise NebOrientationError(
                    "selected witness local primary does not reconstruct "
                    f"mapping_RP at fragment {fragment_index}, R{r}")

        raw_choices = [{
            "mapping": primary,
            "fragment_index": fragment_index,
            "choice": "primary",
            "alternate_index": None,
            "multiplicity": int(symmetry.get("multiplicity") or 1),
        }]
        for alternate_index, record in enumerate(
                symmetry.get("alternates") or ()):
            alternate = _int_mapping(record.get("witness") or {})
            if set(alternate) != set(primary):
                raise NebOrientationError(
                    "fragment alternate R domain differs from its primary: "
                    f"fragment {fragment_index}, alternate {alternate_index}")
            if set(alternate.values()) != set(primary.values()):
                raise NebOrientationError(
                    "fragment alternate P image set differs from its primary: "
                    f"fragment {fragment_index}, alternate {alternate_index}")
            raw_choices.append({
                "mapping": alternate,
                "fragment_index": fragment_index,
                "choice": "alternate",
                "alternate_index": alternate_index,
                "multiplicity": int(record.get("multiplicity") or 1),
            })

        unique: list[dict] = []
        seen = set()
        for choice in raw_choices:
            key = tuple(sorted(choice["mapping"].items()))
            if key in seen:
                continue
            seen.add(key)
            unique.append(choice)
        choice_sets.append(tuple(unique))
    return tuple(choice_sets)


def _materialize_encoded_candidates(
    source_mapping: Mapping[int, int],
    witness_index: int,
    local_symmetry: Mapping,
    blocks: Sequence[ShuffleBlock],
    elements_R: Sequence[str],
    elements_P: Sequence[str],
    *,
    max_candidate_mappings: int,
) -> tuple[EncodedCandidate, ...]:
    """Expand only concrete alternates and complete blocks in one witness."""
    source = dict(source_mapping)
    fragment_sets = _fragment_choice_sets(local_symmetry, source)
    fragment_count = math.prod(
        len(choices) for choices in fragment_sets) if fragment_sets else 1
    block_count = math.prod(
        math.factorial(len(block.r_atoms)) for block in blocks
    ) if blocks else 1
    encoded_count = fragment_count * block_count
    if encoded_count > max_candidate_mappings:
        raise NebOrientationError(
            "selected AAM witness candidate family exceeds cap: "
            f"{encoded_count} > {max_candidate_mappings}")

    fragment_products = (
        itertools.product(*fragment_sets)
        if fragment_sets else [()]
    )
    block_options = [
        tuple(itertools.permutations(block.p_atoms))
        for block in blocks
    ]
    block_products = (
        tuple(itertools.product(*block_options))
        if block_options else ((),)
    )

    by_mapping: dict[tuple[int, ...], dict] = {}
    for fragment_choices in fragment_products:
        base = dict(source)
        choice_records = []
        for choice in fragment_choices:
            base.update(choice["mapping"])
            choice_records.append({
                "fragment_index": choice["fragment_index"],
                "choice": choice["choice"],
                "alternate_index": choice["alternate_index"],
                "multiplicity": choice["multiplicity"],
            })
        validate_mapping(base, elements_R, elements_P)

        for assignment_tuple in block_products:
            candidate = dict(base)
            assignment_records = []
            for block, p_values in zip(blocks, assignment_tuple):
                candidate.update(zip(block.r_atoms, p_values))
                assignment_records.append({
                    "block_id": block.block_id,
                    "mapping": {
                        str(r): int(p)
                        for r, p in zip(block.r_atoms, p_values)
                    },
                })
            candidate = validate_mapping(
                candidate, elements_R, elements_P)
            key = tuple(candidate[r] for r in range(len(candidate)))
            path = {
                "top_level_witness_index": witness_index,
                "fragment_choices": choice_records,
                "block_assignments": assignment_records,
            }
            entry = by_mapping.setdefault(key, {
                "mapping": candidate,
                "provenance_paths": [],
            })
            entry["provenance_paths"].append(path)

    candidates = tuple(
        EncodedCandidate(
            mapping=entry["mapping"],
            provenance_paths=tuple(entry["provenance_paths"]),
        )
        for key, entry in sorted(by_mapping.items())
    )
    if not any(candidate.mapping == source for candidate in candidates):
        raise NebOrientationError(
            "materialized witness family does not contain mapping_RP")
    return candidates


def build_candidate_family(
    mechanism: Mapping,
    elements_R: Sequence[str],
    elements_P: Sequence[str],
    *,
    max_candidate_mappings: int = 100_000,
) -> CandidateFamily:
    """Materialize the encoded family of the exact selected AAM witness.

    Aggregate ``color_groups`` and lossy ``island_automorph`` /
    ``interbranch`` summaries are never candidate authorization.  Concrete
    fragment alternates remain correlated.  Closed factorial blocks are the
    only compressed permutations that are expanded.  ``extendable`` is a
    growth-state flag and is deliberately ignored.
    """
    source = validate_mapping(
        mechanism.get("mapping_RP") or {}, elements_R, elements_P)
    branch = mechanism.get("branch_symmetry") or {}
    native_index_chirality = branch.get("index_chirality") or {}
    native_status = native_index_chirality.get("status")
    if (
        native_status in {"applied", "conflict"}
        and native_index_chirality.get("selected_mapping_sha256")
        == mapping_sha256(source)
    ):
        provenance_source = (
            "native_core_index_chirality_assignment"
            if native_status == "applied"
            else "native_core_index_chirality_conflict_diagnostic"
        )
        return CandidateFamily(
            source_mapping=source,
            witness_index=-1,
            candidates=(EncodedCandidate(
                mapping=source,
                provenance_paths=({
                    "source": provenance_source,
                    "native_index_chirality_status": native_status,
                    "selected_candidate_id": (
                        native_index_chirality.get(
                            "selected_candidate_id")),
                },),
            ),),
            blocks=(),
            fixed_r_atoms=tuple(sorted(source)),
            discarded_blocks=({
                "reason": (
                    "core_mapping_already_selected; downstream_must_not_"
                    "reselect_atom_assignment"),
            },),
        )
    witnesses = branch.get("witnesses") or ()
    exact = [
        (index, witness)
        for index, witness in enumerate(witnesses)
        if _int_mapping(witness.get("mapping") or {}) == source
    ]
    if len(exact) != 1:
        raise NebOrientationError(
            "expected exactly one branch-symmetry witness matching mapping_RP; "
            f"found {len(exact)}")
    witness_index, witness = exact[0]
    local = witness.get("local_symmetry") or {}
    hard_fixed = _exact_fixed_atoms(local)
    discarded: list[dict] = []
    raw_blocks: list[ShuffleBlock] = []

    flattened = local.get("blocks") or ()
    if flattened:
        discarded.append({
            "reason": "flattened_local_blocks_are_not_candidate_authorization",
            "record_count": len(flattened),
        })

    nested_blocks = []
    for fragment_position, fragment in enumerate(
            local.get("fragments") or ()):
        fragment_index = int(fragment.get(
            "fragment_index", fragment_position))
        symmetry = fragment.get("symmetry") or {}
        for block_position, raw in enumerate(symmetry.get("blocks") or ()):
            nested_blocks.append((
                fragment_index,
                block_position,
                raw,
            ))

    for fragment_index, position, raw in nested_blocks:
        source_kind = raw.get("source")
        if source_kind is not None:
            discarded.append(_block_record(
                raw, "lossy_summary_is_not_candidate_authorization"))
            continue
        r_all = tuple(sorted({int(v) for v in raw.get("r_atoms") or ()}))
        p_all = tuple(sorted({int(v) for v in raw.get("p_atoms") or ()}))
        if bool(raw.get("open")):
            discarded.append(_block_record(
                raw, "open_symmetry_state_is_not_a_complete_shuffle"))
            continue
        if len(r_all) < 2 or len(r_all) != len(p_all):
            discarded.append(_block_record(
                raw, "shuffle_pool_is_not_square_and_nontrivial"))
            continue
        if not _is_complete_factorial(raw, len(r_all)):
            discarded.append(_block_record(
                raw, "shuffle_pool_is_not_complete_factorial"))
            continue
        try:
            images = {source[r] for r in r_all}
        except KeyError:
            discarded.append(_block_record(
                raw, "shuffle_pool_contains_unmapped_reactant_atom"))
            continue
        if images != set(p_all):
            discarded.append(_block_record(
                raw, "shuffle_pool_does_not_match_current_mapping"))
            continue

        fixed_in_block = set(r_all) & hard_fixed
        fixed_p = {source[r] for r in fixed_in_block}
        r_atoms = tuple(r for r in r_all if r not in fixed_in_block)
        p_atoms = tuple(p for p in p_all if p not in fixed_p)
        if len(r_atoms) < 2:
            discarded.append(_block_record(
                raw, "fewer_than_two_mutable_pairs_after_fixed_atoms"))
            continue
        if {source[r] for r in r_atoms} != set(p_atoms):
            discarded.append(_block_record(
                raw, "fixed_atom_reduction_breaks_pool_domain"))
            continue
        for r in r_atoms:
            if any(str(elements_R[r]) != str(elements_P[p])
                   for p in p_atoms):
                raise NebOrientationError(
                    "AAM shuffle pool mixes nuclear elements")
        block_index = str(raw.get("block_index", position))
        block_id = (
            f"w{witness_index}:f{fragment_index}:b{block_index}")
        raw_blocks.append(ShuffleBlock(
            block_id=block_id,
            r_atoms=r_atoms,
            p_atoms=tuple(sorted(p_atoms)),
            witness_index=witness_index,
            fragment_index=fragment_index,
            block_index=block_index,
        ))

    # Remove duplicates and nested residual pools.  A crossing overlap is
    # ambiguous and therefore cannot be silently split into independent pools.
    unique: list[ShuffleBlock] = []
    seen = set()
    for block in sorted(
            raw_blocks,
            key=lambda item: (-len(item.r_atoms), item.r_atoms, item.p_atoms)):
        key = (block.r_atoms, block.p_atoms)
        if key in seen:
            discarded.append({
                **block.to_dict(),
                "reason": "duplicate_shuffle_pool",
            })
            continue
        seen.add(key)
        nested = next((
            kept for kept in unique
            if set(block.r_atoms) < set(kept.r_atoms)
            and set(block.p_atoms) < set(kept.p_atoms)
        ), None)
        if nested is not None:
            discarded.append({
                **block.to_dict(),
                "reason": "nested_residual_pool",
                "retained_parent": nested.block_id,
            })
            continue
        for kept in unique:
            r_overlap = set(block.r_atoms) & set(kept.r_atoms)
            p_overlap = set(block.p_atoms) & set(kept.p_atoms)
            if r_overlap or p_overlap:
                raise NebOrientationError(
                    "crossing AAM shuffle pools cannot be treated as "
                    f"independent: {block.block_id} and {kept.block_id}")
        unique.append(block)

    candidates = _materialize_encoded_candidates(
        source,
        witness_index,
        local,
        unique,
        elements_R,
        elements_P,
        max_candidate_mappings=max_candidate_mappings,
    )
    fixed = tuple(
        r for r in range(len(elements_R))
        if len({candidate.mapping[r] for candidate in candidates}) == 1
    )
    if not hard_fixed.issubset(fixed):
        changed = sorted(hard_fixed - set(fixed))
        raise NebOrientationError(
            "exact-fixed atoms vary in encoded candidate family: "
            + ", ".join(f"R{r}" for r in changed))
    return CandidateFamily(
        source_mapping=source,
        witness_index=witness_index,
        candidates=candidates,
        blocks=tuple(sorted(unique, key=lambda item: item.block_id)),
        fixed_r_atoms=fixed,
        discarded_blocks=tuple(discarded),
    )


def _adjacency(
    wbo: np.ndarray,
    graph_floor: float,
    atom_count: int | None = None,
) -> tuple[tuple[int, ...], ...]:
    matrix = np.asarray(wbo, dtype=float)
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise NebOrientationError("WBO matrix must be square")
    if atom_count is not None and matrix.shape != (atom_count, atom_count):
        raise NebOrientationError(
            f"WBO shape {matrix.shape} does not match {atom_count} atoms")
    if not np.all(np.isfinite(matrix)):
        raise NebOrientationError("WBO matrix contains non-finite values")
    return tuple(
        tuple(int(j) for j in np.flatnonzero(matrix[i] >= graph_floor)
              if int(j) != i)
        for i in range(matrix.shape[0])
    )


def normalized_orientation(
    coords: np.ndarray,
    center: int,
    neighbors: Sequence[int],
) -> float:
    """Signed normalized determinant of three center-to-neighbor vectors."""
    return _orientation_measure(coords, center, neighbors).normalized


@dataclass(frozen=True)
class _OrientationMeasure:
    normalized: float
    determinant: float
    determinant_error_bound: float
    defined: bool
    zero_length: bool


def _orientation_measure(
    coords: np.ndarray,
    center: int,
    neighbors: Sequence[int],
) -> _OrientationMeasure:
    """Evaluate orientation against a scale-aware roundoff bound.

    This is the same determinant predicate used by the native core selector:
    explicit long-double arithmetic and a forward-error bound derived only
    from machine epsilon.  It has no chemistry- or dataset-tuned cutoff.
    """
    neighbors = tuple(int(v) for v in neighbors)
    if len(neighbors) != 3:
        raise NebOrientationError("orientation frame needs three neighbors")
    xyz = np.asarray(coords, dtype=np.longdouble)
    vectors = np.stack(
        [xyz[n] - xyz[int(center)] for n in neighbors], axis=0)
    squared_lengths = np.sum(vectors * vectors, axis=1)
    zero_length = bool(np.any(squared_lengths == 0))
    denominator = np.sqrt(
        squared_lengths[0] * squared_lengths[1] * squared_lengths[2])

    a, b, c = vectors
    positive_terms = (
        a[0] * b[1] * c[2],
        a[1] * b[2] * c[0],
        a[2] * b[0] * c[1],
    )
    negative_terms = (
        a[2] * b[1] * c[0],
        a[1] * b[0] * c[2],
        a[0] * b[2] * c[1],
    )
    determinant = sum(positive_terms) - sum(negative_terms)
    permanent = sum(
        abs(term) for term in (*positive_terms, *negative_terms))
    eps = np.longdouble(np.finfo(np.longdouble).eps)
    operation_count = np.longdouble(16)
    gamma = (
        operation_count * eps
        / (np.longdouble(1) - operation_count * eps)
    )
    error_bound = gamma * max(permanent, denominator)
    defined = bool(
        not zero_length and abs(determinant) > error_bound)
    normalized = (
        np.longdouble(0)
        if zero_length
        else determinant / denominator
    )
    return _OrientationMeasure(
        normalized=float(normalized),
        determinant=float(determinant),
        determinant_error_bound=float(error_bound),
        defined=defined,
        zero_length=zero_length,
    )


def build_orientation_frames(
    family: CandidateFamily,
    coords_R: np.ndarray,
    coords_P: np.ndarray,
    wbo_R: np.ndarray,
    wbo_P: np.ndarray,
    *,
    graph_floor: float = 0.2,
) -> tuple[tuple[OrientationFrame, ...], int]:
    """Build every defined local frame persistent under ``mapping_RP``."""
    xyz_R = np.asarray(coords_R, dtype=float)
    xyz_P = np.asarray(coords_P, dtype=float)
    if xyz_R.shape != xyz_P.shape or xyz_R.shape != (
            len(family.source_mapping), 3):
        raise NebOrientationError("endpoint coordinate shapes do not match")
    atom_count = len(family.source_mapping)
    adjacency_R = _adjacency(wbo_R, graph_floor, atom_count)
    adjacency_P = _adjacency(wbo_P, graph_floor, atom_count)
    block_by_r: dict[int, set[str]] = {}
    for block in family.blocks:
        for r in block.r_atoms:
            block_by_r.setdefault(r, set()).add(block.block_id)

    frames: list[OrientationFrame] = []
    undefined = 0
    mapping = family.source_mapping
    for center, neighbors in enumerate(adjacency_R):
        if len(neighbors) < 3:
            continue
        p_center = mapping[center]
        p_shell = set(adjacency_P[p_center])
        persistent = tuple(
            neighbor for neighbor in sorted(neighbors)
            if mapping[neighbor] in p_shell
        )
        if len(persistent) < 3:
            continue
        for triple in itertools.combinations(persistent, 3):
            atoms = (center, *triple)
            block_ids = sorted({
                block_id
                for atom in atoms
                for block_id in block_by_r.get(atom, ())
            })
            measure_R = _orientation_measure(
                xyz_R, center, triple)
            mapped_triple = tuple(mapping[r] for r in triple)
            measure_P = _orientation_measure(
                xyz_P, p_center, mapped_triple)
            if not measure_R.defined or not measure_P.defined:
                undefined += 1
                continue
            frames.append(OrientationFrame(
                center=int(center),
                neighbors=tuple(int(v) for v in triple),
                reactant_orientation=measure_R.normalized,
                block_ids=tuple(block_ids),
            ))
    return tuple(frames), undefined


def _native_orientation_frames(
    mechanism: Mapping,
    family: CandidateFamily,
) -> tuple[tuple[OrientationFrame, ...], int] | None:
    """Reuse the exact hard-frame set certified by native chirality.

    Native v3 records one affine four-neighbor tetrahedron per hard frame.
    Native v2 used center-to-three-neighbor frames and remains readable for
    older artifacts. Rebuilding either set here would create a second policy,
    so an applied native record is authoritative.
    """
    record = (
        (mechanism.get("branch_symmetry") or {})
        .get("index_chirality") or {}
    )
    if (
        record.get("schema_version") not in {
            "rxn_core.index_chirality/v2",
            "rxn_core.index_chirality/v3",
        }
        or record.get("status") != "applied"
        or record.get("selected_mapping_sha256")
        != mapping_sha256(family.source_mapping)
    ):
        return None
    if "frames" not in record or "undefined_frame_count" not in record:
        raise NebOrientationError(
            "applied native index-chirality metadata is missing its "
            "authoritative frames or undefined-frame count")

    atom_count = len(family.source_mapping)
    expected_neighbor_count = (
        4
        if record.get("schema_version") == "rxn_core.index_chirality/v3"
        else 3
    )
    frames = []
    seen = set()
    for position, raw in enumerate(record.get("frames") or ()):
        try:
            center = int(raw["center_R"])
            neighbors = tuple(
                int(value)
                for value in raw["neighbors_R_index_order"]
            )
            orientation = float(
                raw["reactant_normalized_orientation"])
        except (KeyError, TypeError, ValueError) as exc:
            raise NebOrientationError(
                f"invalid native index-chirality frame {position}") from exc
        if (
            len(neighbors) != expected_neighbor_count
            or len(set(neighbors)) != expected_neighbor_count
            or center in neighbors
            or not (0 <= center < atom_count)
            or any(value < 0 or value >= atom_count for value in neighbors)
            or not math.isfinite(orientation)
            or orientation == 0.0
        ):
            raise NebOrientationError(
                f"invalid native index-chirality frame {position}")
        key = (center, neighbors)
        if key in seen:
            raise NebOrientationError(
                f"duplicate native index-chirality frame {position}")
        seen.add(key)
        frames.append(OrientationFrame(
            center=center,
            neighbors=neighbors,
            reactant_orientation=orientation,
            block_ids=(),
        ))

    try:
        undefined = int(record["undefined_frame_count"])
    except (TypeError, ValueError) as exc:
        raise NebOrientationError(
            "invalid native undefined-frame count") from exc
    if undefined < 0:
        raise NebOrientationError(
            "native undefined-frame count cannot be negative")
    return tuple(frames), undefined


def evaluate_orientations(
    mapping: Mapping[int, int],
    frames: Sequence[OrientationFrame],
    coords_P: np.ndarray,
    *,
    adjacency_P: Sequence[Sequence[int]] | None = None,
) -> tuple[int, tuple[dict, ...]]:
    """Return the number and details of local orientation violations."""
    normalized = _int_mapping(mapping)
    xyz_P = np.asarray(coords_P, dtype=float)
    violations = []
    for frame in frames:
        p_center = normalized[frame.center]
        p_neighbors = tuple(normalized[r] for r in frame.neighbors)
        if adjacency_P is not None:
            shell = set(adjacency_P[p_center])
            missing = [p for p in p_neighbors if p not in shell]
            if missing:
                violations.append({
                    "center": frame.center,
                    "neighbors": list(frame.neighbors),
                    "reason": "mapped_frame_is_not_persistent",
                    "missing_product_neighbors": missing,
                })
                continue
        if len(p_neighbors) == 4:
            # Native v3: affine orientation of the four ligand points.
            measure = _orientation_measure(
                xyz_P, p_neighbors[0], p_neighbors[1:])
        elif len(p_neighbors) == 3:
            # Legacy non-native/v2 frame.
            measure = _orientation_measure(
                xyz_P, p_center, p_neighbors)
        else:
            raise NebOrientationError(
                "orientation frame must have three legacy neighbors or "
                "four native tetrahedral neighbors")
        if not measure.defined:
            violations.append({
                "center": frame.center,
                "neighbors": list(frame.neighbors),
                "reason": "product_frame_degenerate",
                "product_normalized_orientation": measure.normalized,
                "product_determinant": measure.determinant,
                "product_determinant_error_bound": (
                    measure.determinant_error_bound),
            })
            continue
        value = measure.normalized
        if math.copysign(1.0, value) != math.copysign(
                1.0, frame.reactant_orientation):
            violations.append({
                "center": frame.center,
                "neighbors": list(frame.neighbors),
                "reason": "index_orientation_inverted",
                "reactant_normalized_orientation": (
                    frame.reactant_orientation),
                "product_normalized_orientation": value,
            })
    return len(violations), tuple(violations)


def proper_kabsch(
    reference: np.ndarray,
    mobile: np.ndarray,
    *,
    weights: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray, float, int]:
    """Fit ``mobile @ rotation + translation`` to reference without reflection."""
    reference = np.asarray(reference, dtype=float)
    mobile = np.asarray(mobile, dtype=float)
    if reference.shape != mobile.shape or reference.ndim != 2 or (
            reference.shape[1] != 3):
        raise NebOrientationError("Kabsch inputs must have equal (n, 3) shape")
    if len(reference) == 0:
        raise NebOrientationError("at least one fixed anchor is required")
    if weights is None:
        weights = np.ones(len(reference), dtype=float)
    weights = np.asarray(weights, dtype=float)
    if weights.shape != (len(reference),) or np.any(weights <= 0):
        raise NebOrientationError("Kabsch weights must be positive")
    weights = weights / float(weights.sum())
    ref_center = np.sum(reference * weights[:, None], axis=0)
    mob_center = np.sum(mobile * weights[:, None], axis=0)
    ref_zero = reference - ref_center
    mob_zero = mobile - mob_center
    covariance = (mob_zero * weights[:, None]).T @ ref_zero
    u, _singular, vt = np.linalg.svd(covariance)
    rotation = u @ vt
    if np.linalg.det(rotation) < 0:
        u[:, -1] *= -1
        rotation = u @ vt
    translation = ref_center - mob_center @ rotation
    fitted = mobile @ rotation + translation
    rmsd = float(np.sqrt(np.sum(
        weights * np.sum((reference - fitted) ** 2, axis=1))))
    rank = int(np.linalg.matrix_rank(mob_zero, tol=1e-10))
    return rotation, translation, rmsd, rank


def build_fixed_transform(
    family: CandidateFamily,
    elements_R: Sequence[str],
    coords_R: np.ndarray,
    coords_P: np.ndarray,
    *,
    core_atoms: Iterable[int] = (),
) -> RigidTransform:
    """Build one candidate-independent proper fit from fixed mapped atoms."""
    core = {int(v) for v in core_atoms}
    stable = [r for r in family.fixed_r_atoms if r not in core]
    preferred = [r for r in stable if str(elements_R[r]) != "H"]
    choices = (preferred, stable, list(family.fixed_r_atoms))
    anchors: list[int] = []
    for choice in choices:
        if not choice:
            continue
        anchors = list(choice)
        p_atoms = [family.source_mapping[r] for r in anchors]
        centered_R = np.asarray(coords_R)[anchors]
        centered_R = centered_R - centered_R.mean(axis=0)
        centered_P = np.asarray(coords_P)[p_atoms]
        centered_P = centered_P - centered_P.mean(axis=0)
        if (
            len(anchors) >= 3
            and np.linalg.matrix_rank(centered_R, tol=1e-10) >= 2
            and np.linalg.matrix_rank(centered_P, tol=1e-10) >= 2
        ):
            break
    else:
        raise NebOrientationError(
            "at least three non-collinear fixed mapped atoms are required "
            "for the geometry tie-break")
    p_atoms = [family.source_mapping[r] for r in anchors]
    weights = np.array([
        0.1 if str(elements_R[r]) == "H" else 1.0
        for r in anchors
    ], dtype=float)
    rotation, translation, rmsd, rank = proper_kabsch(
        np.asarray(coords_R)[anchors],
        np.asarray(coords_P)[p_atoms],
        weights=weights,
    )
    if np.linalg.det(rotation) <= 0:
        raise NebOrientationError("proper Kabsch unexpectedly returned reflection")
    return RigidTransform(
        rotation=rotation,
        translation=translation,
        anchor_r_atoms=tuple(anchors),
        anchor_rank=rank,
        rmsd=rmsd,
    )


def _event_signature(
    mapping: Mapping[int, int],
    wbo_R: np.ndarray,
    wbo_P: np.ndarray,
    elements_R: Sequence[str],
    *,
    dwbo_threshold: float,
    metal_dwbo_threshold: float | None,
) -> tuple[tuple[tuple[int, int], ...], tuple[tuple[int, int], ...]]:
    from rxn_core.frag import bond_event_threshold

    mapping = _int_mapping(mapping)
    broken = []
    formed = []
    for left in range(len(elements_R)):
        for right in range(left + 1, len(elements_R)):
            threshold = bond_event_threshold(
                elements_R, left, right,
                default_threshold=dwbo_threshold,
                metal_threshold=metal_dwbo_threshold,
            )
            delta = float(
                wbo_R[left, right]
                - wbo_P[mapping[left], mapping[right]]
            )
            if delta >= threshold:
                broken.append((left, right))
            elif -delta >= threshold:
                formed.append((left, right))
    return tuple(broken), tuple(formed)


class _UnionFind:
    def __init__(self, size: int):
        self.parent = list(range(size))

    def find(self, value: int) -> int:
        while self.parent[value] != value:
            self.parent[value] = self.parent[self.parent[value]]
            value = self.parent[value]
        return value

    def union(self, left: int, right: int) -> None:
        a, b = self.find(left), self.find(right)
        if a != b:
            self.parent[b] = a


def _block_components(
    blocks: Sequence[ShuffleBlock],
    frames: Sequence[OrientationFrame],
    adjacency_R: Sequence[Sequence[int]],
    adjacency_P: Sequence[Sequence[int]],
) -> tuple[tuple[int, ...], ...]:
    if not blocks:
        return ()
    by_id = {block.block_id: index for index, block in enumerate(blocks)}
    uf = _UnionFind(len(blocks))
    for frame in frames:
        indexes = [by_id[value] for value in frame.block_ids if value in by_id]
        for index in indexes[1:]:
            uf.union(indexes[0], index)
    for left in range(len(blocks)):
        r_left = set(blocks[left].r_atoms)
        p_left = set(blocks[left].p_atoms)
        for right in range(left + 1, len(blocks)):
            r_right = set(blocks[right].r_atoms)
            p_right = set(blocks[right].p_atoms)
            r_connected = any(
                r_right.intersection(adjacency_R[r]) for r in r_left)
            p_connected = any(
                p_right.intersection(adjacency_P[p]) for p in p_left)
            if r_connected or p_connected:
                uf.union(left, right)
    groups: dict[int, list[int]] = {}
    for index in range(len(blocks)):
        groups.setdefault(uf.find(index), []).append(index)
    return tuple(
        tuple(values)
        for _, values in sorted(groups.items(), key=lambda item: min(item[1]))
    )


def _mapping_displacement(
    mapping: Mapping[int, int],
    r_atoms: Iterable[int],
    coords_R: np.ndarray,
    aligned_P: np.ndarray,
) -> tuple[float, float]:
    atoms = tuple(sorted({int(v) for v in r_atoms}))
    if not atoms:
        return 0.0, 0.0
    distances = np.array([
        np.linalg.norm(coords_R[r] - aligned_P[int(mapping[r])])
        for r in atoms
    ])
    return float(distances.max()), float(np.sqrt(np.mean(distances ** 2)))


def select_neb_mapping(
    family: CandidateFamily,
    frames: Sequence[OrientationFrame],
    elements_R: Sequence[str],
    elements_P: Sequence[str],
    coords_R: np.ndarray,
    coords_P: np.ndarray,
    wbo_R: np.ndarray,
    wbo_P: np.ndarray,
    *,
    core_atoms: Iterable[int] = (),
    dwbo_threshold: float = 0.5,
    metal_dwbo_threshold: float | None = 0.3,
    graph_floor: float = 0.2,
    max_component_candidates: int = 100_000,
    allow_orientation_conflict: bool = False,
) -> NebOrientationResult:
    """Select one orientation-consistent mapping from recorded AAM pools."""
    source = dict(family.source_mapping)
    core = {int(v) for v in core_atoms}
    transform = build_fixed_transform(
        family, elements_R, coords_R, coords_P, core_atoms=core)
    aligned_P = transform.apply(np.asarray(coords_P, dtype=float))
    frame_tuple = tuple(frames)
    atom_count = len(source)
    adjacency_R = _adjacency(wbo_R, graph_floor, atom_count)
    adjacency_P = _adjacency(wbo_P, graph_floor, atom_count)
    source_violations, _ = evaluate_orientations(
        source,
        frame_tuple,
        coords_P,
        adjacency_P=adjacency_P,
    )
    source_signature = _event_signature(
        source, wbo_R, wbo_P, elements_R,
        dwbo_threshold=dwbo_threshold,
        metal_dwbo_threshold=metal_dwbo_threshold,
    )
    mutable_atoms = set(source) - set(family.fixed_r_atoms)
    best = None
    selected_candidate = None
    mechanism_rejections = 0
    orientation_rejections = 0
    for encoded in family.candidates:
        candidate = encoded.mapping
        if _event_signature(
                candidate, wbo_R, wbo_P, elements_R,
                dwbo_threshold=dwbo_threshold,
                metal_dwbo_threshold=metal_dwbo_threshold,
        ) != source_signature:
            mechanism_rejections += 1
            continue
        violation_count, _details = evaluate_orientations(
            candidate,
            frame_tuple,
            coords_P,
            adjacency_P=adjacency_P,
        )
        if violation_count:
            orientation_rejections += 1
        maximum, rmsd = _mapping_displacement(
            candidate, mutable_atoms, coords_R, aligned_P)
        changed = sum(
            candidate[r] != source[r] for r in mutable_atoms)
        assignment_key = tuple(candidate[r] for r in sorted(candidate))
        rank = (
            violation_count,
            round(rmsd, 12),
            round(maximum, 12),
            changed,
            assignment_key,
        )
        if best is None or rank < best:
            best = rank
            selected_candidate = encoded
    if best is None or selected_candidate is None:
        raise OrientationConflict(
            "no encoded selected-witness candidate preserves the mechanism")
    if best[0] != 0 and not allow_orientation_conflict:
        raise OrientationConflict(
            "no encoded selected-witness candidate preserves every defined "
            f"persistent orientation frame; best has {best[0]} violations")
    selected = dict(selected_candidate.mapping)
    choices = [{
        "candidate_source": (
            "selected top-level witness concrete fragment alternates "
            "and closed complete factorial blocks"
        ),
        "encoded_candidate_count": len(family.candidates),
        "mechanism_rejected_candidate_count": mechanism_rejections,
        "orientation_rejected_candidate_count": orientation_rejections,
        "selected_candidate_provenance": list(
            selected_candidate.provenance_paths),
        "orientation_frame_count": len(frame_tuple),
        "rmsd_angstrom": best[1],
        "maximum_displacement_angstrom": best[2],
        "mapping_change_count": best[3],
    }]

    validate_mapping(selected, elements_R, elements_P)
    for r in family.fixed_r_atoms:
        if selected[r] != source[r]:
            raise NebOrientationError(f"fixed mapping changed at R{r}")
    for block in family.blocks:
        if {selected[r] for r in block.r_atoms} != set(block.p_atoms):
            raise NebOrientationError(
                f"selected mapping left AAM pool {block.block_id}")

    final_violations, final_details = evaluate_orientations(
        selected,
        frame_tuple,
        coords_P,
        adjacency_P=adjacency_P,
    )
    if final_violations and not allow_orientation_conflict:
        raise OrientationConflict(
            f"final mapping has {final_violations} orientation violations: "
            f"{final_details[:3]}")
    selected_signature = _event_signature(
        selected, wbo_R, wbo_P, elements_R,
        dwbo_threshold=dwbo_threshold,
        metal_dwbo_threshold=metal_dwbo_threshold,
    )
    if selected_signature != source_signature:
        raise NebOrientationError(
            "selected shuffle changes the mechanism bond-event signature")

    r_order = np.array([selected[r] for r in range(len(selected))], dtype=int)
    product_r_order = aligned_P[r_order]
    maximum, rmsd = _mapping_displacement(
        selected, mutable_atoms, coords_R, aligned_P)
    broken_R, formed_R = selected_signature
    formed_P = tuple(sorted(
        tuple(sorted((selected[left], selected[right])))
        for left, right in formed_R
    ))
    return NebOrientationResult(
        source_mapping=source,
        selected_mapping=selected,
        selected_candidate=selected_candidate,
        family=family,
        frames=frame_tuple,
        transform=transform,
        aligned_product_native_order=aligned_P,
        aligned_product_r_order=product_r_order,
        block_choices=tuple(choices),
        undefined_frame_count=0,
        source_violation_count=source_violations,
        final_violation_count=final_violations,
        max_mutable_displacement=maximum,
        mutable_rmsd=rmsd,
        event_signature_unchanged=True,
        broken_bonds_R=broken_R,
        formed_bonds_R=formed_R,
        formed_bonds_P=formed_P,
    )


def optimize_neb_orientation(
    mechanism: Mapping,
    elements_R: Sequence[str],
    coords_R: np.ndarray,
    wbo_R: np.ndarray,
    elements_P: Sequence[str],
    coords_P: np.ndarray,
    wbo_P: np.ndarray,
    *,
    graph_floor: float = 0.2,
    dwbo_threshold: float = 0.5,
    metal_dwbo_threshold: float | None = 0.3,
    max_component_candidates: int = 100_000,
    allow_orientation_conflict: bool = False,
) -> NebOrientationResult:
    """Convenience entry point for one selected AAM mechanism."""
    family = build_candidate_family(
        mechanism,
        elements_R,
        elements_P,
        max_candidate_mappings=max_component_candidates,
    )
    native_frames = _native_orientation_frames(mechanism, family)
    if native_frames is None:
        frames, undefined = build_orientation_frames(
            family,
            coords_R,
            coords_P,
            wbo_R,
            wbo_P,
            graph_floor=graph_floor,
        )
    else:
        frames, undefined = native_frames
    result = select_neb_mapping(
        family,
        frames,
        elements_R,
        elements_P,
        coords_R,
        coords_P,
        wbo_R,
        wbo_P,
        core_atoms=mechanism.get("core_atoms") or (),
        dwbo_threshold=dwbo_threshold,
        metal_dwbo_threshold=metal_dwbo_threshold,
        graph_floor=graph_floor,
        max_component_candidates=max_component_candidates,
        allow_orientation_conflict=allow_orientation_conflict,
    )
    return replace(result, undefined_frame_count=int(undefined))


def write_result(
    result: NebOrientationResult,
    elements_R: Sequence[str],
    coords_R: np.ndarray,
    output_dir: str | Path,
) -> dict[str, Path]:
    """Write one immutable-source NEB orientation result directory."""
    from rxn_core.chemistry_computations import write_xyz_str

    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    r_path = out / "R.xyz"
    p_path = out / "P_neb_ordered.xyz"
    endpoints_path = out / "neb_endpoints.xyz"
    record_path = out / "neb_orientation.json"
    csv_path = out / "mapping_R_to_P.csv"
    r_text = write_xyz_str(
        elements_R, coords_R, "Reactant endpoint; original R order")
    p_text = write_xyz_str(
        elements_R,
        result.aligned_product_r_order,
        "Product endpoint; AAM-constrained NEB order; proper global fit",
    )
    r_path.write_text(r_text, encoding="utf-8")
    p_path.write_text(p_text, encoding="utf-8")
    endpoints_path.write_text(r_text + p_text, encoding="utf-8")
    record_path.write_text(
        json.dumps(result.to_dict(), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    csv_lines = ["reactant_index,product_index"]
    csv_lines.extend(
        f"{r},{result.selected_mapping[r]}"
        for r in sorted(result.selected_mapping)
    )
    csv_path.write_text("\n".join(csv_lines) + "\n", encoding="utf-8")
    return {
        "reactant": r_path,
        "product": p_path,
        "endpoints": endpoints_path,
        "record": record_path,
        "mapping": csv_path,
    }
