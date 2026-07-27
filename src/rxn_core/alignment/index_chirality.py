"""Post-process one selected AAM witness for index-orientation consensus."""
from __future__ import annotations

from collections import deque
from dataclasses import dataclass
import math
from typing import Mapping, Sequence

import numpy as np

from ..frag import bond_event_threshold, build_graph
from ..matcher import _nauty_atom_generators


INDEX_CHIRALITY_SCHEMA = "rxn_core.index_chirality/v1"


class IndexChiralityError(ValueError):
    """Invalid index-chirality input or symmetry metadata."""


class IndexChiralityConflict(IndexChiralityError):
    """No allowed final automorphism preserves every mutable frame."""


@dataclass(frozen=True)
class IndexFrame:
    center_R: int
    neighbors_R: tuple[int, int, int, int]
    orientation_R: float
    orientation_P_source: float

    @property
    def frame_id(self):
        shell = "-".join(str(value) for value in self.neighbors_R)
        return f"f:{self.center_R}:{shell}"


@dataclass(frozen=True)
class IndexChiralitySelection:
    source_mapping: dict[int, int]
    selected_mapping: dict[int, int]
    metadata: dict


@dataclass(frozen=True)
class _OrientationMeasure:
    normalized: float
    determinant: float
    determinant_error_bound: float
    defined: bool
    zero_length: bool


def _int_mapping(mapping: Mapping[int, int]) -> dict[int, int]:
    return {int(r): int(p) for r, p in dict(mapping).items()}


def validate_mapping(mapping, elements_R, elements_P):
    mapping = _int_mapping(mapping)
    atom_count = len(elements_R)
    if len(elements_P) != atom_count:
        raise IndexChiralityError("endpoint atom counts differ")
    if set(mapping) != set(range(atom_count)):
        raise IndexChiralityError("mapping is not complete")
    if set(mapping.values()) != set(range(atom_count)):
        raise IndexChiralityError("mapping is not bijective")
    for r, p in mapping.items():
        if str(elements_R[r]) != str(elements_P[p]):
            raise IndexChiralityError(
                f"element mismatch at R{r}->P{p}")
    return mapping


def _adjacency(wbo, graph_floor, atom_count):
    matrix = np.asarray(wbo, dtype=float)
    if matrix.shape != (atom_count, atom_count):
        raise IndexChiralityError("WBO shape does not match atom count")
    return tuple(
        tuple(int(j) for j in np.flatnonzero(matrix[i] >= graph_floor)
              if int(j) != i)
        for i in range(atom_count)
    )


def _orientation_measure(coords, origin, other_points):
    xyz = np.asarray(coords, dtype=np.longdouble)
    vectors = np.stack(
        [xyz[int(point)] - xyz[int(origin)] for point in other_points],
        axis=0,
    )
    squared_lengths = np.sum(vectors * vectors, axis=1)
    zero_length = bool(np.any(squared_lengths == 0))
    denominator = np.sqrt(np.prod(squared_lengths))
    a, b, c = vectors
    positive = (
        a[0] * b[1] * c[2],
        a[1] * b[2] * c[0],
        a[2] * b[0] * c[1],
    )
    negative = (
        a[2] * b[1] * c[0],
        a[1] * b[0] * c[2],
        a[0] * b[2] * c[1],
    )
    determinant = sum(positive) - sum(negative)
    permanent = sum(abs(term) for term in (*positive, *negative))
    eps = np.longdouble(np.finfo(np.longdouble).eps)
    gamma = 16 * eps / (1 - 16 * eps)
    error_bound = gamma * max(permanent, denominator)
    defined = bool(not zero_length and abs(determinant) > error_bound)
    normalized = 0 if zero_length else determinant / denominator
    return _OrientationMeasure(
        normalized=float(normalized),
        determinant=float(determinant),
        determinant_error_bound=float(error_bound),
        defined=defined,
        zero_length=zero_length,
    )


def normalized_index_orientation(coords, neighbors):
    """Orientation of four index-ordered neighbor points."""
    neighbors = tuple(int(value) for value in neighbors)
    if len(neighbors) != 4:
        raise IndexChiralityError("an index frame needs four neighbors")
    return _orientation_measure(
        coords, neighbors[0], neighbors[1:]).normalized


def build_index_frames(source_mapping, coords_R, coords_P, wbo_R, wbo_P,
                       *, graph_floor=0.2):
    """Return persistent, numerically defined degree-four endpoint frames."""
    source = _int_mapping(source_mapping)
    atom_count = len(source)
    adjacency_R = _adjacency(wbo_R, graph_floor, atom_count)
    adjacency_P = _adjacency(wbo_P, graph_floor, atom_count)
    frames = []
    undefined = []
    for center_R, neighbors in enumerate(adjacency_R):
        if len(neighbors) != 4:
            continue
        neighbors_R = tuple(sorted(neighbors))
        center_P = source[center_R]
        neighbors_P = tuple(source[r] for r in neighbors_R)
        if not set(neighbors_P) <= set(adjacency_P[center_P]):
            continue
        measure_R = _orientation_measure(
            coords_R, neighbors_R[0], neighbors_R[1:])
        measure_P = _orientation_measure(
            coords_P, neighbors_P[0], neighbors_P[1:])
        if not measure_R.defined or not measure_P.defined:
            undefined.append({
                "center_R": int(center_R),
                "neighbors_R_index_order": list(neighbors_R),
                "reactant_normalized_orientation": measure_R.normalized,
                "source_product_normalized_orientation": measure_P.normalized,
                "reactant_determinant": measure_R.determinant,
                "source_product_determinant": measure_P.determinant,
                "reason": "numerically_indeterminate_orientation",
            })
            continue
        frames.append(IndexFrame(
            center_R=int(center_R),
            neighbors_R=neighbors_R,
            orientation_R=measure_R.normalized,
            orientation_P_source=measure_P.normalized,
        ))
    return tuple(frames), tuple(undefined)


def mapping_event_signature(mapping, wbo_R, wbo_P, elements_R, *,
                            dwbo_threshold=0.5,
                            metal_dwbo_threshold=0.3):
    """Exact broken/formed R-index pairs for one complete mapping."""
    mapping = _int_mapping(mapping)
    broken = []
    formed = []
    for left in range(len(elements_R)):
        for right in range(left + 1, len(elements_R)):
            threshold = bond_event_threshold(
                elements_R, left, right,
                default_threshold=float(dwbo_threshold),
                metal_threshold=metal_dwbo_threshold)
            delta = float(
                wbo_R[left, right]
                - wbo_P[mapping[left], mapping[right]])
            if delta >= threshold:
                broken.append((left, right))
            elif -delta >= threshold:
                formed.append((left, right))
    return tuple(broken), tuple(formed)


def _selected_blocks(branch_symmetry):
    blocks = []
    seen = set()
    for block in dict(branch_symmetry or {}).get("blocks") or ():
        source = block.get("source") or "sym_block"
        if source not in {
                "sym_block", "alternate_witness",
                "chosen_candidate_automorph",
                "chosen_fragment_automorph"}:
            continue
        r_atoms = tuple(sorted(int(r) for r in block.get("r_atoms") or ()))
        p_atoms = tuple(sorted(int(p) for p in block.get("p_atoms") or ()))
        if not r_atoms or len(p_atoms) < len(r_atoms):
            continue
        key = (r_atoms, p_atoms)
        if key in seen:
            continue
        seen.add(key)
        blocks.append((r_atoms, frozenset(p_atoms), source))
    return tuple(blocks)


def _allowed_automorphism_mappings(source, branch_symmetry, g_P, *,
                                   symmetry_wbo_tol, max_variants):
    """Enumerate the chosen candidate's strict final automorphism orbit."""
    blocks = _selected_blocks(branch_symmetry)
    if not blocks:
        return (dict(source),), (), 0, False
    mutable = sorted({r for r_atoms, _, _ in blocks for r in r_atoms})
    mutable_set = set(mutable)
    fixed = {r: p for r, p in source.items() if r not in mutable_set}
    tag_parts = {}

    def add_tag(atom, tag):
        tag_parts.setdefault(int(atom), []).append(tag)

    for r, p in fixed.items():
        add_tag(p, ("fixed", int(r), int(p)))
    for block_index, (_r_atoms, p_atoms, _source) in enumerate(blocks):
        for p in p_atoms:
            add_tag(p, ("block", int(block_index)))
    generators = _nauty_atom_generators(
        g_P, wbo_tol=float(symmetry_wbo_tol),
        atom_color_tags={
            atom: tuple(parts) for atom, parts in tag_parts.items()
        })
    seed = tuple(source[r] for r in mutable)
    seen = {seed}
    queue = deque([seed])
    truncated = False
    while queue:
        state = queue.popleft()
        for generator in generators:
            candidate = tuple(generator.get(p, p) for p in state)
            if candidate in seen:
                continue
            if len(seen) >= int(max_variants):
                truncated = True
                queue.clear()
                break
            seen.add(candidate)
            queue.append(candidate)
    mappings = []
    for state in sorted(seen):
        candidate = dict(source)
        candidate.update(zip(mutable, state))
        if any(candidate[r] not in p_atoms
               for r_atoms, p_atoms, _ in blocks for r in r_atoms):
            continue
        mappings.append(candidate)
    return tuple(mappings), blocks, len(generators), truncated


def _frame_violations(mapping, frames, coords_P, wbo_P, *, graph_floor):
    adjacency_P = _adjacency(wbo_P, graph_floor, len(mapping))
    violations = []
    for frame in frames:
        center_P = mapping[frame.center_R]
        neighbors_P = tuple(mapping[r] for r in frame.neighbors_R)
        if not set(neighbors_P) <= set(adjacency_P[center_P]):
            violations.append({
                "frame_id": frame.frame_id,
                "reason": "mapped_frame_is_not_persistent",
            })
            continue
        measure = _orientation_measure(
            coords_P, neighbors_P[0], neighbors_P[1:])
        if (not measure.defined or math.copysign(1.0, measure.normalized)
                != math.copysign(1.0, frame.orientation_R)):
            violations.append({
                "frame_id": frame.frame_id,
                "reason": (
                    "index_orientation_reversed" if measure.defined
                    else "mapped_product_orientation_is_undefined"),
                "reactant_normalized_orientation": frame.orientation_R,
                "product_normalized_orientation": measure.normalized,
            })
    return tuple(violations)


def select_index_chirality_assignment(
        source_mapping, branch_symmetry,
        elements_R: Sequence[str], coords_R, wbo_R,
        elements_P: Sequence[str], coords_P, wbo_P, *,
        graph_floor=0.2, symmetry_wbo_tol=0.2,
        dwbo_threshold=0.5, metal_dwbo_threshold=0.3,
        max_variants=20000, require_zero=True, anchor_map=None):
    """Choose an event-preserving final automorphism with matching parity."""
    source = validate_mapping(source_mapping, elements_R, elements_P)
    frames, undefined = build_index_frames(
        source, coords_R, coords_P, wbo_R, wbo_P,
        graph_floor=graph_floor)
    g_P = build_graph(elements_P, wbo_P, bond_cut=graph_floor)
    mappings, blocks, generator_count, orbit_truncated = (
        _allowed_automorphism_mappings(
        source, branch_symmetry, g_P,
        symmetry_wbo_tol=symmetry_wbo_tol,
        max_variants=max_variants))
    source_signature = mapping_event_signature(
        source, wbo_R, wbo_P, elements_R,
        dwbo_threshold=dwbo_threshold,
        metal_dwbo_threshold=metal_dwbo_threshold)
    mappings = tuple(
        mapping for mapping in mappings
        if mapping_event_signature(
            mapping, wbo_R, wbo_P, elements_R,
            dwbo_threshold=dwbo_threshold,
            metal_dwbo_threshold=metal_dwbo_threshold
        ) == source_signature
        and all(mapping[int(r)] == int(p)
                for r, p in dict(anchor_map or {}).items())
    )
    if not mappings:
        raise IndexChiralityConflict(
            "no final automorphism preserves the selected bond events")
    switchable = sorted(
        r for r in source
        if len({mapping[r] for mapping in mappings}) > 1)
    switchable_set = set(switchable)
    active_frames = tuple(
        frame for frame in frames
        if switchable_set.intersection(
            (frame.center_R, *frame.neighbors_R))
    )
    immutable_frames = tuple(
        frame for frame in frames if frame not in active_frames)
    evaluated = []
    for mapping in mappings:
        violations = _frame_violations(
            mapping, active_frames, coords_P, wbo_P,
            graph_floor=graph_floor)
        changed = sum(mapping[r] != source[r] for r in source)
        evaluated.append((
            len(violations), changed,
            tuple(mapping[r] for r in sorted(mapping)),
            mapping, violations,
        ))
    evaluated.sort(key=lambda item: item[:3])
    best = evaluated[0]
    if require_zero and best[0] != 0:
        suffix = (
            f" within max_variants={max_variants}"
            if orbit_truncated else ""
        )
        raise IndexChiralityConflict(
            f"no allowed final automorphism preserves index chirality; "
            f"minimum violations={best[0]}{suffix}")
    selected = dict(best[3])
    source_violations = _frame_violations(
        source, active_frames, coords_P, wbo_P,
        graph_floor=graph_floor)
    immutable_mismatches = _frame_violations(
        source, immutable_frames, coords_P, wbo_P,
        graph_floor=graph_floor)
    metadata = {
        "schema_version": INDEX_CHIRALITY_SCHEMA,
        "policy": "preserve",
        "status": "applied" if selected != source else "already_consistent",
        "candidate_source": "selected_final_candidate_pynauty_automorphisms",
        "candidate_count": len(mappings),
        "candidate_orbit_truncated": bool(orbit_truncated),
        "pynauty_generator_count": int(generator_count),
        "allowed_block_count": len(blocks),
        "switchable_r_atoms": switchable,
        "all_defined_frame_count": len(frames),
        "defined_frame_count": len(active_frames),
        "immutable_frame_count": len(immutable_frames),
        "undefined_frame_count": len(undefined),
        "source_index_chirality_violation_count": len(source_violations),
        "selected_index_chirality_violation_count": int(best[0]),
        "immutable_source_mismatch_count": len(immutable_mismatches),
        "mapping_changes": [
            {
                "r_atom": int(r),
                "source_p_atom": int(source[r]),
                "selected_p_atom": int(selected[r]),
            }
            for r in sorted(source) if source[r] != selected[r]
        ],
        "active_frames": [
            {
                "id": frame.frame_id,
                "center_R": frame.center_R,
                "neighbors_R_index_order": list(frame.neighbors_R),
                "reactant_normalized_orientation": frame.orientation_R,
                "source_product_normalized_orientation": (
                    frame.orientation_P_source),
            }
            for frame in active_frames
        ],
        "undefined_frames": list(undefined),
        "immutable_mismatches": list(immutable_mismatches),
        "event_signature_unchanged": True,
    }
    return IndexChiralitySelection(
        source_mapping=source,
        selected_mapping=selected,
        metadata=metadata,
    )
