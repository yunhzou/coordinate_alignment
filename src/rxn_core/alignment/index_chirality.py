"""Exact tetrahedral index-parity over native AAM symmetry witnesses.

The selector never discovers product automorphisms or expands factorial image
domains.  Complete AAM witnesses and nested alternates remain atomic.  A live,
closed symmetry block contributes only its odd/even permutation parity, and
coexisting disjoint blocks from one fragment state are solved exactly over
GF(2).  Historical fragment states are never composed with one another.
"""
from __future__ import annotations

import hashlib
import itertools
import json
import math
from dataclasses import dataclass
from typing import Mapping, Sequence

import numpy as np

from ..frag import bond_event_threshold


INDEX_CHIRALITY_SCHEMA = "rxn_core.index_chirality/v3"


class IndexChiralityError(ValueError):
    """Invalid input or inconsistent native AAM symmetry metadata."""


class IndexChiralityConflict(IndexChiralityError):
    """No existing AAM-authorized assignment preserves index chirality."""


@dataclass(frozen=True)
class IndexFrame:
    """One persistent four-neighbor frame ordered by reactant indices."""

    center_R: int
    neighbors_R: tuple[int, int, int, int]
    orientation_R: float
    orientation_P_source: float

    @property
    def frame_id(self) -> str:
        values = "-".join(str(value) for value in self.neighbors_R)
        return f"f:{self.center_R}:{values}"

    def to_dict(self) -> dict:
        return {
            "id": self.frame_id,
            "center_R": self.center_R,
            "neighbors_R_index_order": list(self.neighbors_R),
            "reactant_normalized_orientation": self.orientation_R,
            "source_product_normalized_orientation": (
                self.orientation_P_source),
            "required_sign": 1 if self.orientation_R > 0 else -1,
        }


@dataclass(frozen=True)
class IndexChiralitySelection:
    """Selected assignment plus the bounded candidate audit record."""

    source_mapping: dict[int, int]
    selected_mapping: dict[int, int]
    allowed_assignments: tuple["_AuthorizedAssignment", ...]
    metadata: dict


@dataclass(frozen=True)
class _AuthorizedAssignment:
    """One explicit or fragment-local parity-authorized assignment."""

    candidate_id: str
    mapping: dict[int, int]
    provenance: tuple[dict, ...]
    metadata: dict

    @property
    def action_id(self) -> str:
        """Backward-compatible identifier for callers inspecting selections."""
        return self.candidate_id

    @property
    def cycles_P(self) -> tuple:
        """There is no product-group action in the bounded selector."""
        return ()

    @property
    def generator_word(self) -> tuple:
        """There are no automorphism generators in the bounded selector."""
        return ()


def _int_mapping(mapping: Mapping[int, int]) -> dict[int, int]:
    return {int(r): int(p) for r, p in dict(mapping).items()}


def mapping_sha256(mapping: Mapping[int, int]) -> str:
    """Stable hash of a complete R->P mapping."""
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
    """Return a complete element-preserving bijection or raise."""
    normalized = _int_mapping(mapping)
    atom_count = len(elements_R)
    if len(elements_P) != atom_count:
        raise IndexChiralityError(
            "reactant and product atom counts differ")
    if set(normalized) != set(range(atom_count)):
        raise IndexChiralityError(
            "mapping must cover every reactant atom exactly once")
    if set(normalized.values()) != set(range(atom_count)):
        raise IndexChiralityError("mapping is not bijective")
    for r, p in normalized.items():
        if str(elements_R[r]) != str(elements_P[p]):
            raise IndexChiralityError(
                f"element mismatch at R{r}->P{p}: "
                f"{elements_R[r]} != {elements_P[p]}")
    return normalized


def _add_domain_pair(
    domains: dict[int, set[int]],
    r: int,
    p: int,
    elements_R: Sequence[str],
    elements_P: Sequence[str],
) -> None:
    r = int(r)
    p = int(p)
    if r not in domains or not (0 <= p < len(elements_P)):
        raise IndexChiralityError(
            f"AAM symmetry metadata contains invalid pair R{r}->P{p}")
    if str(elements_R[r]) != str(elements_P[p]):
        raise IndexChiralityError(
            f"AAM symmetry domain mixes elements at R{r}->P{p}")
    domains[r].add(p)


def _nested_fragments(witness: Mapping) -> tuple[Mapping, ...]:
    """Return only fragments owned by one complete branch witness."""
    local = witness.get("local_symmetry") or {}
    return tuple(local.get("fragments") or ())


def _patched_alternates(
    owner_mapping: Mapping[int, int],
    witness: Mapping,
    elements_R: Sequence[str],
    elements_P: Sequence[str],
) -> tuple[tuple[int, int, dict[int, int], tuple[int, ...]], ...]:
    """Patch each nested alternate atomically onto its complete owner."""
    owner = _int_mapping(owner_mapping)
    patched = []
    for fragment_position, fragment in enumerate(_nested_fragments(witness)):
        symmetry = fragment.get("symmetry") or {}
        exact_fixed = {
            int(value) for value in symmetry.get("exact_fixed") or ()}
        for alternate_position, alternate in enumerate(
                symmetry.get("alternates") or ()):
            patch = _int_mapping(alternate.get("witness") or {})
            if not patch:
                continue
            if any(
                r in patch and patch[r] != owner.get(r)
                for r in exact_fixed
            ):
                # Exact-fixed labels belong to the owner witness and cannot be
                # changed by an atomic alternate from the same fragment state.
                continue
            candidate = dict(owner)
            candidate.update(patch)
            try:
                candidate = validate_mapping(
                    candidate, elements_R, elements_P)
            except IndexChiralityError:
                # A fragment record is captured when that island is committed.
                # A later island can consume one of the same product atoms, so
                # replaying an old local alternate onto the final complete
                # witness is not necessarily bijective.  Such a stale row is
                # provenance, not final-assignment authority.
                continue
            patched.append((
                int(fragment_position),
                int(alternate_position),
                candidate,
                tuple(sorted(patch)),
            ))
    return tuple(patched)


def aam_image_domains(
    source_mapping: Mapping[int, int],
    branch_symmetry: Mapping | None,
    elements_R: Sequence[str],
    elements_P: Sequence[str],
    *,
    anchor_map: Mapping[int, int] | None = None,
) -> dict[int, tuple[int, ...]]:
    """Return images observed in explicit complete/atomic AAM witnesses.

    This compatibility helper intentionally ignores ``image_domains`` and
    symmetry blocks.  Those records cannot be multiplied independently into
    complete assignments.
    """
    source = validate_mapping(source_mapping, elements_R, elements_P)
    domains = {r: {p} for r, p in source.items()}
    branch = dict(branch_symmetry or {})
    witnesses = branch.get("witnesses") or ()

    for witness_position, witness in enumerate(witnesses):
        raw_mapping = _int_mapping(witness.get("mapping") or {})
        if not raw_mapping:
            raise IndexChiralityError(
                f"branch witness {witness_position} has no complete mapping")
        try:
            complete = validate_mapping(
                raw_mapping, elements_R, elements_P)
        except IndexChiralityError as exc:
            raise IndexChiralityError(
                f"invalid branch witness {witness_position}: {exc}") from exc
        complete_candidates = [complete]
        complete_candidates.extend(
            patched
            for _, _, patched, _ in _patched_alternates(
                complete, witness, elements_R, elements_P)
        )
        for candidate in complete_candidates:
            for r, p in candidate.items():
                _add_domain_pair(
                    domains, r, p, elements_R, elements_P)

    for raw_r, raw_p in dict(anchor_map or {}).items():
        r, p = int(raw_r), int(raw_p)
        if source.get(r) != p:
            raise IndexChiralityError(
                f"selected AAM mapping violates anchor R{r}->P{p}")
        domains[r] = {p}

    return {
        r: tuple(sorted(values))
        for r, values in sorted(domains.items())
    }


def _adjacency(
    wbo: np.ndarray,
    graph_floor: float,
    atom_count: int,
) -> tuple[tuple[int, ...], ...]:
    matrix = np.asarray(wbo, dtype=float)
    if matrix.shape != (atom_count, atom_count):
        raise IndexChiralityError(
            f"WBO shape {matrix.shape} does not match {atom_count} atoms")
    if not np.all(np.isfinite(matrix)):
        raise IndexChiralityError("WBO matrix contains non-finite values")
    return tuple(
        tuple(
            int(j)
            for j in np.flatnonzero(matrix[i] >= float(graph_floor))
            if int(j) != i
        )
        for i in range(atom_count)
    )


@dataclass(frozen=True)
class _OrientationMeasure:
    normalized: float
    determinant: float
    determinant_error_bound: float
    defined: bool
    zero_length: bool


def _index_orientation_measure(
    coords: np.ndarray,
    center: int,
    neighbors: Sequence[int],
) -> _OrientationMeasure:
    """Evaluate orientation with a scale-aware roundoff bound.

    The determinant is evaluated explicitly in long-double arithmetic.  Its
    error bound is the absolute six-term determinant permanent multiplied by
    a standard ``gamma_n`` factor derived only from machine epsilon.  There is
    no chemistry- or dataset-tuned volume cutoff.
    """
    neighbors = tuple(int(value) for value in neighbors)
    if len(neighbors) != 3:
        raise IndexChiralityError(
            "an index-chirality frame needs exactly three neighbors")
    xyz = np.asarray(coords, dtype=np.longdouble)
    vectors = np.stack(
        [xyz[neighbor] - xyz[int(center)] for neighbor in neighbors],
        axis=0,
    )
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
    # The norm product supplies a nonzero scale even when every determinant
    # term is exactly zero (for example, an exactly planar xy frame).
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


def normalized_index_orientation(
    coords: np.ndarray,
    center: int,
    neighbors: Sequence[int],
) -> float:
    """Signed normalized determinant for three index-ordered vectors."""
    return _index_orientation_measure(
        coords, center, neighbors).normalized


def build_index_frames(
    source_mapping: Mapping[int, int],
    coords_R: np.ndarray,
    coords_P: np.ndarray,
    wbo_R: np.ndarray,
    wbo_P: np.ndarray,
    *,
    graph_floor: float = 0.2,
) -> tuple[tuple[IndexFrame, ...], tuple[dict, ...]]:
    """Build one affine tetrahedral frame per persistent degree-four center."""
    source = _int_mapping(source_mapping)
    atom_count = len(source)
    xyz_R = np.asarray(coords_R, dtype=float)
    xyz_P = np.asarray(coords_P, dtype=float)
    if xyz_R.shape != (atom_count, 3) or xyz_P.shape != (atom_count, 3):
        raise IndexChiralityError(
            "endpoint coordinates must both have shape (n_atoms, 3)")
    if not np.all(np.isfinite(xyz_R)) or not np.all(np.isfinite(xyz_P)):
        raise IndexChiralityError(
            "endpoint coordinates must contain only finite values")
    adjacency_R = _adjacency(wbo_R, graph_floor, atom_count)
    adjacency_P = _adjacency(wbo_P, graph_floor, atom_count)
    frames = []
    undefined = []
    for center_R, neighbors_R in enumerate(adjacency_R):
        # Exactly four graph neighbors define one tetrahedral permutation
        # parity.  Three-coordinate centers can invert without exchanging
        # labels; higher-coordinate centers have no unique tetrahedral parity.
        if len(neighbors_R) != 4:
            continue
        center_P = source[center_R]
        shell_P = set(adjacency_P[center_P])
        neighbors_R = tuple(sorted(int(value) for value in neighbors_R))
        if not all(source[neighbor_R] in shell_P for neighbor_R in neighbors_R):
            continue
        neighbors_P = tuple(source[r] for r in neighbors_R)
        # Affine orientation of the four neighbor points.  The central atom is
        # deliberately not the determinant origin: relabeling the four points
        # then changes the sign exactly by permutation parity.
        measure_R = _index_orientation_measure(
            xyz_R, neighbors_R[0], neighbors_R[1:])
        measure_P = _index_orientation_measure(
            xyz_P, neighbors_P[0], neighbors_P[1:])
        if not measure_R.defined or not measure_P.defined:
            zero_length = measure_R.zero_length or measure_P.zero_length
            undefined.append({
                "center_R": int(center_R),
                "neighbors_R_index_order": [
                    int(value) for value in neighbors_R],
                "reactant_normalized_orientation": measure_R.normalized,
                "source_product_normalized_orientation": (
                    measure_P.normalized),
                "reactant_determinant": measure_R.determinant,
                "source_product_determinant": measure_P.determinant,
                "reactant_determinant_error_bound": (
                    measure_R.determinant_error_bound),
                "source_product_determinant_error_bound": (
                    measure_P.determinant_error_bound),
                "raw_determinant": measure_R.determinant,
                "determinant_error_bound": max(
                    measure_R.determinant_error_bound,
                    measure_P.determinant_error_bound,
                ),
                "reason": (
                    "zero_length_frame_vector"
                    if zero_length
                    else "numerically_indeterminate_orientation"
                ),
            })
            continue
        frames.append(IndexFrame(
            center_R=int(center_R),
            neighbors_R=neighbors_R,
            orientation_R=measure_R.normalized,
            orientation_P_source=measure_P.normalized,
        ))
    return tuple(frames), tuple(undefined)


def index_chirality_violations(
    mapping: Mapping[int, int],
    frames: Sequence[IndexFrame],
    coords_P: np.ndarray,
    wbo_P: np.ndarray,
    *,
    graph_floor: float = 0.2,
) -> tuple[int, tuple[dict, ...]]:
    """Evaluate one complete assignment against precomputed index frames."""
    mapping = _int_mapping(mapping)
    adjacency_P = _adjacency(wbo_P, graph_floor, len(mapping))
    xyz_P = np.asarray(coords_P, dtype=float)
    violations = []
    for frame in frames:
        center_P = mapping[frame.center_R]
        neighbors_P = tuple(mapping[r] for r in frame.neighbors_R)
        shell_P = set(adjacency_P[center_P])
        missing = [p for p in neighbors_P if p not in shell_P]
        if missing:
            violations.append({
                "frame_id": frame.frame_id,
                "reason": "mapped_frame_is_not_persistent",
                "missing_product_neighbors": missing,
            })
            continue
        measure_P = _index_orientation_measure(
            xyz_P, neighbors_P[0], neighbors_P[1:])
        if not measure_P.defined:
            violations.append({
                "frame_id": frame.frame_id,
                "reason": (
                    "mapped_product_frame_has_zero_length_vector"
                    if measure_P.zero_length
                    else "mapped_product_orientation_is_"
                         "numerically_indeterminate"
                ),
                "product_normalized_orientation": measure_P.normalized,
                "product_determinant": measure_P.determinant,
                "product_determinant_error_bound": (
                    measure_P.determinant_error_bound),
            })
            continue
        if math.copysign(1.0, measure_P.normalized) != math.copysign(
                1.0, frame.orientation_R):
            violations.append({
                "frame_id": frame.frame_id,
                "reason": "index_orientation_reversed",
                "reactant_normalized_orientation": frame.orientation_R,
                "product_normalized_orientation": measure_P.normalized,
            })
    return len(violations), tuple(violations)


def mapping_event_signature(
    mapping: Mapping[int, int],
    wbo_R: np.ndarray,
    wbo_P: np.ndarray,
    elements_R: Sequence[str],
    *,
    dwbo_threshold: float = 0.5,
    metal_dwbo_threshold: float | None = 0.3,
) -> tuple[tuple[tuple[int, int], ...], tuple[tuple[int, int], ...]]:
    """Exact broken/formed R-index pairs under one complete assignment."""
    mapping = _int_mapping(mapping)
    broken = []
    formed = []
    for left in range(len(elements_R)):
        for right in range(left + 1, len(elements_R)):
            threshold = bond_event_threshold(
                elements_R,
                left,
                right,
                default_threshold=float(dwbo_threshold),
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


@dataclass(frozen=True)
class _ParityBlock:
    """One independent closed block represented only by permutation parity."""

    block_id: str
    r_atoms: tuple[int, ...]
    p_atoms: tuple[int, ...]
    unanchored_r_atoms: tuple[int, ...]
    odd_swap_R: tuple[int, int]


@dataclass(frozen=True)
class _SeedRoute:
    """One atomic mapping, optionally with one live fragment block state."""

    mapping: dict[int, int]
    provenance: dict
    blocks: tuple[_ParityBlock, ...] = ()
    exact_fixed_r_atoms: tuple[int, ...] = ()


def _mapping_tuple(mapping: Mapping[int, int]) -> tuple[int, ...]:
    normalized = _int_mapping(mapping)
    return tuple(normalized[r] for r in sorted(normalized))


def _provenance_key(record: Mapping) -> str:
    return json.dumps(dict(record), sort_keys=True, separators=(",", ":"))


def _canonical_odd_swap(
    base_mapping: Mapping[int, int],
    r_atoms: Sequence[int],
    fixed_r_atoms: Sequence[int],
    source_mapping: Mapping[int, int],
    elements_R: Sequence[str],
    elements_P: Sequence[str],
) -> tuple[int, int] | None:
    """Return the canonical one-transposition odd representative."""
    base = _int_mapping(base_mapping)
    source = _int_mapping(source_mapping)
    fixed = {int(r) for r in fixed_r_atoms}
    unanchored = tuple(sorted(
        int(r) for r in r_atoms if int(r) not in fixed))
    ranked = []
    for left_R, right_R in itertools.combinations(unanchored, 2):
        left_P = base[left_R]
        right_P = base[right_R]
        if (
            str(elements_R[left_R]) != str(elements_P[right_P])
            or str(elements_R[right_R]) != str(elements_P[left_P])
        ):
            continue
        candidate = dict(base)
        candidate[left_R] = right_P
        candidate[right_R] = left_P
        rank = (
            sum(candidate[r] != source[r] for r in source),
            _mapping_tuple(candidate),
            int(left_R),
            int(right_R),
        )
        ranked.append((rank, (int(left_R), int(right_R))))
    return min(ranked)[1] if ranked else None


def _fragment_parity_blocks(
    fragment: Mapping,
    witness_position: int,
    fragment_position: int,
    base_mapping: Mapping[int, int],
    source_mapping: Mapping[int, int],
    fixed_r_atoms: Sequence[int],
    elements_R: Sequence[str],
    elements_P: Sequence[str],
    atom_count: int,
) -> tuple[tuple[_ParityBlock, ...], tuple[dict, ...], bool]:
    """Validate one fragment-local family of independent parity blocks."""
    symmetry = fragment.get("symmetry") or {}
    records = []
    diagnostics = []
    seen_r: set[int] = set()
    seen_p: set[int] = set()
    base = _int_mapping(base_mapping)
    fixed = {int(r) for r in fixed_r_atoms}
    for block_position, block in enumerate(symmetry.get("blocks") or ()):
        block_id = (
            f"w:{int(witness_position)}/"
            f"f:{int(fragment_position)}/"
            f"b:{int(block_position)}"
        )
        if block.get("source") is not None:
            continue
        if bool(block.get("open")):
            diagnostics.append({
                "fragment_state": (
                    f"w:{int(witness_position)}/f:{int(fragment_position)}"),
                "block_id": block_id,
                "reason": "open_block_has_no_complete_parity",
            })
            continue
        raw_r = tuple(int(value) for value in block.get("r_atoms") or ())
        raw_p = tuple(int(value) for value in block.get("p_atoms") or ())
        invalid = (
            len(raw_r) < 2
            or len(raw_r) != len(raw_p)
            or len(set(raw_r)) != len(raw_r)
            or len(set(raw_p)) != len(raw_p)
            or any(value < 0 or value >= atom_count for value in raw_r)
            or any(value < 0 or value >= atom_count for value in raw_p)
        )
        if invalid:
            diagnostics.append({
                "fragment_state": (
                    f"w:{int(witness_position)}/f:{int(fragment_position)}"),
                "block_id": block_id,
                "reason": "invalid_closed_block",
            })
            return (), tuple(diagnostics), False
        r_atoms = tuple(sorted(raw_r))
        p_atoms = tuple(sorted(raw_p))
        if seen_r.intersection(r_atoms) or seen_p.intersection(p_atoms):
            diagnostics.append({
                "fragment_state": (
                    f"w:{int(witness_position)}/f:{int(fragment_position)}"),
                "block_id": block_id,
                "reason": "fragment_blocks_are_not_pairwise_disjoint",
            })
            return (), tuple(diagnostics), False
        seen_r.update(r_atoms)
        seen_p.update(p_atoms)
        if {base[r] for r in r_atoms} != set(p_atoms):
            diagnostics.append({
                "fragment_state": (
                    f"w:{int(witness_position)}/f:{int(fragment_position)}"),
                "block_id": block_id,
                "reason": "base_mapping_does_not_occupy_block_pool",
            })
            return (), tuple(diagnostics), False
        unanchored = tuple(r for r in r_atoms if r not in fixed)
        odd_swap = _canonical_odd_swap(
            base,
            r_atoms,
            fixed,
            source_mapping,
            elements_R,
            elements_P,
        )
        if len(unanchored) < 2 or odd_swap is None:
            diagnostics.append({
                "fragment_state": (
                    f"w:{int(witness_position)}/f:{int(fragment_position)}"),
                "block_id": block_id,
                "reason": "fewer_than_two_compatible_unanchored_atoms",
            })
            continue
        records.append(_ParityBlock(
            block_id=block_id,
            r_atoms=r_atoms,
            p_atoms=p_atoms,
            unanchored_r_atoms=unanchored,
            odd_swap_R=odd_swap,
        ))
    return tuple(records), tuple(diagnostics), True


def _seed_routes(
    source: Mapping[int, int],
    branch_symmetry: Mapping | None,
    elements_R: Sequence[str],
    elements_P: Sequence[str],
    anchors: Mapping[int, int],
) -> tuple[tuple[_SeedRoute, ...], tuple[dict, ...]]:
    """Build atomic seeds and fragment-local parity routes."""
    routes = [_SeedRoute(
        mapping=dict(source),
        provenance={"kind": "source"},
    )]
    diagnostics = []
    branch = dict(branch_symmetry or {})
    for witness_position, witness in enumerate(
            branch.get("witnesses") or ()):
        raw_mapping = _int_mapping(witness.get("mapping") or {})
        if not raw_mapping:
            raise IndexChiralityError(
                f"branch witness {witness_position} has no complete mapping")
        try:
            complete = validate_mapping(
                raw_mapping, elements_R, elements_P)
        except IndexChiralityError as exc:
            raise IndexChiralityError(
                f"invalid branch witness {witness_position}: {exc}") from exc
        routes.append(_SeedRoute(
            mapping=complete,
            provenance={
                "kind": "branch_witness",
                "witness_index": int(witness_position),
            },
        ))
        for fragment_position, fragment in enumerate(
                _nested_fragments(witness)):
            symmetry = fragment.get("symmetry") or {}
            exact_fixed = tuple(sorted({
                *anchors,
                *(
                    int(value)
                    for value in symmetry.get("exact_fixed") or ()
                ),
            }))
            blocks, block_diagnostics, supported = _fragment_parity_blocks(
                fragment,
                witness_position,
                fragment_position,
                complete,
                source,
                exact_fixed,
                elements_R,
                elements_P,
                len(elements_R),
            )
            diagnostics.extend(block_diagnostics)
            if supported and blocks:
                routes.append(_SeedRoute(
                    mapping=complete,
                    provenance={
                        "kind": "fragment_parity_seed",
                        "witness_index": int(witness_position),
                        "fragment_index": int(fragment_position),
                    },
                    blocks=blocks,
                    exact_fixed_r_atoms=exact_fixed,
                ))
        for (
            fragment_position,
            alternate_position,
            patched,
            patched_r_atoms,
        ) in _patched_alternates(
            complete, witness, elements_R, elements_P
        ):
            routes.append(_SeedRoute(
                mapping=patched,
                provenance={
                    "kind": "nested_alternate",
                    "witness_index": int(witness_position),
                    "fragment_index": int(fragment_position),
                    "alternate_index": int(alternate_position),
                    "patched_r_atoms": [
                        int(value) for value in patched_r_atoms],
                },
            ))
            fragment = _nested_fragments(witness)[fragment_position]
            symmetry = fragment.get("symmetry") or {}
            exact_fixed = tuple(sorted({
                *anchors,
                *(
                    int(value)
                    for value in symmetry.get("exact_fixed") or ()
                ),
            }))
            blocks, block_diagnostics, supported = _fragment_parity_blocks(
                fragment,
                witness_position,
                fragment_position,
                patched,
                source,
                exact_fixed,
                elements_R,
                elements_P,
                len(elements_R),
            )
            diagnostics.extend(block_diagnostics)
            if supported and blocks:
                routes.append(_SeedRoute(
                    mapping=patched,
                    provenance={
                        "kind": "nested_alternate_fragment_parity_seed",
                        "witness_index": int(witness_position),
                        "fragment_index": int(fragment_position),
                        "alternate_index": int(alternate_position),
                        "patched_r_atoms": [
                            int(value) for value in patched_r_atoms],
                    },
                    blocks=blocks,
                    exact_fixed_r_atoms=exact_fixed,
                ))
    return tuple(routes), tuple(diagnostics)


def _solve_gf2(
    equations: Sequence[tuple[int, int]],
    variable_count: int,
) -> tuple[int, ...] | None:
    """Solve a deterministic GF(2) system with every free variable set to 0."""
    rows = [
        int(mask) | ((int(rhs) & 1) << variable_count)
        for mask, rhs in equations
    ]
    pivot_row = 0
    pivots: dict[int, int] = {}
    for column in range(variable_count):
        pivot = next(
            (
                row_index
                for row_index in range(pivot_row, len(rows))
                if (rows[row_index] >> column) & 1
            ),
            None,
        )
        if pivot is None:
            continue
        rows[pivot_row], rows[pivot] = rows[pivot], rows[pivot_row]
        for row_index in range(len(rows)):
            if row_index != pivot_row and ((rows[row_index] >> column) & 1):
                rows[row_index] ^= rows[pivot_row]
        pivots[column] = pivot_row
        pivot_row += 1
    coefficient_mask = (1 << variable_count) - 1
    if any(
        (row & coefficient_mask) == 0
        and ((row >> variable_count) & 1)
        for row in rows
    ):
        return None
    solution = [0] * variable_count
    for column, row_index in pivots.items():
        solution[column] = int(
            (rows[row_index] >> variable_count) & 1)
    return tuple(solution)


def _solve_route_parity(
    route: _SeedRoute,
    frames: Sequence[IndexFrame],
    coords_P: np.ndarray,
    wbo_P: np.ndarray,
    *,
    graph_floor: float,
) -> tuple[dict[int, int] | None, dict]:
    """Solve all frame signs for one coexisting fragment block state."""
    violation_count, details = index_chirality_violations(
        route.mapping,
        frames,
        coords_P,
        wbo_P,
        graph_floor=graph_floor,
    )
    details_by_frame = {
        str(detail["frame_id"]): dict(detail)
        for detail in details
    }
    non_parity = [
        detail
        for detail in details
        if detail.get("reason") != "index_orientation_reversed"
    ]
    audit = {
        "fragment_state": {
            key: value
            for key, value in route.provenance.items()
            if key in {
                "kind",
                "witness_index",
                "fragment_index",
                "alternate_index",
            }
        },
        "exact_fixed_r_atoms": list(route.exact_fixed_r_atoms),
        "variable_blocks": [
            {
                "block_id": block.block_id,
                "r_atoms": list(block.r_atoms),
                "p_atoms": list(block.p_atoms),
                "unanchored_r_atoms": list(block.unanchored_r_atoms),
                "canonical_odd_swap_R": list(block.odd_swap_R),
            }
            for block in route.blocks
        ],
        "base_violation_count": int(violation_count),
    }
    if non_parity:
        audit.update({
            "status": "unsupported_non_parity_frame_failure",
            "details": non_parity,
        })
        return None, audit

    equations = []
    equation_records = []
    for frame in frames:
        neighbor_set = set(frame.neighbors_R)
        mask = 0
        partial = []
        for variable_index, block in enumerate(route.blocks):
            block_set = set(block.r_atoms)
            if frame.center_R in block_set:
                partial.append(block.block_id)
                continue
            intersection = block_set & neighbor_set
            if not intersection:
                continue
            if block_set <= neighbor_set:
                mask |= 1 << variable_index
            else:
                partial.append(block.block_id)
        if partial:
            audit.update({
                "status": "unsupported_partial_frame_block_intersection",
                "frame_id": frame.frame_id,
                "partial_block_ids": sorted(partial),
            })
            return None, audit
        rhs = int(frame.frame_id in details_by_frame)
        equations.append((mask, rhs))
        equation_records.append({
            "frame_id": frame.frame_id,
            "block_ids": [
                route.blocks[index].block_id
                for index in range(len(route.blocks))
                if (mask >> index) & 1
            ],
            "required_parity": rhs,
        })

    solution = _solve_gf2(equations, len(route.blocks))
    audit["equations"] = equation_records
    if solution is None:
        audit["status"] = "inconsistent_parity_equations"
        return None, audit
    mapping = dict(route.mapping)
    applied = []
    for bit, block in zip(solution, route.blocks):
        if not bit:
            continue
        left_R, right_R = block.odd_swap_R
        left_P, right_P = mapping[left_R], mapping[right_R]
        mapping[left_R], mapping[right_R] = right_P, left_P
        applied.append({
            "block_id": block.block_id,
            "parity": 1,
            "swapped_r_atoms": [left_R, right_R],
            "swapped_p_atoms": [left_P, right_P],
        })
    changed_exact_fixed = [
        int(r)
        for r in route.exact_fixed_r_atoms
        if mapping[r] != route.mapping[r]
    ]
    if changed_exact_fixed:
        raise IndexChiralityError(
            "fragment parity solution changed exact-fixed atoms: "
            + ", ".join(f"R{r}" for r in changed_exact_fixed)
        )
    audit.update({
        "status": "solved",
        "solution_bits": list(solution),
        "applied_odd_blocks": applied,
    })
    return mapping, audit


def select_index_chirality_assignment(
    source_mapping: Mapping[int, int],
    branch_symmetry: Mapping | None,
    elements_R: Sequence[str],
    coords_R: np.ndarray,
    wbo_R: np.ndarray,
    elements_P: Sequence[str],
    coords_P: np.ndarray,
    wbo_P: np.ndarray,
    *,
    anchor_map: Mapping[int, int] | None = None,
    graph_floor: float = 0.2,
    dwbo_threshold: float = 0.5,
    metal_dwbo_threshold: float | None = 0.3,
) -> IndexChiralitySelection:
    """Select an exact tetrahedral-parity member of native AAM choices."""
    source = validate_mapping(source_mapping, elements_R, elements_P)
    anchors = _int_mapping(anchor_map or {})
    for r, p in anchors.items():
        if source.get(r) != p:
            raise IndexChiralityError(
                f"selected AAM mapping violates anchor R{r}->P{p}")

    source_signature = mapping_event_signature(
        source,
        wbo_R,
        wbo_P,
        elements_R,
        dwbo_threshold=dwbo_threshold,
        metal_dwbo_threshold=metal_dwbo_threshold,
    )
    event_atoms = {
        int(atom)
        for pair in (*source_signature[0], *source_signature[1])
        for atom in pair
    }
    routes, fragment_state_diagnostics = _seed_routes(
        source, branch_symmetry, elements_R, elements_P, anchors)
    authorized_routes = tuple(
        route
        for route in routes
        if (
            all(route.mapping.get(r) == p for r, p in anchors.items())
            and mapping_event_signature(
                route.mapping,
                wbo_R,
                wbo_P,
                elements_R,
                dwbo_threshold=dwbo_threshold,
                metal_dwbo_threshold=metal_dwbo_threshold,
            ) == source_signature
        )
    )
    switchable_atoms = {
        int(r)
        for route in authorized_routes
        for r in source
        if route.mapping[r] != source[r]
    }
    switchable_atoms.update(
        int(r)
        for route in authorized_routes
        for block in route.blocks
        for r in block.unanchored_r_atoms
    )
    all_frames, raw_undefined_frames = build_index_frames(
        source,
        coords_R,
        coords_P,
        wbo_R,
        wbo_P,
        graph_floor=graph_floor,
    )
    frames = tuple(
        frame
        for frame in all_frames
        if {frame.center_R, *frame.neighbors_R} & switchable_atoms
    )
    immutable_frames = tuple(
        frame
        for frame in all_frames
        if not ({frame.center_R, *frame.neighbors_R} & switchable_atoms)
    )
    undefined_frames = tuple(
        {
            **record,
            "reaction_event_incident": bool(
                {
                    int(record["center_R"]),
                    *(
                        int(value)
                        for value in record["neighbors_R_index_order"]
                    ),
                } & event_atoms
            ),
            "aam_switchable_incident": bool(
                {
                    int(record["center_R"]),
                    *(
                        int(value)
                        for value in record["neighbors_R_index_order"]
                    ),
                } & switchable_atoms
            ),
        }
        for record in raw_undefined_frames
    )
    source_violation_count, _ = index_chirality_violations(
        source,
        frames,
        coords_P,
        wbo_P,
        graph_floor=graph_floor,
    )
    (
        immutable_source_mismatch_count,
        immutable_source_details,
    ) = index_chirality_violations(
        source,
        immutable_frames,
        coords_P,
        wbo_P,
        graph_floor=graph_floor,
    )
    immutable_details_by_frame: dict[str, list[dict]] = {}
    for detail in immutable_source_details:
        immutable_details_by_frame.setdefault(
            str(detail.get("frame_id")), []).append(dict(detail))
    immutable_frame_records = []
    for frame in immutable_frames:
        details = immutable_details_by_frame.get(frame.frame_id, [])
        immutable_frame_records.append({
            **frame.to_dict(),
            "reason": "no_AAM_authorized_switchable_atom_in_frame",
            "reaction_event_incident": bool(
                {frame.center_R, *frame.neighbors_R} & event_atoms),
            "source_index_chirality_mismatch": bool(details),
            "source_mismatch_details": details,
        })

    evaluation_cache: dict[tuple[int, ...], dict] = {}

    def evaluate_candidate(
        raw_mapping: Mapping[int, int],
        provenance: Mapping,
    ) -> dict:
        mapping = validate_mapping(
            raw_mapping, elements_R, elements_P)
        key = _mapping_tuple(mapping)
        existing = evaluation_cache.get(key)
        if existing is not None:
            provenance_key = _provenance_key(provenance)
            known = {
                _provenance_key(record)
                for record in existing["provenance"]
            }
            if provenance_key not in known:
                existing["provenance"].append(dict(provenance))
            return existing

        anchors_match = all(mapping.get(r) == p for r, p in anchors.items())
        signature_matches = mapping_event_signature(
            mapping,
            wbo_R,
            wbo_P,
            elements_R,
            dwbo_threshold=dwbo_threshold,
            metal_dwbo_threshold=metal_dwbo_threshold,
        ) == source_signature
        violation_count, details = index_chirality_violations(
            mapping,
            frames,
            coords_P,
            wbo_P,
            graph_floor=graph_floor,
        )
        violated_frames = {
            detail.get("frame_id")
            for detail in details
            if detail.get("frame_id") is not None
        }
        if not anchors_match:
            status = "anchor_rejected"
        elif not signature_matches:
            status = "mechanism_rejected"
        elif violation_count:
            status = "chirality_rejected"
        else:
            status = "eligible"
        digest = mapping_sha256(mapping)
        row = {
            "candidate_id": f"candidate:{digest}",
            "mapping_sha256": digest,
            "mapping_RP": [
                [int(r), int(mapping[r])] for r in sorted(mapping)],
            "changes_from_source": int(sum(
                mapping[r] != source[r] for r in source)),
            "mapping_valid": True,
            "anchors_match": bool(anchors_match),
            "mechanism_signature_matches": bool(signature_matches),
            "status": status,
            "index_chirality": {
                "violation_count": int(violation_count),
                "frame_flip_bits": [
                    int(frame.frame_id in violated_frames)
                    for frame in frames
                ],
                "violations": list(details),
            },
            "provenance": [dict(provenance)],
            "_mapping": mapping,
        }
        evaluation_cache[key] = row
        return row

    parity_solves = []
    for route in authorized_routes:
        evaluate_candidate(route.mapping, route.provenance)
        if not route.blocks:
            continue
        parity_mapping, parity_audit = _solve_route_parity(
            route,
            frames,
            coords_P,
            wbo_P,
            graph_floor=graph_floor,
        )
        if parity_mapping is not None:
            provenance = {
                "kind": "fragment_parity_solution",
                "seed": dict(route.provenance),
                "solution_bits": list(parity_audit["solution_bits"]),
                "applied_odd_blocks": list(
                    parity_audit["applied_odd_blocks"]),
            }
            evaluated = evaluate_candidate(parity_mapping, provenance)
            parity_audit["candidate_id"] = evaluated["candidate_id"]
            parity_audit["candidate_status"] = evaluated["status"]
        parity_solves.append(parity_audit)

    evaluation_rows = []
    for key in sorted(evaluation_cache):
        row = evaluation_cache[key]
        row["provenance"].sort(key=_provenance_key)
        evaluation_rows.append({
            field: value
            for field, value in row.items()
            if not field.startswith("_")
        })
    zero_rows = sorted(
        (
            row for row in evaluation_cache.values()
            if (
                row["anchors_match"]
                and row["mechanism_signature_matches"]
                and row["index_chirality"]["violation_count"] == 0
            )
        ),
        key=lambda row: (
            int(row["changes_from_source"]),
            _mapping_tuple(row["_mapping"]),
        ),
    )
    if not zero_rows:
        raise IndexChiralityConflict(
            "no existing AAM-authorized atomic or fragment-parity candidate "
            "preserves every defined tetrahedral index frame")

    allowed = tuple(
        _AuthorizedAssignment(
            candidate_id=str(row["candidate_id"]),
            mapping=dict(row["_mapping"]),
            provenance=tuple(dict(value) for value in row["provenance"]),
            metadata={
                "mechanism_signature_matches": True,
                "index_chirality": dict(row["index_chirality"]),
            },
        )
        for row in zero_rows
    )
    selected = allowed[0]
    allowed_rows = [
        {
            field: value
            for field, value in row.items()
            if not field.startswith("_")
        }
        for row in zero_rows
    ]
    selected_changes = [
        {
            "r_atom": int(r),
            "source_p_atom": int(source[r]),
            "selected_p_atom": int(selected.mapping[r]),
        }
        for r in sorted(source)
        if selected.mapping[r] != source[r]
    ]
    metadata = {
        "schema_version": INDEX_CHIRALITY_SCHEMA,
        "policy": "preserve",
        "status": "applied",
        "constraint": (
            "all_defined_persistent_AAM_switchable_tetrahedral_"
            "frame_signs_match_R"),
        "source_mapping_sha256": mapping_sha256(source),
        "selected_mapping_sha256": mapping_sha256(selected.mapping),
        "source_index_chirality_violation_count": source_violation_count,
        "selected_index_chirality_violation_count": 0,
        "switchable_r_atoms": sorted(switchable_atoms),
        "all_defined_frame_count": len(all_frames),
        "persistent_defined_frame_count": (
            len(frames) + len(immutable_frames)),
        "defined_frame_count": len(frames),
        "immutable_frame_count": len(immutable_frames),
        "immutable_source_mismatch_count": int(
            immutable_source_mismatch_count),
        "reaction_event_incident_frame_count": sum(
            bool({frame.center_R, *frame.neighbors_R} & event_atoms)
            for frame in all_frames
        ),
        "undefined_frame_count": len(undefined_frames),
        "frames": [
            {
                **frame.to_dict(),
                "reaction_event_incident": bool(
                    {frame.center_R, *frame.neighbors_R} & event_atoms),
            }
            for frame in frames
        ],
        "immutable_frames": immutable_frame_records,
        "undefined_frames": list(undefined_frames),
        "reaction_event_atoms_R": sorted(event_atoms),
        "reaction_event_signature_R": {
            "broken": [list(pair) for pair in source_signature[0]],
            "formed": [list(pair) for pair in source_signature[1]],
        },
        "candidate_search": {
            "semantics": (
                "explicit_complete_AAM_witnesses_atomic_nested_alternates_"
                "and_exact_fragment_local_closed_block_GF2_parity"
            ),
            "seed_route_count": len(authorized_routes),
            "rejected_seed_route_count": (
                len(routes) - len(authorized_routes)),
            "explicit_witness_seed_count": sum(
                route.provenance.get("kind") == "branch_witness"
                for route in authorized_routes
            ),
            "nested_alternate_seed_count": sum(
                route.provenance.get("kind") == "nested_alternate"
                for route in authorized_routes
            ),
            "fragment_parity_seed_count": sum(
                bool(route.blocks)
                for route in authorized_routes
            ),
            "nested_alternate_fragment_parity_seed_count": sum(
                route.provenance.get("kind")
                == "nested_alternate_fragment_parity_seed"
                for route in authorized_routes
            ),
            "parity_variable_count": sum(
                len(route.blocks) for route in authorized_routes),
            "gf2_equation_count": sum(
                len(record.get("equations") or ())
                for record in parity_solves
            ),
            "gf2_solved_route_count": sum(
                record.get("status") == "solved"
                for record in parity_solves
            ),
            "fragment_state_diagnostics": [
                *fragment_state_diagnostics,
                *parity_solves,
            ],
            "unique_candidate_evaluation_count": len(evaluation_rows),
            "candidate_evaluations": evaluation_rows,
        },
        "allowed_candidate_count": len(allowed),
        "allowed_candidates": allowed_rows,
        "selected_candidate_id": selected.candidate_id,
        "selection_rule": (
            "zero AAM-switchable tetrahedral index-parity violations, then "
            "minimum changes from the native source witness, then canonical "
            "mapping order"
        ),
        "mapping_changes": selected_changes,
        "invariants": {
            "selected_mapping_is_complete_bijection": True,
            "selected_mapping_is_explicit_or_fragment_parity_authorized": True,
            "exact_mechanism_signature_unchanged": True,
            "product_automorphism_generation_used": False,
            "group_closure_used": False,
            "independent_image_domain_expansion_used": False,
            "permutation_or_factorial_expansion_used": False,
            "cross_fragment_block_composition_used": False,
            "fragment_local_GF2_solve_used": bool(parity_solves),
            "free_GF2_variables_are_zero": True,
            "immutable_frames_are_diagnostics_not_constraints": True,
        },
    }
    return IndexChiralitySelection(
        source_mapping=source,
        selected_mapping=dict(selected.mapping),
        allowed_assignments=allowed,
        metadata=metadata,
    )
