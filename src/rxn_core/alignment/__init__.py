"""Molecule-level WBO graph alignment.

This package owns the public alignment API and branch-level scheduling. Lower
layers handle fragment growth (`rxn_core.growth`) and symmetry-compressed
candidate state (`rxn_core.matcher`).
"""
from __future__ import annotations

from .api import (
    MatchCandidate,
    MatchResult,
    align_from_arrays,
    analyze_alignment,
    cut_edges_above_floor,
    match_wbo_graphs,
)
from .sweep import (
    cut_sweep_items,
    cut_sweep,
    merge_cut_sweep_pools,
    run_cut_sweep_chunk,
    select_min_mechanisms,
)
from .branch import (
    SYM_REPAIR_MAX_EVALS,
    _Branch,
    _alignment_state_signature,
    _chemistry_orbit_signature,
    _chirality_violations,
    _generate_seed_orders,
    _mapping_change_score,
    _orbit_pair,
    find_islands,
    symmetry_repair_mapping,
)
from .index_chirality import (
    IndexChiralityConflict,
    IndexChiralityError,
    IndexChiralitySelection,
    build_index_frames,
    mapping_event_signature,
    normalized_index_orientation,
    select_index_chirality_assignment,
)
