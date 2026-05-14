"""rxn_core: WBO-graph atom alignment + per-mode reactive scoring +
ranked HTML viewer for transition-state initial guesses.

Public API:

    # Atom alignment R <-> X (symmetry-aware WBO graph alignment)
    from rxn_core import align_from_arrays, find_islands, cut_sweep, ts_core_pool

    # xtb single-point + WBO graph + bond classification
    from rxn_core import (
        run_xtb, parse_xyz, write_xyz_str,
        build_graph, classify_bonds, expand_mapping,
    )

    # Hessian + per-mode reactive features
    from rxn_core import (
        parse_g98_modes,
        core_atoms_in_R_frame, reindex_modes_to_R,
        reaction_coord_delta, bond_reaction_vector,
        bond_overlap_per_mode, rxn_overlap_per_mode,
    )

    # IO helpers
    from rxn_core import (
        load_cached_xtb, reindex_to_R_frame,
    )

    # Reusable pipeline stages and BGCP full-view pipeline
    from rxn_core.pipeline import (
        load_endpoint_from_xyz, alignment_inputs_from_xyz,
        load_step_inputs, step_inputs_from_arrays,
        ts_target_from_xyz, ts_target_from_arrays, load_ts_targets,
        discover_mechanisms_from_xyz, discover_mechanisms_from_arrays,
        run_rp_stage, write_rp_alignment_files, write_ts_alignment_files,
        run_ts_stage, write_view_stage,
        load_ts_targets_from_specs, process_xyz_stage,
        run_full_pipeline_stage, process_step, main,
    )
"""
from .alignment import (
    align_from_arrays, find_islands, match_wbo_graphs,
    MatchCandidate, MatchResult, cut_edges_above_floor,
    cut_sweep, select_min_mechanisms, ts_core_pool,
)
from .matcher import expand_chemistry_relevant_atoms
from .chemistry_computations import (
    run_xtb, parse_xyz, write_xyz_str,
    load_cached_xtb, reindex_to_R_frame,
)
from .frag import build_graph, classify_bonds, expand_mapping
from .modes import (
    parse_g98_modes,
    core_atoms_in_R_frame, reindex_modes_to_R,
    reaction_coord_delta, bond_reaction_vector,
    bond_overlap_per_mode, rxn_overlap_per_mode,
    kabsch,
)

_PIPELINE_EXPORTS = {
    "load_endpoint_from_xyz", "alignment_inputs_from_xyz",
    "load_step_inputs", "step_inputs_from_arrays",
    "ts_target_from_xyz", "ts_target_from_arrays", "load_ts_targets",
    "discover_mechanisms_from_xyz", "discover_mechanisms_from_arrays",
    "run_rp_stage", "write_rp_alignment_files", "write_ts_alignment_files",
    "run_ts_stage", "write_view_stage",
    "load_ts_targets_from_specs", "process_xyz_stage",
    "run_full_pipeline_stage", "process_step",
}


def __getattr__(name):
    if name in _PIPELINE_EXPORTS:
        from . import pipeline
        return getattr(pipeline, name)
    raise AttributeError(name)


__all__ = [
    "align_from_arrays", "find_islands", "expand_chemistry_relevant_atoms",
    "match_wbo_graphs", "MatchCandidate", "MatchResult",
    "cut_edges_above_floor", "cut_sweep", "select_min_mechanisms",
    "ts_core_pool",
    "run_xtb", "parse_xyz", "write_xyz_str",
    "build_graph", "classify_bonds", "expand_mapping",
    "parse_g98_modes",
    "core_atoms_in_R_frame", "reindex_modes_to_R",
    "reaction_coord_delta", "bond_reaction_vector",
    "bond_overlap_per_mode", "rxn_overlap_per_mode",
    "kabsch",
    "load_cached_xtb", "reindex_to_R_frame",
    "load_endpoint_from_xyz", "alignment_inputs_from_xyz",
    "load_step_inputs", "step_inputs_from_arrays",
    "ts_target_from_xyz", "ts_target_from_arrays", "load_ts_targets",
    "discover_mechanisms_from_xyz", "discover_mechanisms_from_arrays",
    "run_rp_stage", "write_rp_alignment_files", "write_ts_alignment_files",
    "run_ts_stage", "write_view_stage",
    "load_ts_targets_from_specs", "process_xyz_stage",
    "run_full_pipeline_stage", "process_step",
]

__version__ = "0.1.0"
