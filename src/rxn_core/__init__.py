"""rxn_core: WBO-graph atom alignment + per-mode reactive scoring +
ranked HTML viewer for transition-state initial guesses.

Public API:

    # Atom alignment R <-> X (symmetry-aware WBO graph alignment)
    from rxn_core import align_from_arrays, find_islands, cut_sweep

    # xtb single-point + WBO graph + bond classification
    from rxn_core import (
        run_xtb, parse_xyz, write_xyz_str,
        build_graph, classify_bonds, expand_mapping,
    )

    # Hessian + per-mode reactive features
    from rxn_core import (
        parse_g98_modes,
        core_atoms_in_R_frame, reindex_modes_to_R,
        bond_reaction_vector,
        bond_overlap_per_mode, rxn_overlap_per_mode,
    )

    # IO helpers
    from rxn_core import (
        load_cached_xtb, reindex_to_R_frame,
    )

    # Reusable pipeline stages and BGCP full-view pipeline
    from rxn_core.pipeline import (
        load_endpoint_from_xyz, alignment_inputs_from_xyz,
        smiles_inputs_from_strings,
        load_step_inputs, step_inputs_from_arrays,
        ts_target_from_xyz, ts_target_from_arrays, load_ts_targets,
        discover_mechanisms_from_xyz, discover_mechanisms_from_arrays,
        discover_mechanisms_from_smiles,
        run_rp_stage, write_rp_alignment_files, write_ts_alignment_files,
        run_ts_stage, write_view_stage,
        load_ts_targets_from_specs, process_xyz_stage, process_smiles_stage,
        run_full_pipeline_stage, process_step, main,
    )
"""
from .alignment import (
    align_from_arrays, find_islands, match_wbo_graphs,
    MatchCandidate, MatchResult, cut_edges_above_floor,
    cut_sweep, cut_sweep_items, merge_cut_sweep_pools,
    run_cut_sweep_chunk, select_min_mechanisms,
)
from .chemistry_computations import (
    run_xtb, parse_xyz, write_xyz_str,
    load_cached_xtb, reindex_to_R_frame,
)
from .frag import build_graph, classify_bonds, expand_mapping
from .frag import WeightedGraph, WeightedNode, build_weighted_graph
from .smiles import (
    FormalWBOEndpoint,
    smiles_to_formal_wbo,
    smiles_to_weighted_graph,
)
from .modes import (
    parse_g98_modes,
    core_atoms_in_R_frame, reindex_modes_to_R,
    bond_reaction_vector,
    bond_overlap_per_mode, rxn_overlap_per_mode,
)
from .subgraph import SubgraphMatch, match_weighted_subgraph

_PIPELINE_EXPORTS = {
    "load_endpoint_from_xyz", "alignment_inputs_from_xyz",
    "smiles_inputs_from_strings",
    "load_step_inputs", "step_inputs_from_arrays",
    "ts_target_from_xyz", "ts_target_from_arrays", "load_ts_targets",
    "discover_mechanisms_from_xyz", "discover_mechanisms_from_arrays",
    "discover_mechanisms_from_smiles",
    "rp_cut_work_items", "run_rp_cut_chunk", "merge_rp_cut_chunks",
    "run_rp_stage", "write_rp_alignment_files", "write_ts_alignment_files",
    "run_ts_stage", "merge_ts_stage_chunks", "write_view_stage",
    "load_ts_targets_from_specs", "process_xyz_stage", "process_smiles_stage",
    "run_subgraph_cli", "run_full_pipeline_stage", "process_step",
}


def __getattr__(name):
    if name in _PIPELINE_EXPORTS:
        from . import pipeline
        return getattr(pipeline, name)
    raise AttributeError(name)


__all__ = [
    "align_from_arrays", "find_islands",
    "match_wbo_graphs", "MatchCandidate", "MatchResult",
    "cut_edges_above_floor", "cut_sweep", "cut_sweep_items",
    "merge_cut_sweep_pools", "run_cut_sweep_chunk",
    "select_min_mechanisms",
    "run_xtb", "parse_xyz", "write_xyz_str",
    "build_graph", "classify_bonds", "expand_mapping",
    "WeightedGraph", "WeightedNode", "build_weighted_graph",
    "FormalWBOEndpoint", "smiles_to_formal_wbo", "smiles_to_weighted_graph",
    "SubgraphMatch", "match_weighted_subgraph",
    "parse_g98_modes",
    "core_atoms_in_R_frame", "reindex_modes_to_R",
    "bond_reaction_vector",
    "bond_overlap_per_mode", "rxn_overlap_per_mode",
    "load_cached_xtb", "reindex_to_R_frame",
    "load_endpoint_from_xyz", "alignment_inputs_from_xyz",
    "smiles_inputs_from_strings",
    "load_step_inputs", "step_inputs_from_arrays",
    "ts_target_from_xyz", "ts_target_from_arrays", "load_ts_targets",
    "discover_mechanisms_from_xyz", "discover_mechanisms_from_arrays",
    "discover_mechanisms_from_smiles",
    "rp_cut_work_items", "run_rp_cut_chunk", "merge_rp_cut_chunks",
    "run_rp_stage", "write_rp_alignment_files", "write_ts_alignment_files",
    "run_ts_stage", "merge_ts_stage_chunks", "write_view_stage",
    "load_ts_targets_from_specs", "process_xyz_stage", "process_smiles_stage",
    "run_subgraph_cli", "run_full_pipeline_stage", "process_step",
]

__version__ = "0.1.0"
