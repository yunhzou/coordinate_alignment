"""rxn_core: WBO-graph atom alignment + per-mode reactive scoring +
ranked HTML viewer for transition-state initial guesses.

Public API:

    # Atom alignment R <-> X (priority-queue subgraph iso, multi-branch)
    from rxn_core import align_from_arrays, find_islands_pq

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
        load_cached_xtb, fill_unmapped_greedy, reindex_to_R_frame,
    )

    # End-to-end pipeline (matching the rxn-core-pipeline CLI)
    from rxn_core.pipeline import process_step, main
"""
from .pq import align_from_arrays, find_islands_pq, expand_chemistry_relevant_atoms
from .frag import (
    run_xtb, parse_xyz, write_xyz_str,
    build_graph, classify_bonds, expand_mapping,
)
from .modes import (
    parse_g98_modes,
    core_atoms_in_R_frame, reindex_modes_to_R,
    reaction_coord_delta, bond_reaction_vector,
    bond_overlap_per_mode, rxn_overlap_per_mode,
    kabsch,
)
from .align import (
    load_cached_xtb, fill_unmapped_greedy, reindex_to_R_frame,
)

__all__ = [
    "align_from_arrays", "find_islands_pq", "expand_chemistry_relevant_atoms",
    "run_xtb", "parse_xyz", "write_xyz_str",
    "build_graph", "classify_bonds", "expand_mapping",
    "parse_g98_modes",
    "core_atoms_in_R_frame", "reindex_modes_to_R",
    "reaction_coord_delta", "bond_reaction_vector",
    "bond_overlap_per_mode", "rxn_overlap_per_mode",
    "kabsch",
    "load_cached_xtb", "fill_unmapped_greedy", "reindex_to_R_frame",
]

__version__ = "0.1.0"
