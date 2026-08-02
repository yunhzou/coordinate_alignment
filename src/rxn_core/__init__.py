"""Typed, composable atom-mapping and transition-state analysis.

The public API has one direction of dependency::

    search_aam -> compile_mapping_families -> select_rp_mappings
                                         -> analyze_transition_state

Serialization, command-line orchestration, and viewers are typed boundary
adapters outside these computational contracts.
"""
from .aam import search_aam
from .analytical import compile_mapping_families
from .core_aam import search_core_assignments
from .domain import (
    AAMMechanism,
    AAMProblem,
    AAMResult,
    AAMSearchConfig,
    AAMSearchMetrics,
    AnalyticalAAMResult,
    AnalyticalBranch,
    AnalyticalMechanism,
    AtomAssignment,
    CoreAAMBranch,
    CoreAAMResult,
    MolecularEndpoint,
    RPMechanism,
    RPResult,
    ReactionContext,
    ResolvedMechanism,
    TSMechanismResult,
    TSResult,
    TSScore,
    TSScoringConfig,
    TransitionStateTarget,
    VibrationalModes,
)
from .rp import align_reaction, select_rp_mappings
from .ts import analyze_transition_state, reaction_context_from_rp
from .artifacts import (
    reaction_from_record, reaction_record, rp_record, ts_record,
    write_rp_bundle, write_ts_record,
)
from .alignment.post_aam import (
    AAMBranch, AAMHierarchy, AtomBijection, AtomPermutation,
    FragmentMatch, PermutationGroup, SymmetryDomain,
)

# Stable graph and chemistry primitives are reusable public components.
from .frag import (
    WeightedGraph,
    WeightedNode,
    build_graph,
    build_weighted_graph,
    classify_bonds,
)
from .modes import bond_overlap_per_mode, rxn_overlap_per_mode
from .subgraph import SubgraphMatch, match_weighted_subgraph

__all__ = [
    "AAMMechanism", "AAMProblem", "AAMResult", "AAMSearchConfig",
    "AAMSearchMetrics", "AnalyticalAAMResult", "AnalyticalBranch",
    "AnalyticalMechanism", "AtomAssignment", "CoreAAMBranch",
    "CoreAAMResult", "MolecularEndpoint", "RPMechanism", "RPResult",
    "ReactionContext", "ResolvedMechanism",
    "TSMechanismResult", "TSResult", "TSScore", "TSScoringConfig",
    "TransitionStateTarget", "VibrationalModes",
    "AAMBranch", "AAMHierarchy", "AtomBijection", "AtomPermutation",
    "FragmentMatch", "PermutationGroup", "SymmetryDomain",
    "search_aam", "compile_mapping_families", "search_core_assignments",
    "select_rp_mappings", "align_reaction", "analyze_transition_state",
    "reaction_context_from_rp",
    "rp_record", "ts_record", "write_rp_bundle", "write_ts_record",
    "reaction_record", "reaction_from_record",
    "WeightedGraph", "WeightedNode", "build_graph", "build_weighted_graph",
    "classify_bonds", "bond_overlap_per_mode", "rxn_overlap_per_mode",
    "SubgraphMatch", "match_weighted_subgraph",
]

__version__ = "0.2.0"
