"""Symmetry-aware WBO graph matching internals.

The package is intentionally split by responsibility:
- primitives: WBO access and tolerance primitives
- orbits: exact nauty orbit grouping plus WL fallback/debug helper
- state: _SymCand/_SymBlock compressed matching state
- support: witness search inside unresolved symmetry blocks
- dedupe: one-hop boundary-aware state dedupe
- extend: one-atom candidate extension
"""
from __future__ import annotations

from .dedupe import _boundary_signature, _dedup_sym_cands, _p_relation_signature
from .extend import _extend_sym_cands
from .policy import (
    DEFAULT_NODE_POLICY,
    AttributeNodeMatchPolicy,
    CallableNodeMatchPolicy,
    ElementNodeMatchPolicy,
    as_node_match_policy,
)
from .orbits import (
    _atom_tuple_orbit,
    _cand_canon_signature,
    _color_refine_orbits,
    _dedup_cands_by_orbit,
    _group_nodes_by_signature,
    _nauty_atom_generators,
    _nauty_orbits,
    _orbit_wbo_bucket,
)
from .primitives import (
    SYM_SUPPORT_MAX_STATES,
    _edge_wbo,
    _growth_edge_supported,
    _orbit_id,
    _wbo_bucket,
)
from .state import (
    _SymBlock,
    _SymCand,
    _cand_has_open_choice,
    _cand_map,
    _cand_possible_p_atoms,
    _sym_block_assignment_expr,
    _sym_block_indexes,
    _sym_cand_variants,
    _symmetry_state,
)
from .support import (
    _force_sym_value,
    _refine_sym_assignments,
    _r_compatible_with_block,
    _support_witness_for_value,
)
