"""Single-atom extension for symmetry-compressed fragment matches.

This module is the elementary matching step used by island growth:

    current fragment F  ->  enlarged fragment F union {n}

The input candidate list is *local* to one growing island.  Each `_SymCand`
is one canonical compressed state: it has one concrete witness mapping, plus
symmetry blocks/multiplicity that represent many injective concrete mappings.
The outer ``list[_SymCand]`` stores multiple canonical-distinct local states.

The main function is intentionally private because it operates on internal
growth objects.  The public molecule-level APIs are in ``rxn_core.alignment``.
"""
from __future__ import annotations

from collections import defaultdict
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from typing import Any

from .dedupe import _dedup_sym_cands, _p_relation_signature
from .policy import DEFAULT_NODE_POLICY, as_node_match_policy
from .primitives import _edge_wbo, _growth_edge_supported
from .state import _SymCand, _sym_block_indexes, _sym_cand_variants
from .support import (
    _force_sym_value,
    _refine_sym_assignments,
    _r_compatible_with_block,
    _support_witness_for_value,
)


Node = int
Wbo = float
Support = dict[Node, Node]
SymCandidate = _SymCand | Mapping[Node, Node]
OrbitMap = Mapping[Node, int] | None
EdgeKey = tuple[Node, Node]
TargetEntry = tuple[Node, Support, bool]




def _extend_sym_cands(
    cands: Iterable[SymCandidate],
    fragment_old: set[Node],
    n: Node,
    g_R,
    g_P,
    mapping: Mapping[Node, Node],
    iso_tol: float,
    islands_R: Mapping[Node, int] | None,
    p_orbits: OrbitMap = None,
    r_orbits: OrbitMap = None,
    deferred_edges: Iterable[EdgeKey] = (),
    anchor_u: Node | None = None,
    anchor_wbo: Wbo | None = None,
    dedupe_edges: Iterable[EdgeKey] | None = None,
    node_policy=None,
) -> list[_SymCand]:
    """Symmetry-compressed incremental extension.

    This is the elementary one-atom fragment-matching step:

        old fragment F  ->  new fragment F union {n}

    Parameters
    ----------
    cands
        Parallel canonical states for the old fragment.  Each `_SymCand`
        stores one witness plus compact automorphism data; the outer list stores
        multiple canonical-distinct local possibilities for this island.
    fragment_old
        R atom set for the old fragment F.
    n
        New R atom being added to F.
    g_R, g_P
        Reactant/target WBO graphs.  They carry both node labels and a
        full WBO matrix at `graph["wbo_matrix"]`.
    mapping
        Already locked global R->P mapping from other islands.  Free extension
        derives the inverse internally only to prevent reusing locked P atoms;
        merge extension checks the pre-locked image of `n`.
    iso_tol
        Active R-pair WBO tolerance.  For every atom `u` in the old fragment
        with an active R-side graph edge to `n`, a proposed image `n -> v` is
        valid only if some supported witness satisfies
        `WBO_P[v,image(u)] >= graph_floor` and
        `abs(WBO_R[n,u] - WBO_P[v,image(u)]) <= iso_tol`.
    islands_R
        Locked-island labels for R atoms.  If `n` is already mapped and belongs
        to an island, extension merges the whole not-yet-fragment part of that
        island, not just one atom.
    p_orbits, r_orbits
        Symmetry orbit maps, usually from `_nauty_orbits(...)` with the
        configured symmetry WBO tolerance.
        They are compression keys only; active R-pair validity is still checked
        against exact WBO values by `iso_tol`.
    deferred_edges
        Previously failed frontier edges kept as one-hop boundary evidence for
        dedupe.  They are not chemistry cuts and are not removed from the WBO
        graph.
    anchor_u, anchor_wbo
        The popped growth edge `anchor_u -> n` and its WBO.  The edge is not a
        special chemistry rule; this just keeps trace/debug and support checks
        tied to the exact heap edge that triggered the extension.
    dedupe_edges
        Optional boundary-edge set to use for deduping the child states.  The
        growth loop passes current deferred edges plus current one-hop frontier
        edges so symmetry dedupe does not erase distinguishable boundary states.
    node_policy
        Node-level compatibility policy.  The default is same element.  Custom
        policies can admit target nodes by electronic or other descriptors
        before the normal WBO/``iso_tol`` edge verifier is applied.

    Returns
    -------
    list[_SymCand]
        Canonical-distinct compressed matches for `fragment_old union {n}`.
        Each returned `_SymCand` may represent many concrete injective
        assignments through symmetry blocks and multiplicity/alternates.
    """
    node_policy = as_node_match_policy(node_policy)
    ctx = _make_extension_context(
        fragment_old, n, g_R, g_P, mapping, iso_tol, islands_R,
        p_orbits, r_orbits, deferred_edges, anchor_u, anchor_wbo,
        dedupe_edges, node_policy)
    if ctx is None:
        return []

    primary_cands, alternate_cands = _expand_primary_and_alternates(cands)
    children: list[_SymCand] = []
    for cand in primary_cands:
        children.extend(_extend_one_candidate(cand, ctx))

    if not children and alternate_cands:
        return _extend_sym_cands(
            alternate_cands,
            set(ctx.fragment_old),
            ctx.n,
            ctx.g_R,
            ctx.g_P,
            ctx.mapping,
            ctx.iso_tol,
            ctx.islands_R,
            p_orbits=ctx.p_orbits,
            r_orbits=ctx.r_orbits,
            deferred_edges=ctx.deferred_edges,
            anchor_u=ctx.anchor_u,
            anchor_wbo=ctx.anchor_wbo,
            dedupe_edges=ctx.dedupe_edges,
            node_policy=ctx.node_policy,
        )
    return _dedupe_children(children, ctx)



@dataclass(frozen=True)
class _ExtensionContext:
    """All shared state for one ``F -> F union {n}`` extension.

    Keeping this as an explicit object makes the rest of the file easier to
    audit: helper functions either transform one `_SymCand`, or group/dedupe
    the resulting children.  Nothing here changes matching semantics.
    """

    fragment_old: frozenset[Node]
    n: Node
    n_element: str
    g_R: Any
    g_P: Any
    mapping: Mapping[Node, Node]
    locked_p_atoms: frozenset[Node]
    iso_tol: float
    islands_R: Mapping[Node, int] | None
    p_orbits: OrbitMap
    r_orbits: OrbitMap
    deferred_edges: tuple[EdgeKey, ...]
    anchor_u: Node | None
    anchor_wbo: Wbo | None
    dedupe_edges: tuple[EdgeKey, ...] | None
    bonded_in_frag: tuple[Node, ...]
    r_wbos: tuple[tuple[Node, Wbo], ...]
    strict_r_wbos: Mapping[Node, Wbo]
    island_atoms: tuple[Node, ...]
    node_policy: Any = DEFAULT_NODE_POLICY

    @property
    def is_merge(self) -> bool:
        """True when the new atom already has a locked global image."""
        return self.n in self.mapping

    @property
    def sig_fragment(self) -> set[Node]:
        """Fragment used in dedupe signatures after this extension attempt."""
        return set(self.fragment_old) | {self.n, *self.island_atoms}

    @property
    def boundary_edges(self) -> Iterable[EdgeKey]:
        """Boundary evidence used when deduping child canonical states."""
        return self.deferred_edges if self.dedupe_edges is None else self.dedupe_edges


def _locked_island_atoms(
    n: Node,
    fragment_old: set[Node],
    mapping: Mapping[Node, Node],
    islands_R: Mapping[Node, int] | None,
) -> tuple[Node, ...]:
    """Atoms that must be merged when ``n`` is already locked.

    If ``n`` belongs to an existing island, growth should absorb the whole
    not-yet-fragment part of that island.  Otherwise the merge is just ``n``.
    For a free atom, the only new atom is ``n``.
    """
    if n not in mapping:
        return (n,)
    if islands_R is None or n not in islands_R:
        return (n,)
    iid = islands_R[n]
    return tuple(
        r for r, island_id in islands_R.items()
        if island_id == iid and r not in fragment_old
    )


def _strict_growth_wbos(
    fragment_old: set[Node],
    anchor_u: Node | None,
    anchor_wbo: Wbo | None,
    g_R,
    n: Node,
) -> dict[Node, Wbo]:
    """Exact WBO for the heap edge that triggered this extension.

    The heap edge is not a special chemistry rule.  It is one active R-side
    pair in the extension vector.  Keeping it explicit makes traces and support
    checks use the exact value popped from the heap.
    """
    if anchor_u is None or anchor_u not in fragment_old:
        return {}
    return {
        anchor_u: _edge_wbo(g_R, anchor_u, n)
        if anchor_wbo is None else anchor_wbo
    }


def _active_fragment_neighbors(
    fragment_old: set[Node],
    g_R,
    n: Node,
) -> tuple[Node, ...]:
    """Fragment atoms whose R-side pair to ``n`` is in the active graph.

    Extension growth is local on the R graph.  The validity vector is therefore
    scoped to active R-side pairs, rather than every zero-WBO pair to the
    current fragment.
    """
    return tuple(sorted(u for u in fragment_old if g_R.has_edge(u, n)))


def _make_extension_context(
    fragment_old: set[Node],
    n: Node,
    g_R,
    g_P,
    mapping: Mapping[Node, Node],
    iso_tol: float,
    islands_R: Mapping[Node, int] | None,
    p_orbits: OrbitMap,
    r_orbits: OrbitMap,
    deferred_edges: Iterable[EdgeKey],
    anchor_u: Node | None,
    anchor_wbo: Wbo | None,
    dedupe_edges: Iterable[EdgeKey] | None,
    node_policy,
) -> _ExtensionContext | None:
    """Collect repeated extension inputs into one typed context."""
    bonded_in_frag = _active_fragment_neighbors(fragment_old, g_R, n)
    if not bonded_in_frag:
        return None
    locked_p_atoms = frozenset(mapping.values())
    return _ExtensionContext(
        fragment_old=frozenset(fragment_old),
        n=n,
        n_element=g_R.nodes[n].get('element'),
        g_R=g_R,
        g_P=g_P,
        mapping=mapping,
        locked_p_atoms=locked_p_atoms,
        iso_tol=iso_tol,
        islands_R=islands_R,
        p_orbits=p_orbits,
        r_orbits=r_orbits,
        deferred_edges=tuple(tuple(sorted(e)) for e in deferred_edges),
        anchor_u=anchor_u,
        anchor_wbo=anchor_wbo,
        dedupe_edges=(None if dedupe_edges is None
                      else tuple(tuple(sorted(e)) for e in dedupe_edges)),
        bonded_in_frag=bonded_in_frag,
        r_wbos=tuple((u, _edge_wbo(g_R, u, n)) for u in bonded_in_frag),
        strict_r_wbos=_strict_growth_wbos(
            fragment_old, anchor_u, anchor_wbo, g_R, n),
        island_atoms=_locked_island_atoms(n, fragment_old, mapping, islands_R),
        node_policy=node_policy,
    )


def _expand_primary_and_alternates(
    cands: Iterable[SymCandidate],
) -> tuple[list[_SymCand], list[_SymCand]]:
    """Normalize raw mappings and split alternate witnesses.

    `_SymCand.alternates` holds concrete witnesses collapsed into the same
    canonical state.  Normally the primary witness is enough.  If all primary
    witnesses fail to extend, `_extend_sym_cands` retries the alternates before
    declaring the extension impossible.
    """
    expanded: list[_SymCand] = []
    alternates: list[_SymCand] = []
    for raw_cand in cands:
        cand = raw_cand if isinstance(raw_cand, _SymCand) else _SymCand(raw_cand)
        if cand.alternates:
            if cand.multiplicity:
                expanded.append(cand.without_alternates(cand.multiplicity))
            alternates.extend(_sym_cand_variants(cand)[1:])
        else:
            expanded.append(cand)
    return expanded, alternates


def _candidate_covers_fragment(cand: _SymCand, ctx: _ExtensionContext) -> bool:
    """A candidate can extend only if it maps every old fragment atom."""
    return all(u in cand.mapping for u in ctx.bonded_in_frag)


def _supported_value(
    cand: _SymCand,
    ctx: _ExtensionContext,
    v: Node,
    join_block_idx: int | None,
) -> Support | None:
    """Return support assignments proving ``ctx.n -> v`` is valid.

    This is where active R-pair WBO validity enters the extension step.  The
    helper searches inside unresolved symmetry blocks instead of trusting the
    stored witness, so an arbitrary witness cannot incorrectly reject a valid
    correlated assignment.
    """
    return _support_witness_for_value(
        cand,
        ctx.n,
        v,
        list(ctx.bonded_in_frag),
        list(ctx.r_wbos),
        ctx.g_P,
        ctx.iso_tol,
        join_block_idx=join_block_idx,
        strict_r_wbos=ctx.strict_r_wbos,
    )


def _force_required_image(
    cand: _SymCand,
    ctx: _ExtensionContext,
    r: Node,
    p: Node,
):
    """Force one R atom to one P atom while respecting symmetry blocks."""
    return _force_sym_value(
        cand, r, p, set(ctx.fragment_old),
        ctx.g_R, ctx.r_orbits, ctx.p_orbits,
    )


def _island_merge_wbo_consistent(cand: _SymCand, ctx: _ExtensionContext) -> bool:
    """Check full WBO consistency after merging locked island atoms.

    The merge path may add more than just ``ctx.n``.  Every added island atom
    must be WBO-compatible with every materialized atom already in the candidate.
    """
    base = cand.materialize()
    graph_floor = float(ctx.g_P.graph.get("bond_cut", 0.2))
    check_set = set(base.keys())
    for r in ctx.island_atoms:
        for r2 in sorted(check_set):
            if r2 == r:
                continue
            if r >= r2 and r2 in ctx.island_atoms:
                continue
            if not ctx.g_R.has_edge(r, r2):
                continue
            w_r = _edge_wbo(ctx.g_R, r, r2)
            p, p2 = base[r], base[r2]
            if not _growth_edge_supported(
                    w_r, _edge_wbo(ctx.g_P, p, p2),
                    ctx.iso_tol, graph_floor):
                return False
    return True


def _extend_locked_merge(cand: _SymCand, ctx: _ExtensionContext) -> _SymCand | None:
    """Extend into an atom that already has a locked branch-level image.

    This is the ``merge_island`` path.  The target image is fixed by the
    current branch, but the current island candidate may still need to reshuffle
    a symmetry block to support that fixed image.
    """
    v_n = ctx.mapping[ctx.n]
    _, p_to_block = _sym_block_indexes(cand)
    join_idx = p_to_block.get(v_n)
    if join_idx is not None and not _r_compatible_with_block(
            cand, join_idx, ctx.n, set(ctx.fragment_old), ctx.g_R, ctx.r_orbits):
        return None

    support = _supported_value(cand, ctx, v_n, join_idx)
    if support is None:
        return None

    child = _force_required_image(cand, ctx, ctx.n, v_n)
    if child is None:
        return None
    support = dict(support)
    support[ctx.n] = v_n
    child = child.with_witness(support) if isinstance(child, _SymCand) else child
    if child is None:
        return None

    for r in ctx.island_atoms:
        if r == ctx.n:
            continue
        child = _force_required_image(child, ctx, r, ctx.mapping[r])
        if child is None:
            return None

    return child if _island_merge_wbo_consistent(child, ctx) else None


def _target_join_info(
    cand: _SymCand,
    ctx: _ExtensionContext,
    v: Node,
) -> tuple[int | None, bool]:
    """Return the open block containing ``v`` and whether it can grow freely.

    The compatibility predicate is a compression rule, not the validity rule.
    If it fails, exact support checking may still prove that ``ctx.n -> v`` is
    valid.  In that case the child must refine/fix the assignment instead of
    enlarging the block as a symmetric set-to-set choice.
    """
    _, p_to_block = _sym_block_indexes(cand)
    join_idx = p_to_block.get(v)
    if join_idx is None:
        return None, False
    block = cand.blocks[join_idx]
    if not block.open:
        return None, False
    can_extend = _r_compatible_with_block(
        cand, join_idx, ctx.n, set(ctx.fragment_old),
        ctx.g_R, ctx.r_orbits)
    return join_idx, can_extend


def _collect_free_target_entries(
    cand: _SymCand,
    ctx: _ExtensionContext,
) -> tuple[dict[int, list[TargetEntry]], dict[Any, list[TargetEntry]]]:
    """Find every target atom that can support free ``ctx.n`` extension.

    Results are split into two groups:

    - ``block_join``: target atom lies inside an existing open symmetry block.
      The child either extends that block or refines it under a correlated
      support assignment.
    - ``by_group``: target atom is outside existing blocks.  Entries are
      grouped by element/orbit/context signature so equivalent target atoms can
      become one new `_SymBlock`.
    """
    block_join: dict[int, list[TargetEntry]] = defaultdict(list)
    by_group: dict[Any, list[TargetEntry]] = defaultdict(list)
    for v in ctx.g_P.nodes():
        if v in ctx.locked_p_atoms:
            continue
        if not ctx.node_policy.compatible(ctx.g_R, ctx.n, ctx.g_P, v):
            continue
        join_idx, can_extend = _target_join_info(cand, ctx, v)
        support = _supported_value(cand, ctx, v, join_idx)
        if support is None:
            continue
        if join_idx is not None:
            block_join[join_idx].append((v, support, can_extend))
        else:
            sig = _p_relation_signature(
                cand, v, ctx.g_P, ctx.p_orbits,
                node_policy=ctx.node_policy)
            by_group[sig].append((v, support, True))
    return block_join, by_group


def _children_from_block_join(
    cand: _SymCand,
    n: Node,
    join_idx: int,
    entries: list[TargetEntry],
) -> list[_SymCand]:
    """Build children where ``n`` joins an existing open symmetry block."""
    children: list[_SymCand] = []

    free_entries = [
        (v, support)
        for v, support, can_extend in entries
        if not support and can_extend
    ]
    if free_entries:
        # No old-block witness needs to be fixed.  The existing block expands
        # from, for example, "one H occupies methyl-H pool" to "two H atoms
        # occupy methyl-H pool".
        witness_v = sorted(v for v, _ in free_entries)[0]
        child = cand.with_extended_block(join_idx, n)
        if child is not None:
            child = child.with_witness({n: witness_v})
        if child is not None:
            children.append(child)

    for v, support, can_extend in entries:
        if not support and can_extend:
            continue
        # The join is either correlated with an old block assignment, or exact
        # support exists while the source atom is not block-compatible enough
        # to enlarge the block.  Refine/freeze instead of treating it as an
        # independent symmetric choice.
        fixed = dict(support)
        fixed[n] = v
        child = _refine_sym_assignments(cand, fixed)
        if child is not None:
            children.append(child)
    return children


def _children_from_context_group(
    cand: _SymCand,
    n: Node,
    entries: list[TargetEntry],
) -> list[_SymCand]:
    """Build children for target atoms outside existing symmetry blocks."""
    children: list[_SymCand] = []

    correlated = [(v, support) for v, support, _can_extend in entries
                  if support]
    for v, support in correlated:
        # This target image is outside old blocks, but an old block must still
        # reshuffle to support it.  The child therefore carries a refined
        # correlated witness.
        fixed = dict(support)
        fixed[n] = v
        child = _refine_sym_assignments(cand, fixed)
        if child is not None:
            children.append(child)

    independent = [(v, support) for v, support, _can_extend in entries
                   if not support]
    if not independent:
        return children

    group = tuple(sorted(v for v, _ in independent))
    witness_v, support = sorted(independent, key=lambda item: item[0])[0]
    if len(group) > 1:
        # Symmetric independent choice: one new block represents all injective
        # choices into this target pool.
        child = cand.with_new_block(n, group, extendable=True)
    else:
        child = cand.with_fixed(n, witness_v)
    if child is not None and isinstance(child, _SymCand):
        support = dict(support)
        support[n] = witness_v
        child = child.with_witness(support)
    if child is not None:
        children.append(child)
    return children


def _extend_free_atom(cand: _SymCand, ctx: _ExtensionContext) -> list[_SymCand]:
    """Extend by an unmapped atom using all valid same-element target atoms."""
    block_join, by_group = _collect_free_target_entries(cand, ctx)
    children: list[_SymCand] = []
    for join_idx, entries in sorted(block_join.items()):
        children.extend(_children_from_block_join(cand, ctx.n, join_idx, entries))
    for _, entries in sorted(by_group.items(), key=lambda kv: str(kv[0])):
        children.extend(_children_from_context_group(cand, ctx.n, entries))
    return children


def _extend_one_candidate(cand: _SymCand, ctx: _ExtensionContext) -> list[_SymCand]:
    """Extend one canonical state; return zero or more canonical child states."""
    if not _candidate_covers_fragment(cand, ctx):
        return []
    if ctx.is_merge:
        child = _extend_locked_merge(cand, ctx)
        return [] if child is None else [child]
    return _extend_free_atom(cand, ctx)


def _dedupe_children(children: list[_SymCand], ctx: _ExtensionContext) -> list[_SymCand]:
    """Collapse child states that are equivalent under orbit/boundary context."""
    return _dedup_sym_cands(
        children,
        ctx.g_R,
        ctx.g_P,
        ctx.r_orbits,
        ctx.p_orbits,
        fragment=ctx.sig_fragment,
        deferred_edges=ctx.boundary_edges,
        locked_mapping=ctx.mapping,
        node_policy=ctx.node_policy,
    )
