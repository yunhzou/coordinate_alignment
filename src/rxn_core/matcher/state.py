"""Compressed symmetry candidate state."""
from __future__ import annotations

import os
from collections import defaultdict
from dataclasses import dataclass

from .primitives import _load_fast_kernels, _orbit_id

# ``RXN_CORE_VERIFY_ROLES=1`` recomputes every cached or derived role
# dictionary from scratch and asserts equality (see ``_derive_roles``).
_VERIFY_ROLES = os.environ.get('RXN_CORE_VERIFY_ROLES') == '1'


@dataclass(frozen=True)
class _SymBlock:
    """A local symmetry class inside one partial candidate.

    `r_atoms` are already present in the growing R fragment.  `p_atoms` is
    the target pool they occupy up to local symmetry; the candidate stores a
    deterministic concrete witness, but extension/dedup reason over the pool.
    Initial seed pools are not extendable because they represent "where the
    anchor could be", not a sibling class around an already chosen anchor.
    """
    r_atoms: tuple
    p_atoms: tuple
    extendable: bool = True

    def __post_init__(self):
        object.__setattr__(self, 'r_atoms', tuple(sorted(set(self.r_atoms))))
        object.__setattr__(self, 'p_atoms', tuple(sorted(set(self.p_atoms))))

    @property
    def open(self):
        return len(self.r_atoms) < len(self.p_atoms)

    @property
    def complete(self):
        return len(self.r_atoms) >= len(self.p_atoms)


class _SymCand:
    """Compressed partial mapping with symmetry blocks and one witness.

    The public alignment API still consumes concrete dict mappings.  Inside
    growth, this object keeps target pools grouped so a K-fold symmetric target
    produces one candidate instead of K concrete bijections.  `mapping` is the
    deterministic witness used for cheap WBO existence checks.
    """
    __slots__ = (
        'mapping', 'blocks', 'exact_fixed', 'multiplicity',
        'automorph_blocks', '_roles',
    )

    def __init__(self, mapping=None, blocks=(), exact_fixed=(), multiplicity=1,
                 automorph_blocks=()):
        blocks = tuple(blocks)
        # ``raw`` is only read below and ``m`` is always a fresh dict, so an
        # incoming plain dict needs no defensive copy.
        raw = mapping if type(mapping) is dict else dict(mapping or {})
        if blocks:
            block_r = {r for b in blocks for r in b.r_atoms}
            m = {r: p for r, p in raw.items() if r not in block_r}
            used = set(m.values())
            for b in blocks:
                if len(b.r_atoms) > len(b.p_atoms):
                    raise ValueError("symmetry block has more R atoms than P atoms")
                for r in b.r_atoms:
                    if r not in raw:
                        continue
                    p = raw[r]
                    if p not in b.p_atoms or p in used:
                        raise ValueError("symmetry block witness conflicts with fixed atoms")
                    m[r] = p
                    used.add(p)
                available = [p for p in b.p_atoms if p not in used]
                missing = [r for r in b.r_atoms if r not in m]
                if len(available) < len(missing):
                    raise ValueError("symmetry block witness conflicts with fixed atoms")
                for r, p in zip(missing, available):
                    m[r] = p
                    used.add(p)
        else:
            m = dict(raw)
        self.mapping = m
        self.blocks = blocks
        self.exact_fixed = frozenset(exact_fixed)
        self.multiplicity = int(multiplicity)
        expanded_automorph_blocks = []
        for block in automorph_blocks:
            p_atoms = set(block.p_atoms)
            r_atoms = set(block.r_atoms)
            before = len(r_atoms)
            r_atoms.update(r for r, p in m.items() if p in p_atoms)
            if len(r_atoms) == before and block.extendable is False:
                # Already expanded against this witness: the rebuilt block
                # would compare equal to ``block`` field by field, and its
                # atom tuples are already sorted and unique.
                expanded_automorph_blocks.append(block)
            else:
                expanded_automorph_blocks.append(_SymBlock(
                    tuple(r_atoms), tuple(p_atoms), extendable=False))
        self.automorph_blocks = tuple(expanded_automorph_blocks)
        # Lazily cached ``_cand_roles_from_scratch(self)``; see
        # ``_derive_roles``.  Never part of equality or serialisation.
        self._roles = None

    def __contains__(self, r):
        return r in self.mapping

    def __getitem__(self, r):
        return self.mapping[r]

    def get(self, r, default=None):
        return self.mapping.get(r, default)

    def items(self):
        return self.mapping.items()

    def values(self):
        return self.mapping.values()

    def materialize(self):
        return dict(self.mapping)

    def has_open_choice(self):
        return any(b.open for b in self.blocks)

    def with_fixed(self, r, p):
        if r in self.mapping:
            return self if self.mapping[r] == p else None
        if p in self.mapping.values():
            return None
        m = dict(self.mapping)
        m[r] = p
        try:
            child = _SymCand(m, self.blocks, self.exact_fixed,
                             self.multiplicity, self.automorph_blocks)
        except ValueError:
            return None
        if self._roles is not None:
            # ``r`` was unmapped, hence outside every block, so it becomes a
            # fixed pair; blocks are unchanged.
            child._roles = _derive_roles(
                self, child, added=((p, ('mapped', int(r))),))
        return child

    def with_new_block(self, r, p_atoms, extendable=True):
        p_atoms = tuple(sorted(set(p_atoms)))
        if not p_atoms:
            return None
        if len(p_atoms) == 1:
            return self.with_fixed(r, p_atoms[0])
        if r in self.mapping:
            return None
        b = _SymBlock((r,), p_atoms, extendable=extendable)
        try:
            child = _SymCand(self.mapping, self.blocks + (b,),
                             self.exact_fixed, self.multiplicity,
                             self.automorph_blocks)
        except ValueError:
            return None
        if self._roles is not None:
            # ``r`` was unmapped, so the fixed pairs are unchanged; the new
            # block adds one pool role to each of its target atoms.
            role = _pool_role(b)
            child._roles = _derive_roles(
                self, child, added=[(p, role) for p in b.p_atoms])
        return child

    def with_extended_block(self, block_idx, r):
        if r in self.mapping:
            return None
        blocks = list(self.blocks)
        b = blocks[block_idx]
        if not b.extendable or b.complete:
            return None
        blocks[block_idx] = _SymBlock(b.r_atoms + (r,), b.p_atoms,
                                      extendable=b.extendable)
        try:
            child = _SymCand(self.mapping, blocks, self.exact_fixed,
                             self.multiplicity, self.automorph_blocks)
        except ValueError:
            return None
        if self._roles is not None:
            # Same pool, new member list: swap that block's role on every
            # pool atom.  ``r`` was unmapped, so fixed pairs are unchanged.
            new_b = blocks[block_idx]
            old_role = _pool_role(b)
            new_role = _pool_role(new_b)
            child._roles = _derive_roles(
                self, child,
                removed=[(p, old_role) for p in b.p_atoms],
                added=[(p, new_role) for p in new_b.p_atoms])
        return child

    def with_witness(self, assignments):
        r_to_block = {}
        for idx, block in enumerate(self.blocks):
            for r in block.r_atoms:
                r_to_block[r] = idx
        touched_blocks = {
            r_to_block[r] for r in assignments
            if r in r_to_block
        }
        m = {
            r: p for r, p in self.mapping.items()
            if r not in r_to_block or r_to_block[r] not in touched_blocks
        }
        m.update(assignments)
        try:
            child = _SymCand(m, self.blocks, self.exact_fixed,
                             self.multiplicity, self.automorph_blocks)
        except ValueError:
            return None
        if self._roles is not None:
            # Blocks are unchanged, so only assignments to atoms outside every
            # block can alter the fixed pairs (block witnesses carry no role).
            removed = []
            added = []
            for r, p in assignments.items():
                if r in r_to_block:
                    continue
                old_p = self.mapping.get(r)
                if old_p == p:
                    continue
                role = ('mapped', int(r))
                if old_p is not None:
                    removed.append((old_p, role))
                added.append((p, role))
            child._roles = _derive_roles(
                self, child, removed=removed, added=added)
        return child

    def with_multiplicity(self, multiplicity):
        child = _SymCand(self.mapping, self.blocks, self.exact_fixed,
                         multiplicity, self.automorph_blocks)
        # Roles never depend on the multiplicity; the mapping content,
        # blocks and (already expanded) automorph domains are unchanged.
        child._roles = self._roles
        return child

    def with_automorph_equivalent(self, other):
        """Merge an exactly automorphic witness without storing a bijection."""
        if not isinstance(other, _SymCand):
            other = _SymCand(other)
        varying_r = tuple(sorted(
            r for r in set(self.mapping) & set(other.mapping)
            if self.mapping[r] != other.mapping[r]
        ))
        blocks = list(self.automorph_blocks) + list(other.automorph_blocks)
        if varying_r:
            p_atoms = tuple(sorted({
                p for r in varying_r
                for p in (self.mapping[r], other.mapping[r])
            }))
            blocks.append(_SymBlock(varying_r, p_atoms, extendable=False))

        # The domain is search/display metadata.  The exact automorphism group
        # is reconstructed from the completed AAM fragment relation.
        merged = []
        for block in blocks:
            r_set = set(block.r_atoms)
            p_set = set(block.p_atoms)
            changed = True
            while changed:
                changed = False
                keep = []
                for prior in merged:
                    if (r_set.intersection(prior.r_atoms)
                            or p_set.intersection(prior.p_atoms)):
                        r_set.update(prior.r_atoms)
                        p_set.update(prior.p_atoms)
                        changed = True
                    else:
                        keep.append(prior)
                merged = keep
            merged.append(_SymBlock(
                tuple(r_set), tuple(p_set), extendable=False))
        merged.sort(key=lambda block: (block.r_atoms, block.p_atoms))
        child = _SymCand(
            self.mapping, self.blocks, self.exact_fixed,
            self.multiplicity + other.multiplicity,
            tuple(merged))
        if self._roles is not None:
            # Mapping content and blocks are unchanged; only the automorph
            # domains differ, which ``_derive_roles`` compares itself.
            child._roles = _derive_roles(self, child)
        return child

    def structural_signature(self, g_R, g_P, r_orbits=None, p_orbits=None):
        block_r = {r for b in self.blocks for r in b.r_atoms}
        fixed = tuple(sorted(
            (('exact', r) if r in self.exact_fixed
             else ('orbit', _orbit_id(r_orbits, r)),
             ('exact', p) if r in self.exact_fixed
             else ('orbit', _orbit_id(p_orbits, p)))
            for r, p in self.mapping.items()
            if r not in block_r
        ))
        blocks = tuple(sorted(
            (
                tuple(sorted(_orbit_id(r_orbits, r) for r in b.r_atoms)),
                len(b.r_atoms),
                tuple(sorted(_orbit_id(p_orbits, p) for p in b.p_atoms)),
                len(b.p_atoms),
                bool(b.extendable),
            )
            for b in self.blocks
        ))
        return fixed, blocks


def _pool_role(block):
    return ('pool', tuple(int(r) for r in block.r_atoms),
            bool(block.extendable))


def _automorph_role(block):
    return ('automorph_domain', tuple(int(r) for r in block.r_atoms))


def _cand_roles_from_scratch(cand, group_domains=False):
    """Target atom -> sorted role tuple, computed from the candidate alone.

    This is the reference definition used by
    ``_CandidateAutomorphismCanonicalizer._candidate_roles``: fixed pairs
    contribute ``('mapped', r)`` to their image, every block contributes one
    pool role to each atom of its pool and every automorph domain one group
    role to each of its atoms.  Multi-role atoms sort their roles by ``repr``,
    so the tuple is a function of the role multiset alone.
    """
    if isinstance(cand, _SymCand):
        mapping = cand.mapping
        blocks = cand.blocks
    else:
        mapping = dict(cand)
        blocks = ()

    roles = defaultdict(list)
    block_r = {r for block in blocks for r in block.r_atoms}
    if group_domains and isinstance(cand, _SymCand):
        block_r.update(
            r for block in cand.automorph_blocks
            for r in block.r_atoms)
    for r, p in sorted(mapping.items()):
        if r not in block_r:
            roles[int(p)].append(('mapped', int(r)))
    for block in blocks:
        block_role = _pool_role(block)
        for p in block.p_atoms:
            roles[int(p)].append(block_role)
    if isinstance(cand, _SymCand):
        for block in cand.automorph_blocks:
            group_role = _automorph_role(block)
            for p in block.p_atoms:
                roles[int(p)].append(group_role)
    # A single role needs no ordering; sorted() of one item is that item.
    return {p: (tuple(items) if len(items) == 1
                else tuple(sorted(items, key=repr)))
            for p, items in roles.items()}


def _derive_roles(parent, child, removed=(), added=()):
    """``child``'s role dictionary from ``parent``'s cached one.

    ``removed`` / ``added`` are ``(p, role)`` items describing how the fixed
    mapping and pool roles differ between the two candidates; automorph
    domain roles are diffed here from the two ``automorph_blocks`` tuples.
    Because ``_cand_roles_from_scratch`` sorts each atom's roles by ``repr``
    (an injective key on these tuples), the result depends only on the role
    multiset per atom, so applying the difference to the parent's multisets
    reproduces the from-scratch dictionary exactly.  Returns ``None`` when
    the parent has no cached roles.
    """
    base = parent._roles
    if base is None:
        return None
    roles = dict(base)
    work = {}

    def items_for(p):
        p = int(p)
        items = work.get(p)
        if items is None:
            items = work[p] = list(roles.get(p, ()))
        return items

    try:
        for p, role in removed:
            items_for(p).remove(role)
        for p, role in added:
            items_for(p).append(role)
        old_blocks = parent.automorph_blocks
        new_blocks = child.automorph_blocks
        if (len(old_blocks) == len(new_blocks)
                and all(old.p_atoms == new.p_atoms
                        for old, new in zip(old_blocks, new_blocks))):
            for old, new in zip(old_blocks, new_blocks):
                if old.r_atoms != new.r_atoms:
                    old_role = _automorph_role(old)
                    new_role = _automorph_role(new)
                    for p in old.p_atoms:
                        items = items_for(p)
                        items.remove(old_role)
                        items.append(new_role)
        else:
            for old in old_blocks:
                role = _automorph_role(old)
                for p in old.p_atoms:
                    items_for(p).remove(role)
            for new in new_blocks:
                role = _automorph_role(new)
                for p in new.p_atoms:
                    items_for(p).append(role)
    except ValueError:
        # The parent's roles do not contain an item the delta removes; fall
        # back to the from-scratch computation on demand.
        if _VERIFY_ROLES:
            raise
        return None
    for p, items in work.items():
        if not items:
            roles.pop(p, None)
        elif len(items) == 1:
            roles[p] = tuple(items)
        else:
            roles[p] = tuple(sorted(items, key=repr))
    if _VERIFY_ROLES:
        expected = _cand_roles_from_scratch(child)
        assert roles == expected, (
            "derived candidate roles differ from the from-scratch roles",
            roles, expected)
    return roles


def _cand_map(cand):
    return cand.materialize() if isinstance(cand, _SymCand) else dict(cand)


def _cand_possible_p_atoms(cand):
    atoms = set(_cand_map(cand).values())
    if isinstance(cand, _SymCand):
        for block in cand.blocks:
            atoms.update(block.p_atoms)
    return atoms


def _cand_has_open_choice(cand):
    return isinstance(cand, _SymCand) and cand.has_open_choice()


def _sym_block_assignment_expr(block):
    n = len(block.p_atoms)
    k = len(block.r_atoms)
    if k <= 0 or n <= 1:
        return '1'
    if k == 1:
        return str(n)
    if k == n:
        return f'{n}!'
    return f'P({n},{k})'


def _symmetry_state(cand, r_orbits=None, p_orbits=None):
    """Serialize only symmetry carried by the candidate itself.

    The orbit parameters are accepted for call-site compatibility, but display
    symmetry is intentionally not inferred from global endpoint orbit tables.
    `_SymCand`/island state is the automorphism unit for the matching
    algorithm.
    """
    item = {
        'witness': {int(a): int(b) for a, b in _cand_map(cand).items()},
        'blocks': [],
    }
    if isinstance(cand, _SymCand):
        item['exact_fixed'] = [int(x) for x in sorted(cand.exact_fixed)]
        item['multiplicity'] = int(cand.multiplicity)
        item['automorph_blocks'] = [{
            'r_atoms': [int(x) for x in block.r_atoms],
            'p_atoms': [int(x) for x in block.p_atoms],
            'extendable': False,
            'open': False,
            'assignments': 'exact_group',
            'source': 'exact_automorph_group',
        } for block in cand.automorph_blocks]
        for block in cand.blocks:
            item['blocks'].append({
                'r_atoms': [int(x) for x in block.r_atoms],
                'p_atoms': [int(x) for x in block.p_atoms],
                'extendable': bool(block.extendable),
                'open': bool(block.open),
                'assignments': _sym_block_assignment_expr(block),
            })
        item['blocks'].extend(item['automorph_blocks'])
    return item


def _sym_block_indexes(cand):
    r_to_block = {}
    p_to_block = {}
    if isinstance(cand, _SymCand):
        for idx, block in enumerate(cand.blocks):
            for r in block.r_atoms:
                r_to_block[r] = idx
            for p in block.p_atoms:
                p_to_block[p] = idx
    return r_to_block, p_to_block


# The pure-Python constructor stays reachable (the differential test compares
# against it); ``__init__`` is rebound only when RXN_CORE_FAST=1 selects the
# compiled extension.  This runs last so the extension, which imports names
# from this module, sees a fully initialised module.
_SymCand.__init___py = _SymCand.__init__
_fast = _load_fast_kernels()
if _fast is not None:
    # The compiled constructor predates the roles cache slot and the
    # block-free/automorph fast paths of the Python __init__; keep Python.
    pass
