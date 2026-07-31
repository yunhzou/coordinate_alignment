"""Compressed symmetry candidate state."""
from __future__ import annotations

from dataclasses import dataclass

from .primitives import _orbit_id


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
        'automorph_blocks',
    )

    def __init__(self, mapping=None, blocks=(), exact_fixed=(), multiplicity=1,
                 automorph_blocks=()):
        blocks = tuple(blocks)
        block_r = {r for b in blocks for r in b.r_atoms}
        raw = dict(mapping or {})
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
        self.mapping = m
        self.blocks = blocks
        self.exact_fixed = frozenset(exact_fixed)
        self.multiplicity = int(multiplicity)
        expanded_automorph_blocks = []
        for block in automorph_blocks:
            p_atoms = set(block.p_atoms)
            r_atoms = set(block.r_atoms)
            r_atoms.update(r for r, p in m.items() if p in p_atoms)
            expanded_automorph_blocks.append(_SymBlock(
                tuple(r_atoms), tuple(p_atoms), extendable=False))
        self.automorph_blocks = tuple(expanded_automorph_blocks)

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
            return _SymCand(m, self.blocks, self.exact_fixed,
                            self.multiplicity, self.automorph_blocks)
        except ValueError:
            return None

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
            return _SymCand(self.mapping, self.blocks + (b,),
                            self.exact_fixed, self.multiplicity,
                            self.automorph_blocks)
        except ValueError:
            return None

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
            return _SymCand(self.mapping, blocks, self.exact_fixed,
                            self.multiplicity, self.automorph_blocks)
        except ValueError:
            return None

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
            return _SymCand(m, self.blocks, self.exact_fixed,
                            self.multiplicity, self.automorph_blocks)
        except ValueError:
            return None

    def with_multiplicity(self, multiplicity):
        return _SymCand(self.mapping, self.blocks, self.exact_fixed,
                        multiplicity, self.automorph_blocks)

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
        return _SymCand(
            self.mapping, self.blocks, self.exact_fixed,
            self.multiplicity + other.multiplicity,
            tuple(merged))

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
