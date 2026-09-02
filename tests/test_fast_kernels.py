"""Differential tests for the optional compiled matcher kernels.

Every kernel in ``rxn_core.matcher._fast`` must return exactly what its
pure-Python original returns: same values, same container structure, same
element types (``bool`` is not ``int`` here) and same ordering.  Two sources
of inputs are used: calls recorded from the Python originals while a real
symmetric growth runs, and random candidate/graph/orbit-map inputs that also
exercise the fallback paths (plain orbit dicts, no structural zero bucket,
graphs without a WBO matrix, the support-state cap).

The tests are skipped when the extension has not been built
(``.venv/bin/python bench/build_fast.py``).
"""
from __future__ import annotations

import functools
import os
import random
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

_fast = pytest.importorskip(
    "rxn_core.matcher._fast",
    reason="compiled kernels not built; run bench/build_fast.py")

from rxn_core import build_graph  # noqa: E402
from rxn_core.growth import grow_island  # noqa: E402
from rxn_core.matcher import _SymBlock, _SymCand, _nauty_orbits  # noqa: E402
from rxn_core.matcher import canonical as canonical_mod  # noqa: E402
from rxn_core.matcher import dedupe as dedupe_mod  # noqa: E402
from rxn_core.matcher import extend as extend_mod  # noqa: E402
from rxn_core.matcher import support as support_mod  # noqa: E402
from rxn_core.matcher.canonical import _CandidateAutomorphismCanonicalizer  # noqa: E402
from rxn_core.matcher.dedupe import _BoundaryContext  # noqa: E402
from rxn_core.matcher.orbits import _OrbitMap  # noqa: E402
from rxn_core.matcher.primitives import _wbo_bucket  # noqa: E402
from rxn_core.matcher.state import _cand_map, _cand_possible_p_atoms, _sym_block_indexes  # noqa: E402

PY_SYMCAND_INIT = _SymCand.__init___py
PY_SIGNATURE = dedupe_mod._p_relation_signature_from_parts_py
PY_POOL = dedupe_mod._pool_target_signatures_py
PY_WITNESS = support_mod._support_witness_for_value_py
PY_ROLES = _CandidateAutomorphismCanonicalizer._candidate_roles_py
PY_ROLE_KEY = _CandidateAutomorphismCanonicalizer.role_key_from_roles_py
PY_COLORS = _CandidateAutomorphismCanonicalizer._colored_vertices_from_roles_py


# --------------------------------------------------------------------------
# strict structural equality: values, container shapes, element types, order
# --------------------------------------------------------------------------

def _same(a, b):
    if type(a) is not type(b):
        return False
    if isinstance(a, (tuple, list)):
        return len(a) == len(b) and all(_same(x, y) for x, y in zip(a, b))
    if isinstance(a, dict):
        if list(a.keys()) != list(b.keys()):
            return False
        if not all(type(x) is type(y) for x, y in zip(a.keys(), b.keys())):
            return False
        return all(_same(a[k], b[k]) for k in a)
    if isinstance(a, (frozenset, set)):
        return a == b and sorted(map(repr, a)) == sorted(map(repr, b))
    return a == b and repr(a) == repr(b)


def _assert_same(native, python, what):
    assert _same(native, python), f"{what}: native {native!r} != python {python!r}"


# --------------------------------------------------------------------------
# graphs and candidates
# --------------------------------------------------------------------------

def _tetraphenylmethane():
    elements = ["C"]
    bonds = []
    rings = []
    for _ring in range(4):
        base = len(elements)
        carbons = list(range(base, base + 6))
        elements.extend(["C"] * 6)
        for k in range(6):
            bonds.append((carbons[k], carbons[(k + 1) % 6], 1.4))
        bonds.append((0, carbons[0], 1.0))
        hydrogens = []
        for k in range(1, 6):
            h = len(elements)
            elements.append("H")
            bonds.append((carbons[k], h, 1.0))
            hydrogens.append(h)
        rings.append((carbons, hydrogens))
    return elements, bonds, rings


def _matrix(n, bonds):
    w = np.zeros((n, n))
    for a, b, v in bonds:
        w[a, b] = w[b, a] = v
    return w


def _symmetric_case():
    """Tetraphenylmethane and a 1,2-H-shift product: four equivalent rings."""
    elements, bonds, rings = _tetraphenylmethane()
    n = len(elements)
    carbons, hydrogens = rings[0]
    moved = hydrogens[1]
    product_bonds = [b for b in bonds if not (b[0] == carbons[2] and b[1] == moved)]
    product_bonds.append((carbons[3], moved, 1.0))
    g_R = build_graph(elements, _matrix(n, bonds), bond_cut=0.2)
    g_P = build_graph(elements, _matrix(n, product_bonds), bond_cut=0.2)
    return g_R, g_P


def _random_graph(rng, n, elements=("C", "H", "O")):
    els = [rng.choice(elements) for _ in range(n)]
    w = np.zeros((n, n))
    for a in range(n):
        for b in range(a + 1, n):
            roll = rng.random()
            if roll < 0.25:
                value = rng.choice([0.9, 1.0, 1.05, 1.1, 1.5, 2.0, 2.1])
            elif roll < 0.35:
                value = rng.uniform(0.05, 0.45)   # sub-floor contact
            else:
                value = 0.0
            w[a, b] = w[b, a] = value
    return build_graph(els, w, bond_cut=0.5)


def _random_symcand(rng, g_P, n_atoms_R, max_blocks=2):
    """A valid _SymCand with a random witness and up to ``max_blocks`` pools."""
    p_nodes = list(g_P.nodes())
    for _attempt in range(50):
        rng.shuffle(p_nodes)
        r_atoms = list(range(n_atoms_R))
        rng.shuffle(r_atoms)
        k = rng.randint(1, min(len(p_nodes), n_atoms_R))
        chosen_r = r_atoms[:k]
        chosen_p = p_nodes[:k]
        mapping = dict(zip(chosen_r, chosen_p))
        blocks = []
        free_p = [p for p in p_nodes if p not in chosen_p]
        block_r_pool = list(chosen_r)
        for _b in range(rng.randint(0, max_blocks)):
            if len(block_r_pool) < 1 or len(free_p) < 1:
                break
            n_r = rng.randint(1, min(2, len(block_r_pool)))
            br = [block_r_pool.pop() for _ in range(n_r)]
            pool = [mapping[r] for r in br]
            extra = rng.randint(0, min(2, len(free_p)))
            for _ in range(extra):
                pool.append(free_p.pop())
            blocks.append(_SymBlock(tuple(br), tuple(pool),
                                    extendable=rng.random() < 0.7))
        try:
            return _SymCand(mapping, tuple(blocks),
                            exact_fixed=tuple(rng.sample(chosen_r, rng.randint(0, 1))))
        except ValueError:
            continue
    return _SymCand({0: p_nodes[0]})


def _orbit_variants(rng, g_P):
    """Exact orbit map, plain dict, None, partial buckets, no zero bucket."""
    exact = _nauty_orbits(g_P, wbo_tol=0.2)
    partial = _OrbitMap(dict(exact), wbo_buckets={
        key: bucket for key, bucket in exact.wbo_buckets.items()
        if rng.random() < 0.5}, zero_bucket=exact.zero_bucket,
        wbo_tol=exact.wbo_tol)
    no_zero = _OrbitMap(dict(exact), wbo_buckets={
        key: bucket for key, bucket in exact.wbo_buckets.items()
        if rng.random() < 0.5}, zero_bucket=None, wbo_tol=exact.wbo_tol)
    return [exact, dict(exact), None, partial, no_zero]


# --------------------------------------------------------------------------
# recorded calls from a real growth
# --------------------------------------------------------------------------

class _Recorder:
    """Wraps a Python original, recording ``(args, kwargs, result)``.

    Installed either as a module attribute or as a class attribute; in the
    latter case ``__get__`` binds the instance like a plain function would.
    """

    def __init__(self, function):
        self.function = function
        self.calls = []

    def __call__(self, *args, **kwargs):
        result = self.function(*args, **kwargs)
        self.calls.append((args, kwargs, result))
        return result

    def __get__(self, instance, owner=None):
        if instance is None:
            return self
        return functools.partial(self, instance)


def _symcand_state(cand):
    return (cand.mapping, cand.blocks, cand.exact_fixed, cand.multiplicity,
            cand.automorph_blocks)


class _SymCandInitRecorder:
    """Records ``_SymCand.__init__`` arguments and the resulting state, or the
    ValueError the Python constructor raised."""

    def __init__(self):
        self.calls = []

    def __get__(self, instance, owner=None):
        if instance is None:
            return self
        return functools.partial(self, instance)

    def __call__(self, instance, *args, **kwargs):
        try:
            PY_SYMCAND_INIT(instance, *args, **kwargs)
        except ValueError as exc:
            self.calls.append((args, kwargs, ("ValueError", str(exc))))
            raise
        self.calls.append((args, kwargs, _symcand_state(instance)))


@pytest.fixture(scope="module")
def recorded():
    """Inputs and outputs of the Python originals during symmetric growths."""
    g_R, g_P = _symmetric_case()
    p_orbits = _nauty_orbits(g_P, wbo_tol=0.2)
    r_orbits = _nauty_orbits(g_R, wbo_tol=0.2)
    rec = {
        "signature": _Recorder(PY_SIGNATURE),
        "pool": _Recorder(PY_POOL),
        "witness": _Recorder(PY_WITNESS),
        "roles": _Recorder(PY_ROLES),
        "role_key": _Recorder(PY_ROLE_KEY),
        "colors": _Recorder(PY_COLORS),
        "init": _SymCandInitRecorder(),
    }
    saved = (
        extend_mod._p_relation_signature_from_parts,
        dedupe_mod._p_relation_signature_from_parts,
        dedupe_mod._pool_target_signatures,
        extend_mod._support_witness_for_value,
        _CandidateAutomorphismCanonicalizer._candidate_roles,
        _CandidateAutomorphismCanonicalizer.role_key_from_roles,
        _CandidateAutomorphismCanonicalizer._colored_vertices_from_roles,
        _SymCand.__init__,
    )
    extend_mod._p_relation_signature_from_parts = rec["signature"]
    dedupe_mod._p_relation_signature_from_parts = rec["signature"]
    dedupe_mod._pool_target_signatures = rec["pool"]
    extend_mod._support_witness_for_value = rec["witness"]
    _CandidateAutomorphismCanonicalizer._candidate_roles = rec["roles"]
    _CandidateAutomorphismCanonicalizer.role_key_from_roles = rec["role_key"]
    _CandidateAutomorphismCanonicalizer._colored_vertices_from_roles = rec["colors"]
    _SymCand.__init__ = rec["init"]
    # These recordings exercise the Python kernels, so the compiled growth
    # engine (which bypasses them entirely) is disabled for the growth calls.
    native_setting = os.environ.get("RXN_CORE_NATIVE")
    os.environ["RXN_CORE_NATIVE"] = "0"
    try:
        # A free growth from the central carbon and from a ring hydrogen, plus
        # one with a locked prefix so locked roles and merges are exercised.
        for seed, mapping in ((0, {}), (26, {}), (2, {0: 0, 1: 1})):
            grow_island(g_R, g_P, seed, dict(mapping), graph_floor=0.2,
                        iso_tol=0.5, p_orbits=p_orbits, r_orbits=r_orbits)
    finally:
        if native_setting is None:
            os.environ.pop("RXN_CORE_NATIVE", None)
        else:
            os.environ["RXN_CORE_NATIVE"] = native_setting
        (extend_mod._p_relation_signature_from_parts,
         dedupe_mod._p_relation_signature_from_parts,
         dedupe_mod._pool_target_signatures,
         extend_mod._support_witness_for_value,
         _CandidateAutomorphismCanonicalizer._candidate_roles,
         _CandidateAutomorphismCanonicalizer.role_key_from_roles,
         _CandidateAutomorphismCanonicalizer._colored_vertices_from_roles,
         _SymCand.__init__,
         ) = saved
    for name, recorder in rec.items():
        assert recorder.calls, f"no {name} calls were recorded"
    return rec


def test_recorded_relation_signatures_match(recorded):
    calls = recorded["signature"].calls
    assert any(args[0].blocks for args, _k, _r in calls), "no block-carrying call"
    for args, kwargs, expected in calls:
        _assert_same(_fast.p_relation_signature_from_parts(*args, **kwargs),
                     expected, "p_relation_signature_from_parts")


def test_recorded_pool_signatures_match(recorded):
    for args, kwargs, expected in recorded["pool"].calls:
        _assert_same(_fast.pool_target_signatures(*args, **kwargs),
                     expected, "pool_target_signatures")


def test_recorded_support_witnesses_match(recorded):
    calls = recorded["witness"].calls
    assert any(r for _a, _k, r in calls), "no correlated witness was recorded"
    assert any(r is None for _a, _k, r in calls), "no rejection was recorded"
    for args, kwargs, expected in calls:
        _assert_same(_fast.support_witness_for_value(*args, **kwargs),
                     expected, "support_witness_for_value")


def test_recorded_candidate_roles_match(recorded):
    # The Python implementation now serves a per-candidate roles cache whose
    # dictionaries are derived incrementally, so their insertion order can
    # differ from a from-scratch computation.  Every consumer iterates the
    # atom index or sorts/hashes the items, so only dict equality matters.
    # The compiled kernel is not wired in for this function (see
    # canonical.py) but must still agree item for item.
    for args, kwargs, expected in recorded["roles"].calls:
        native = _fast.candidate_roles(*args, **kwargs)
        assert native == expected, f"candidate_roles: {native!r} != {expected!r}"
        assert all(_same(native[k], expected[k]) for k in native)


def test_recorded_role_keys_match(recorded):
    calls = recorded["role_key"].calls
    assert any(not r[1] for _a, _k, r in calls), "no non-singleton key recorded"
    for args, kwargs, expected in calls:
        _assert_same(_fast.role_key_from_roles(*args, **kwargs), expected,
                     "role_key_from_roles")


def test_recorded_colourings_match(recorded):
    for args, kwargs, expected in recorded["colors"].calls:
        _assert_same(_fast.colored_vertices_from_roles(*args, **kwargs),
                     expected, "colored_vertices_from_roles")


def _native_symcand_state(args, kwargs):
    cand = _SymCand.__new__(_SymCand)
    try:
        _fast.symcand_init(cand, *args, **kwargs)
    except ValueError as exc:
        return ("ValueError", str(exc))
    return _symcand_state(cand)


def test_recorded_symcand_constructions_match(recorded):
    calls = recorded["init"].calls
    assert any(kw.get("automorph_blocks") or (len(a) > 4 and a[4])
               for a, kw, _r in calls), "no automorph-block construction recorded"
    for args, kwargs, expected in calls:
        _assert_same(_native_symcand_state(args, kwargs), expected,
                     "_SymCand.__init__")


def test_symcand_init_rejections_match():
    """Each ValueError path of the constructor fires identically."""
    cases = [
        # more R atoms than P atoms in a block
        (({0: 5}, (_SymBlock((1, 2), (6,)),)), {}),
        # witness image outside the block pool
        (({0: 5, 1: 9}, (_SymBlock((1,), (6, 7)),)), {}),
        # witness image already used by a fixed atom
        (({0: 5, 1: 5}, (_SymBlock((1,), (5, 6)),)), {}),
        # not enough free pool atoms for the missing block members
        (({0: 6, 1: 7}, (_SymBlock((2, 3), (6, 7, 8)),)), {}),
        # a valid construction that completes two missing witnesses
        (({0: 5}, (_SymBlock((1, 2), (6, 7, 8)),)), dict(multiplicity=3)),
    ]
    saw_error = False
    for args, kwargs in cases:
        cand = _SymCand.__new__(_SymCand)
        try:
            PY_SYMCAND_INIT(cand, *args, **kwargs)
            expected = _symcand_state(cand)
        except ValueError as exc:
            expected = ("ValueError", str(exc))
            saw_error = True
        _assert_same(_native_symcand_state(args, kwargs), expected,
                     f"_SymCand.__init__ {args}")
    assert saw_error


@pytest.mark.parametrize("seed", range(20))
def test_random_symcand_constructions(seed):
    rng = random.Random(1300 + seed)
    n_p = rng.randint(4, 9)
    n_r = rng.randint(3, 8)
    for _trial in range(40):
        mapping = {}
        for r in rng.sample(range(n_r), rng.randint(0, n_r)):
            mapping[r] = rng.randrange(n_p)          # may repeat images
        blocks = []
        for _b in range(rng.randint(0, 3)):
            r_atoms = tuple(rng.sample(range(n_r), rng.randint(1, 3)))
            p_atoms = tuple(rng.sample(range(n_p), rng.randint(1, 4)))
            blocks.append(_SymBlock(r_atoms, p_atoms, extendable=rng.random() < 0.5))
        automorph = []
        for _b in range(rng.randint(0, 2)):
            automorph.append(_SymBlock(
                tuple(rng.sample(range(n_r), rng.randint(1, 2))),
                tuple(rng.sample(range(n_p), rng.randint(1, 3))),
                extendable=False))
        args = (dict(mapping), tuple(blocks))
        kwargs = dict(exact_fixed=tuple(rng.sample(range(n_r), rng.randint(0, 2))),
                      multiplicity=rng.randint(1, 4),
                      automorph_blocks=tuple(automorph))
        if rng.random() < 0.2:
            args = (None,)
            kwargs = {}
        cand = _SymCand.__new__(_SymCand)
        try:
            PY_SYMCAND_INIT(cand, *args, **kwargs)
            expected = _symcand_state(cand)
        except ValueError as exc:
            expected = ("ValueError", str(exc))
        _assert_same(_native_symcand_state(args, kwargs), expected,
                     "_SymCand.__init__ random")


# --------------------------------------------------------------------------
# random inputs, including the fallback paths
# --------------------------------------------------------------------------

@pytest.mark.parametrize("seed", range(12))
def test_random_relation_signatures(seed):
    rng = random.Random(seed)
    g_P = _random_graph(rng, rng.randint(5, 11))
    cand = _random_symcand(rng, g_P, n_atoms_R=rng.randint(3, 9))
    cm_items = tuple(sorted(_cand_map(cand).items()))
    for orbits in _orbit_variants(rng, g_P):
        for v in g_P.nodes():
            for compact in (False, True):
                for parts in (dict(cm_items=cm_items, blocks=cand.blocks), {}):
                    kwargs = dict(parts, compact=compact)
                    _assert_same(
                        _fast.p_relation_signature_from_parts(
                            cand, v, g_P, orbits, **kwargs),
                        PY_SIGNATURE(cand, v, g_P, orbits, **kwargs),
                        f"signature orbits={type(orbits).__name__} v={v}")
    # a plain dict candidate and a graph without a WBO matrix
    plain = dict(_cand_map(cand))
    g_edges = g_P.copy()
    del g_edges.graph["wbo_matrix"]
    g_edges.graph.pop("_wbo_rows", None)
    for graph in (g_P, g_edges):
        orbits = _nauty_orbits(graph, wbo_tol=0.2)
        for v in graph.nodes():
            for compact in (False, True):
                _assert_same(
                    _fast.p_relation_signature_from_parts(
                        plain, v, graph, orbits, compact=compact),
                    PY_SIGNATURE(plain, v, graph, orbits, compact=compact),
                    "signature plain mapping")


def test_wbo_bucket_fallback_rounding():
    """Signatures on orbit maps without a zero bucket use _wbo_bucket(w)."""
    rng = random.Random(3)
    n = 6
    els = ["C"] * n
    w = np.zeros((n, n))
    for a in range(n):
        for b in range(a + 1, n):
            # values on and near the round-half-even boundaries of x*5
            w[a, b] = w[b, a] = rng.choice([0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 0.0, 0.25])
    g_P = build_graph(els, w, bond_cut=0.05)
    orbits = _OrbitMap({v: v for v in g_P.nodes()}, wbo_buckets={},
                       zero_bucket=None, wbo_tol=0.2)
    cand = _SymCand({0: 0, 1: 1, 2: 2}, (_SymBlock((1, 2), (1, 2, 3)),))
    for v in g_P.nodes():
        native = _fast.p_relation_signature_from_parts(cand, v, g_P, orbits)
        _assert_same(native, PY_SIGNATURE(cand, v, g_P, orbits), "bucket fallback")
        for r, bucket in native[2]:
            assert bucket == _wbo_bucket(w[cand.mapping[r], v])


@pytest.mark.parametrize("seed", range(40))
def test_random_support_witnesses(seed):
    rng = random.Random(100 + seed)
    g_P = _random_graph(rng, rng.randint(5, 12))
    n_R = rng.randint(4, 10)
    cand = _random_symcand(rng, g_P, n_atoms_R=n_R, max_blocks=3)
    mapped = sorted(cand.mapping)
    n = rng.choice([r for r in range(n_R + 1) if r not in cand.mapping] or [n_R])
    block_indexes = _sym_block_indexes(cand)
    for _trial in range(6):
        bonded = sorted(rng.sample(mapped, rng.randint(1, min(4, len(mapped)))))
        if rng.random() < 0.15:
            bonded.append(n_R + 5)   # an atom the candidate does not map
        r_wbos = [(u, rng.choice([0.9, 1.0, 1.1, 2.0, 0.3])) for u in bonded]
        iso_tol = rng.choice([0.2, 0.5, 1.0])
        strict = None
        if rng.random() < 0.5 and bonded:
            strict = {bonded[0]: rng.choice([0.9, 1.0, 2.0])}
        max_states = rng.choice([4096, 1, 2, 3, 5, 4096])
        for v_n in g_P.nodes():
            join = block_indexes[1].get(v_n)
            choices = [None]
            if join is not None:
                choices.append(join)
            if cand.blocks:
                choices.append(rng.randrange(len(cand.blocks)))
            for join_idx in choices:
                for bi in (None, block_indexes):
                    kwargs = dict(join_block_idx=join_idx, strict_r_wbos=strict,
                                  max_states=max_states, block_indexes=bi)
                    try:
                        expected = PY_WITNESS(cand, n, v_n, list(bonded),
                                              list(r_wbos), g_P, iso_tol, **kwargs)
                    except KeyError:
                        with pytest.raises(KeyError):
                            _fast.support_witness_for_value(
                                cand, n, v_n, list(bonded), list(r_wbos),
                                g_P, iso_tol, **kwargs)
                        continue
                    _assert_same(
                        _fast.support_witness_for_value(
                            cand, n, v_n, list(bonded), list(r_wbos), g_P,
                            iso_tol, **kwargs),
                        expected, f"witness v_n={v_n} join={join_idx}")
        # concrete dict candidates take the fixed-image path
        plain = dict(cand.mapping)
        for v_n in g_P.nodes():
            _assert_same(
                _fast.support_witness_for_value(
                    plain, n, v_n, list(bonded), list(r_wbos), g_P, iso_tol,
                    strict_r_wbos=strict),
                PY_WITNESS(plain, n, v_n, list(bonded), list(r_wbos), g_P,
                           iso_tol, strict_r_wbos=strict),
                "witness plain mapping")


def test_support_witness_cap_semantics_and_non_int_cap():
    """The state cap counts backtracking calls identically; a non-int cap is
    delegated to the Python original."""
    n = 6
    els = ["C", "H", "H", "H", "H", "O"]
    w = np.zeros((n, n))
    for h in (1, 2, 3, 4):
        w[0, h] = w[h, 0] = 1.0
    w[1, 5] = w[5, 1] = 0.9
    w[2, 5] = w[5, 2] = 0.9
    g_P = build_graph(els, w, bond_cut=0.5)
    cand = _SymCand({0: 0, 10: 1, 11: 2}, (_SymBlock((10, 11), (1, 2, 3, 4)),))
    bonded = [10, 11]
    r_wbos = [(10, 1.0), (11, 1.0)]
    for cap in (1, 2, 3, 4, 5, 6, 4096):
        _assert_same(
            _fast.support_witness_for_value(cand, 12, 5, bonded, r_wbos, g_P, 0.3,
                                            max_states=cap),
            PY_WITNESS(cand, 12, 5, bonded, r_wbos, g_P, 0.3, max_states=cap),
            f"cap={cap}")
    for cap in (float("inf"), 2.0, 10 ** 30):
        _assert_same(
            _fast.support_witness_for_value(cand, 12, 5, bonded, r_wbos, g_P, 0.3,
                                            max_states=cap),
            PY_WITNESS(cand, 12, 5, bonded, r_wbos, g_P, 0.3, max_states=cap),
            f"cap={cap}")


@pytest.mark.parametrize("seed", range(12))
def test_random_role_kernels(seed):
    rng = random.Random(500 + seed)
    g_P = _random_graph(rng, rng.randint(5, 11))
    orbits = _nauty_orbits(g_P, wbo_tol=0.2)
    nodes = list(g_P.nodes())
    locked = {}
    if rng.random() < 0.6:
        k = rng.randint(1, 2)
        for r, p in zip(range(50, 50 + k), rng.sample(nodes, k)):
            locked[r] = p
    canon = _CandidateAutomorphismCanonicalizer(
        g_P, p_orbits=orbits, locked_mapping=locked)
    for _trial in range(4):
        cand = _random_symcand(rng, g_P, n_atoms_R=rng.randint(3, 8), max_blocks=3)
        if rng.random() < 0.5:
            other = cand.with_witness({})
            cand = cand.with_automorph_equivalent(other)
        for group in (False, True):
            roles_py = PY_ROLES(canon, cand, group_domains=group)
            _assert_same(_fast.candidate_roles(canon, cand, group_domains=group),
                         roles_py, "candidate_roles")
        roles_py = PY_ROLES(canon, cand)
        _assert_same(_fast.role_key_from_roles(canon, roles_py, orbits),
                     PY_ROLE_KEY(canon, roles_py, orbits), "role_key")
        _assert_same(_fast.colored_vertices_from_roles(canon, roles_py),
                     PY_COLORS(canon, roles_py), "colored_vertices")
        plain = dict(cand.mapping)
        _assert_same(_fast.candidate_roles(canon, plain), PY_ROLES(canon, plain),
                     "candidate_roles plain")


@pytest.mark.parametrize("seed", range(10))
def test_random_pool_signatures(seed):
    rng = random.Random(900 + seed)
    g_R = _random_graph(rng, rng.randint(6, 10))
    g_P = _random_graph(rng, rng.randint(6, 11))
    r_orbits = _nauty_orbits(g_R, wbo_tol=0.2)
    for p_orbits in (_nauty_orbits(g_P, wbo_tol=0.2), dict(_nauty_orbits(g_P, wbo_tol=0.2))):
        cand = _random_symcand(rng, g_P, n_atoms_R=g_R.number_of_nodes(), max_blocks=2)
        fragment = set(cand.mapping) | {rng.randrange(g_R.number_of_nodes())}
        outside = [r for r in g_R.nodes() if r not in fragment]
        deferred = [(rng.choice(sorted(fragment)), x) for x in outside[:3]]
        if not deferred:
            continue
        context = _BoundaryContext(g_R, g_P, fragment, deferred, r_orbits,
                                   p_orbits, {}, None)
        cm_items = tuple(sorted(_cand_map(cand).items()))
        inverse = {p: r for r, p in cm_items}
        mapped_r = tuple(r for r, _p in cm_items)
        used_possible = _cand_possible_p_atoms(cand)

        def general(v, cand=cand, p_orbits=p_orbits, cm_items=cm_items):
            return PY_SIGNATURE(cand, v, g_P, p_orbits, cm_items=cm_items,
                                blocks=cand.blocks, compact=True)

        for pool_key in context.pools:
            args = (context.pools[pool_key], used_possible, inverse, mapped_r,
                    cand.blocks, context.target_static, general)
            _assert_same(_fast.pool_target_signatures(*args), PY_POOL(*args),
                         "pool_target_signatures")


# --------------------------------------------------------------------------
# wiring: the extension is selected only by RXN_CORE_FAST=1
# --------------------------------------------------------------------------

_PROBE = (
    "import sys; sys.path.insert(0, %r);"
    "from rxn_core.matcher import support, dedupe, canonical, extend;"
    "C = canonical._CandidateAutomorphismCanonicalizer;"
    "from rxn_core.matcher.state import _SymCand;"
    "print(int(support._support_witness_for_value is support._support_witness_for_value_py),"
    " int(extend._support_witness_for_value is support._support_witness_for_value),"
    " int(dedupe._p_relation_signature_from_parts is dedupe._p_relation_signature_from_parts_py),"
    " int(extend._p_relation_signature_from_parts is dedupe._p_relation_signature_from_parts),"
    " int(dedupe._pool_target_signatures is dedupe._pool_target_signatures_py),"
    " int(C._candidate_roles is C._candidate_roles_py),"
    " int(C.role_key_from_roles is C.role_key_from_roles_py),"
    " int(C._colored_vertices_from_roles is C._colored_vertices_from_roles_py),"
    " int(_SymCand.__init__ is _SymCand.__init___py))"
) % str(ROOT / "src")


def _probe(env_value):
    env = dict(os.environ)
    env.pop("RXN_CORE_FAST", None)
    if env_value is not None:
        env["RXN_CORE_FAST"] = env_value
    out = subprocess.run([sys.executable, "-c", _PROBE], env=env,
                         capture_output=True, text=True, check=True)
    return [int(x) for x in out.stdout.split()]


def test_python_path_is_default_and_fast_path_is_opt_in():
    assert _probe(None) == [1, 1, 1, 1, 1, 1, 1, 1, 1]
    assert _probe("0") == [1, 1, 1, 1, 1, 1, 1, 1, 1]
    # with the variable set the leaf kernels are rebound to the extension and
    # extend's imported names follow the selected implementation; the roles
    # computation and the _SymCand constructor deliberately stay in Python
    # because they carry the per-candidate roles cache (canonical.py,
    # state.py).
    assert _probe("1") == [0, 1, 0, 1, 0, 1, 0, 0, 1]
