"""Differential tests for the native ``_engine.freeze_analytical``.

``rxn_core._engine.freeze_analytical`` (native/src/freeze.cpp) must return
exactly what the pure-Python ``_freeze_analytical_py`` in
``rxn_core.alignment.sweep`` returns: same values, same container structure,
same element types at every level (``bool`` is not ``int``, ``tuple`` is not
``list``), same ordering, and the very same object for non-container leaves.

Skipped when the extension has not been built
(``.venv/bin/python native/build_engine.py``).
"""
from __future__ import annotations

import collections
import os

import pytest

from rxn_core.alignment import sweep
from rxn_core.alignment.sweep import _freeze_analytical_py

try:
    from rxn_core import _engine
except ImportError:  # pragma: no cover - depends on the build
    _engine = None

freeze_native = getattr(_engine, "freeze_analytical", None)

pytestmark = pytest.mark.skipif(
    freeze_native is None, reason="native freeze_analytical is not built")

try:
    import numpy as np
except ImportError:  # pragma: no cover
    np = None


def _assert_same(a, b):
    """Equal values and identical types, recursively through tuples."""
    assert type(a) is type(b), (type(a), type(b))
    if isinstance(a, tuple):
        assert len(a) == len(b)
        for x, y in zip(a, b):
            _assert_same(x, y)
    else:
        assert a == b


class _Leaf:
    """Opaque non-container object; must come back by identity."""


def _symmetry_state(n_blocks=6):
    """Realistic fragment 'symmetry' state dict from the sweep."""
    return {
        'witness': {i: i for i in range(40)},
        'blocks': [{'r_atoms': [1, 2, 3], 'p_atoms': [4, 5, 6],
                    'extendable': True, 'open': False,
                    'assignments': '3!'}] * n_blocks,
        'automorph_blocks': [],
        'exact_fixed': list(range(20)),
        'multiplicity': n_blocks,
    }


LEAF = _Leaf()

CASES = [
    # leaves
    None, True, False, 0, 1, -7, 3.5, 0.0, 'x', '', b'bytes', LEAF,
    # empty containers
    {}, [], (), set(), frozenset(),
    collections.OrderedDict(),
    # lists vs tuples (structure and element types must survive)
    [1, 2, 3], (1, 2, 3), [(1, 2), [3, 4]], ([1, 2], (3, 4)),
    [True, 1, 1.0, '1'], (None, False, 0),
    # dicts with int / str / tuple / bool keys, including 1 and '1' together
    {1: 'a', '1': 'b'},
    {'1': 'b', 1: 'a'},
    {1: {'x': 1}, '1': {'x': 2}},
    {(1, 2): 'a', (1,): 'b', 12: 'c', '(1, 2)': 'd'},
    {'b': [1, {2, (3, 4)}], 'a': ('x',)},
    {'z': 1, 'a': 2, 'm': 3},
    {10: 'ten', 9: 'nine', 100: 'hundred'},  # str keys sort as '10' < '100' < '9'
    {True: 'T', 2: 'two', None: 'none', 2.5: 'f'},
    # nested dicts
    {'outer': {'inner': {'deep': [1, {'k': (2, 3)}]}}, 'list': [{'a': 1}, {'b': 2}]},
    # sets / frozensets of ints, strings and tuples (repr ordering)
    {3, 1, 2}, frozenset({3, 1, 2}),
    {'b', 'a', 'c'}, frozenset({'b', 'a'}),
    {(1, 2), (0, 9), (1, 1)}, frozenset({(2,), (1, 2), (1,)}),
    {1, '1', (1,), 1.5, None, True},          # mixed types: only repr is comparable
    {frozenset({1, 2}), frozenset({3}), frozenset()},
    {10, 9, 100},                             # repr order: '10' < '100' < '9'
    # list / tuple containing sets
    [{2, 1}, frozenset({'b', 'a'})], ({3, 2, 1},),
    # dict subclass
    collections.OrderedDict([('b', 1), ('a', 2), ('c', [3, {4}])]),
    # realistic sweep state
    _symmetry_state(),
    {'hierarchy': _symmetry_state(2), 'cuts': [(1, 2), (3, 4)],
     'mapping': {i: i + 10 for i in range(8)}},
]

if np is not None:
    CASES += [
        np.int64(3),
        [np.int64(1), np.int32(2), 3],
        {np.int64(1): np.int64(2), 'k': [np.float64(0.5)]},
        {np.int64(1), np.int64(2)},
    ]


@pytest.mark.parametrize("value", CASES, ids=lambda v: repr(v)[:60])
def test_native_matches_python(value):
    _assert_same(freeze_native(value), _freeze_analytical_py(value))


def test_leaf_returned_by_identity():
    for leaf in (LEAF, 'string', 12345678901234567890, 3.25, None, b'x'):
        assert freeze_native(leaf) is leaf
    if np is not None:
        v = np.int64(7)
        assert freeze_native(v) is v


def test_containers_are_new_tuples():
    src = [1, [2, 3]]
    out = freeze_native(src)
    assert type(out) is tuple and type(out[1]) is tuple
    assert out == (1, (2, 3))
    src_t = (1, (2, 3))
    assert type(freeze_native(src_t)) is tuple


def test_ordered_dict_after_move_to_end_matches_items_order():
    # OrderedDict.move_to_end reorders .items() but not the underlying dict,
    # so the port has to use .items() for dict subclasses; the stable sort
    # keeps input order between equal-comparing pairs, which is observable
    # only through the element types here (1 vs 1.0 compare equal).
    od = collections.OrderedDict([(1, 1), ('1', 1.0)])
    od.move_to_end(1)
    _assert_same(freeze_native(od), _freeze_analytical_py(od))
    od2 = collections.OrderedDict([('1', 1.0), (1, 1)])
    od2.move_to_end('1')
    _assert_same(freeze_native(od2), _freeze_analytical_py(od2))


def test_stable_sort_keeps_input_order_of_equal_pairs():
    for d in ({1: 1, '1': 1.0}, {'1': 1.0, 1: 1}, {1: True, '1': 1}):
        _assert_same(freeze_native(d), _freeze_analytical_py(d))


def test_unorderable_pairs_raise_same_exception():
    bad = {1: 1, '1': 'a'}          # ('1', 1) vs ('1', 'a') -> int < str
    with pytest.raises(TypeError):
        _freeze_analytical_py(bad)
    with pytest.raises(TypeError):
        freeze_native(bad)
    bad_set_free = {1, 'a'}         # sets are fine: sorted by repr
    _assert_same(freeze_native(bad_set_free), _freeze_analytical_py(bad_set_free))


def test_exceptions_from_str_and_repr_propagate():
    class BadStr:
        def __hash__(self):
            return 1

        def __str__(self):
            raise ValueError("no str")

    class BadRepr:
        def __repr__(self):
            raise KeyError("no repr")

    with pytest.raises(ValueError, match="no str"):
        freeze_native({BadStr(): 1})
    with pytest.raises(KeyError):
        freeze_native({BadRepr(), BadRepr()})


def test_deep_nesting_raises_recursion_error_not_crash():
    value = []
    for _ in range(100000):
        value = [value]
    with pytest.raises(RecursionError):
        _freeze_analytical_py(value)
    with pytest.raises(RecursionError):
        freeze_native(value)


def test_sweep_binds_native_unless_disabled():
    if os.environ.get("RXN_CORE_NATIVE", "1") != "0":
        assert sweep._freeze_analytical is freeze_native
    else:
        assert sweep._freeze_analytical is _freeze_analytical_py
    assert sweep._freeze_analytical_py is _freeze_analytical_py
