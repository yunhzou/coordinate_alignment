"""Differential tests for the compiled symmetry-repair kernel."""
import math
import os
import random
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "bench"))

try:
    from rxn_core import _engine
except ImportError:  # pragma: no cover
    _engine = None

pytestmark = pytest.mark.skipif(
    _engine is None or not hasattr(_engine, "repair_group"),
    reason="native repair kernel not built")


def test_round12_matches_python_round():
    rng = random.Random(1234)
    values = [0.0, 1.0, 0.5, 1e-13, 5e-13, 4.9999999999995e-13, 123456.7890123456789,
              0.1 + 0.2, 1e15 + 0.3, 2.5e-12, 7.5e-12, 1e-300, 1e300]
    for _ in range(200000):
        exponent = rng.uniform(-16, 8)
        mantissa = rng.uniform(0.0, 10.0)
        values.append(mantissa * 10 ** exponent)
        values.append(float(rng.randint(0, 10**6)) * 1e-12 + rng.random() * 1e-13)
    for value in values:
        assert _engine.round12(value) == round(value, 12), value
    for value in (math.inf, -math.inf):
        assert _engine.round12(value) == round(value, 12)
    assert math.isnan(_engine.round12(math.nan))


def test_repair_scores_match_local_scores_on_random_inputs():
    rng = np.random.default_rng(7)
    for _ in range(300):
        n = int(rng.integers(3, 30))
        wbo_P = np.round(rng.random((n, n)) * 3, 3)
        wbo_P = (wbo_P + wbo_P.T) / 2
        np.fill_diagonal(wbo_P, 0.0)
        n_pairs = int(rng.integers(0, 60))
        pair_i = rng.integers(0, n, n_pairs).astype(np.intp)
        pair_j = rng.integers(0, n, n_pairs).astype(np.intp)
        pair_wbo_R = np.round(rng.random(n_pairs) * 3, 3)
        pair_threshold = rng.choice([0.3, 0.5], n_pairs).astype(float)
        pair_r_active = rng.random(n_pairs) > 0.5
        bond_floor = 0.2
        current = rng.permutation(n).astype(np.intp)
        k = int(rng.integers(1, min(6, n) + 1))
        rs = np.sort(rng.choice(n, k, replace=False)).astype(np.intp)
        states = np.array([rng.permutation(current[rs]) for _ in range(20)],
                          dtype=np.intp)

        def local_scores(image_rows):
            pair_wbo_P = wbo_P[image_rows[:, pair_i], image_rows[:, pair_j]]
            difference = pair_wbo_R - pair_wbo_P
            magnitude = np.abs(difference)
            changed_mask = magnitude >= pair_threshold
            contribution = np.where(
                changed_mask, magnitude,
                np.where(pair_r_active | (pair_wbo_P >= bond_floor),
                         magnitude * 0.01, 0.0))
            if contribution.shape[1]:
                totals = np.add.accumulate(contribution, axis=1)[:, -1]
            else:
                totals = np.zeros(contribution.shape[0])
            return [(int(c), round(t, 12)) for c, t in zip(
                np.count_nonzero(changed_mask, axis=1).tolist(), totals.tolist())]

        rows = np.repeat(current[None, :], len(states), axis=0)
        rows[:, rs] = states
        expected = local_scores(rows)
        got = _engine.repair_scores(current, rs, states, pair_i, pair_j, pair_wbo_R,
                                    pair_threshold, pair_r_active, wbo_P, bond_floor)
        assert got == expected


def test_search_aam_agrees_in_verify_mode():
    from cases import CASES
    from rxn_core.aam import search_aam
    from rxn_core.domain import AAMProblem, AAMSearchConfig
    import rxn_core.alignment.branch as branch_mod

    saved = branch_mod._VERIFY_REPAIR
    branch_mod._VERIFY_REPAIR = True
    try:
        for case in ("tempo", "tetraphenyl"):
            R, P = CASES[case]()
            result = search_aam(AAMProblem(R, P, name=case), AAMSearchConfig())
            assert result.mechanisms
    finally:
        branch_mod._VERIFY_REPAIR = saved
