"""Differential oracles for native bookkeeping, independent of chemistry cases."""
from itertools import permutations
import random

import pytest

from rxn_core._group_ops import (occupation_orbit, project_generators,
                                 conjugate_generators, OccupationLimitExceeded)


def reference(witness, degree, stages, attachments, fragments, bonds, limit=-1):
    def key(images):
        return (tuple(sorted(images)), tuple(sorted(images[i] for i in attachments)),
                tuple(sorted((label, tuple(sorted(images[i] for i in positions)))
                             for label, positions in fragments)),
                tuple(sorted(tuple(sorted((images[a], images[b]))) for a, b in bonds)))
    states = {key(witness): (tuple(witness), tuple(range(degree)))}
    for stage in stages:
        queue = list(states.values())
        for images, action in queue:
            for generator in stage:
                moved = tuple(generator[p] for p in images)
                relation = key(moved)
                if relation in states:
                    continue
                if limit >= 0 and len(states) >= limit:
                    raise OccupationLimitExceeded()
                state = moved, tuple(generator[p] for p in action)
                states[relation] = state
                queue.append(state)
    return list(states.values())


def test_exact_orbit_matches_reference_including_first_witness_and_stage_order():
    rng = random.Random(416)
    for degree in range(1, 8):
        for _ in range(40):
            size = rng.randrange(1, degree + 1)
            witness = rng.sample(range(degree), size)
            stages = []
            for _ in range(3):
                stages.append([rng.sample(range(degree), degree) for _ in range(rng.randrange(3))])
            fragments = [(i % 2, list(range(i, size, 2))) for i in range(min(2, size))]
            attachments = list(range(0, size, 3))
            bonds = [(i, i + 1) for i in range(size - 1)]
            actual = occupation_orbit(witness, degree, stages, attachments, fragments, bonds, -1)
            assert [(tuple(images), tuple(action)) for images, action in actual] == reference(
                witness, degree, stages, attachments, fragments, bonds)


def test_limits_and_invalid_native_inputs_are_reported():
    with pytest.raises(OccupationLimitExceeded):
        occupation_orbit([0], 2, [[[1, 0]]], [], [(0, [0])], [], 1)
    with pytest.raises(ValueError):
        occupation_orbit([0], 2, [[[1, 1]]], [], [(0, [0])], [], -1)
    with pytest.raises(ValueError):
        occupation_orbit([0], 2, [], [], [(0, [1])], [], -1)
    with pytest.raises(ValueError):
        occupation_orbit([0], 2, [], [], [(0, [0])], [], -1, [2])


def _observed_relation(images, observed, attachments, fragments, bonds):
    values = tuple(p if p in observed else -1 for p in images)
    return (tuple(sorted(values)), tuple(sorted(values[i] for i in attachments)),
            tuple(sorted((label, tuple(sorted(values[i] for i in positions)))
                         for label, positions in fragments)),
            tuple(sorted(tuple(sorted((values[a], values[b]))) for a, b in bonds)))


def test_permanently_competitor_only_choices_stay_compressed():
    witness = [0, 2, 3]
    stages = [[[0, 1, 3, 2, 4, 5], [0, 1, 2, 4, 3, 5], [0, 1, 2, 3, 5, 4]]]
    fragments = [(0, [0]), (1, [1]), (1, [2])]
    full = occupation_orbit(witness, 6, stages, [], fragments, [], -1)
    compact = occupation_orbit(witness, 6, stages, [], fragments, [], -1, [0, 1])
    assert len(full) == 6
    assert len(compact) == 1
    assert compact[0] == full[0]


def test_future_actions_keep_current_competitors_distinct():
    # The first stage exchanges two competitors; the second sends them to
    # different target atoms. Projecting to the target prematurely loses atom 3.
    stages = [[[0, 2, 1, 3]], [[1, 0, 3, 2]]]
    compact = occupation_orbit([1], 4, stages, [], [(0, [0])], [], -1, [0, 3])
    full = occupation_orbit([1], 4, stages, [], [(0, [0])], [], -1)
    assert compact == full
    assert {images[0] for images, _ in compact} == {0, 1, 2, 3}


def test_observation_quotient_preserves_exhaustive_projected_relations():
    rng = random.Random(1976)
    for degree in range(2, 9):
        for _ in range(50):
            # Some actions mix target/competitor positions; others leave an
            # independent competitor component that can really be compressed.
            active = rng.randrange(1, degree + 1)
            stages = []
            for _ in range(3):
                generators = []
                for _ in range(rng.randrange(4)):
                    generators.append(rng.sample(range(active), active) +
                                      rng.sample(range(active, degree), degree - active))
                stages.append(generators)
            size = rng.randrange(1, degree + 1)
            witness = rng.sample(range(degree), size)
            observed = set(rng.sample(range(active), rng.randrange(active + 1)))
            fragments = [(i % 2, list(range(i, size, 2))) for i in range(min(2, size))]
            attachments = list(range(0, size, 3))
            bonds = [(i, i + 1) for i in range(size - 1)]
            full = reference(witness, degree, stages, attachments, fragments, bonds)
            compact = occupation_orbit(witness, degree, stages, attachments, fragments,
                                       bonds, -1, sorted(observed))
            key = lambda images: _observed_relation(images, observed, attachments, fragments, bonds)
            assert {key(images) for images, _ in compact} == {key(images) for images, _ in full}
            for images, action in compact:
                assert sorted(action) == list(range(degree))
                assert list(images) == [action[p] for p in witness]


def test_projection_preserves_generator_order_and_sparse_atom_labels():
    raw = [(1, 0, 3, 2), (0, 1, 2, 3), (1, 0, 3, 2)]
    assert project_generators(raw, [0, 3], [0, 1]) == ((3, 1, 2, 0),)
    with pytest.raises(ValueError):
        project_generators([(2, 1, 0)], [0, 1], [0, 1])


def test_conjugation_preserves_padded_frames_exhaustively():
    for action in permutations(range(4)):
        for generator in permutations(range(3)):
            expected = list(range(4))
            for i in range(4):
                expected[action[i]] = action[generator[i] if i < 3 else i]
            assert conjugate_generators([generator], action) == (tuple(expected),)
