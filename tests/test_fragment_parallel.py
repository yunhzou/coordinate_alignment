"""Execution buffering must not become either a search cap or an output queue."""
import pytest

from rxn_core.fragment_matching import parallel


def test_augmentation_stream_is_ordered_and_bounds_outstanding_results(monkeypatch):
    submitted, consumed, closed = [], [], []

    class Result:
        def __init__(self, placement):
            self.placement = placement

        def get(self):
            consumed.append(self.placement)
            return self.placement

    class Pool:
        def __init__(self, workers, **kwargs):
            assert workers == 2

        def __enter__(self):
            return self

        def __exit__(self, *args):
            closed.append(True)

        def apply_async(self, function, args):
            submitted.append(args[0])
            assert len(submitted) - len(consumed) <= 2
            return Result(args[0])

    monkeypatch.setattr(parallel.mp, 'get_context', lambda _: type('Context', (), {'Pool': Pool}))
    output = parallel._parallel_augmentation_results(None, None, None, None,
                                                    tuple(range(7)), 2)
    assert not submitted
    assert next(output) == 0
    assert submitted == [0, 1]
    assert list(output) == list(range(1, 7))
    assert submitted == consumed == list(range(7))
    assert closed == [True]


@pytest.mark.parametrize('placements', [(), ('only',)])
def test_small_augmentation_stream_keeps_all_results(monkeypatch, placements):
    monkeypatch.setattr(parallel, '_augment_initial_family', lambda s, t, p, c, r: p)
    assert tuple(parallel._parallel_augmentation_results(
        None, None, None, None, placements, 4)) == placements


def test_merge_freezes_once_and_preserves_every_derivation(monkeypatch):
    from types import SimpleNamespace
    import rxn_core.fragment_matching.detection as detection
    from rxn_core.fragment_matching.models import FragmentCandidate, FragmentDerivation, FragmentDetectionConfig
    placement = SimpleNamespace(retained_atoms=(0,), representative_mapping=((0, 0),), encounter_count=1)
    initial = ((placement,), 0, 0, False, False, 1, 0, False, ())
    a, b = FragmentDerivation((), (), ((0, 0),)), FragmentDerivation((), (), ((0, 1), (1, 0)))
    raw = FragmentCandidate('raw', ((0, 0),), (0,), (0,), (), (), (), (), (), 2,
                            retained_fragments=((0,),), derivations=(a,))
    other = detection.replace(raw, derivations=(b,))
    replacements = []
    original = detection.replace
    def counted(*args, **kwargs):
        replacements.append(args[0])
        return original(*args, **kwargs)
    monkeypatch.setattr(detection, 'replace', counted)
    result = detection._detect_fragments_from_initial(None, SimpleNamespace(graph=None), initial,
        source_id='final', config=FragmentDetectionConfig(), region=None,
        augmentation_runner=lambda _: (((raw, other) * 100, False, 1, ()),))
    assert len(replacements) == len(result.candidates) == 1
    assert result.candidates[0].source_id == 'final'
    assert result.candidates[0].derivations == (a, b) * 100
