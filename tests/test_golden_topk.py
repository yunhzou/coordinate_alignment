import importlib.util
from pathlib import Path

spec=importlib.util.spec_from_file_location('golden_topk',Path(__file__).parents[1]/'bench/golden_topk.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)


def test_event_windows_do_not_trade_coverage_for_fewer_events():
    rows=[dict(rank=1,score=[-10,-20,4],reference_equivalent=False),
          dict(rank=2,score=[-10,-20,5],reference_equivalent=True),
          dict(rank=3,score=[-9,-20,0],reference_equivalent=False)]
    result=module.summarize_classes(rows)
    assert result['reference_event_gap']==1
    assert result['event_windows']['0']==dict(count=1,recovered=False)
    assert result['event_windows']['1']==dict(count=2,recovered=True)


def test_topfive_ties_include_equal_score_without_reference_guidance():
    rows=[dict(rank=i+1,score=[-10,-20,4 if i<7 else 5],reference_equivalent=i==6) for i in range(8)]
    result=module.summarize_classes(rows)
    assert result['reference_class_rank']==7
    assert result['top5_tie_inclusive']==dict(count=7,recovered=True)
