import sys
from pathlib import Path

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'bench'))
from compare_aam_seed_counts import aggregate


def test_accuracy_denominator_stays_full_but_cpu_uses_common_complete_cases():
    variants={}
    for label,cpu in (('seeds1',2.),('seeds2',3.),('seeds10',8.)):
        variants[label]={('golden',i):dict(searches_complete=not(label=='seeds2' and i==2),
            compute_cpu=cpu,recovery='recovered' if i<2 else 'unknown',best_events=None) for i in range(3)}
        variants[label]['holdout',0]=dict(best_events=4,searches_complete=True,compute_cpu=cpu,recovery='unknown')
    slap={i:dict(incomplete_variants=[0] if i==1 else [],mapping_errors=0,invalid_predictions=0,
                 workflow_cpu_excluding_io=1.,recovered=i==0) for i in range(3)}
    result=aggregate(variants,slap)
    assert result['common_case_indices']==[0]
    assert result['methods']['seeds1']['golden_cases']==3
    assert result['methods']['seeds1']['golden_outcomes']=={'recovered':2,'unknown':1}
    assert result['methods']['seeds1']['common_mean_cpu_seconds']==2.
    assert result['methods']['seeds1']['cpu_reduction_factor_vs_ten']==4.
    assert result['methods']['seeds1']['holdout_best_event_comparison']=={'equal':1}
