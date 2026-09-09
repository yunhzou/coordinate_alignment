from pathlib import Path
import sys
import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'bench'))
from compare_elementary_outputs import event_counts, features, certificate
from rxn_core import AAMProblem
from rxn_core.domain import MolecularEndpoint
from rxn_core.family_scoring import bond_events
import json
from types import SimpleNamespace
from compare_elementary_outputs import refine


def test_vector_scores_match_existing_explicit_h_scorer():
    rng=np.random.default_rng(42)
    r=rng.choice([0.,.15,.8,1.,1.7],(6,6));r=np.triu(r,1);r+=r.T
    p=rng.choice([0.,.15,.8,1.,1.7],(6,6));p=np.triu(p,1);p+=p.T
    problem=AAMProblem(MolecularEndpoint(('C',)*6,np.zeros((6,3)),r),
                       MolecularEndpoint(('C',)*6,np.zeros((6,3)),p))
    maps=[rng.permutation(6) for _ in range(20)]
    for m,counts in zip(maps,event_counts(r,p,maps)):
        expected=bond_events(problem,dict(enumerate(map(int,m))))
        assert tuple(counts)==tuple(expected[k] for k in ('broken','formed','order_changed'))


def test_symmetry_equivalence_not_raw_index_equality():
    w=[[0,1,0],[1,0,1],[0,1,0]]
    raw={side:dict(elements=['C','C','C'],wbo=w) for side in ('reactant','product')}
    f=features(raw)
    assert certificate(f,[0,1,2],[0,1,2])==certificate(f,[2,1,0],[0,1,2])
    assert certificate(f,[0,1,2],[0,1,2])!=certificate(f,[1,0,2],[0,1,2])


def test_native_h_refinement_without_changing_heavy_mapping(tmp_path):
    inputs=tmp_path/'inputs/0';inputs.mkdir(parents=True)
    native=tmp_path/'slap_xyz';native.mkdir()
    endpoints={k:dict(elements=['C','H','H'],coordinates=np.zeros((3,3)).tolist(),
        wbo=w) for k,w in [('reactant',[[0,.7,1.4],[.7,0,0],[1.4,0,0]]),
                          ('product',[[0,1.4,.7],[1.4,0,0],[.7,0,0]])]}
    (inputs/'input.json').write_text(json.dumps(endpoints))
    (native/'0.json').write_text(json.dumps(dict(candidates=[dict(graphs=[dict(labels=[1,2,2])]*2)])))
    (tmp_path/'0.json').write_text(json.dumps(dict(scope='',slap=[dict(
        events=dict(broken=0,formed=0,order_changed=2,total=2),mapping=[[0,0],[1,1],[2,2]],
        all_h_label_permutations_score_invariant=False)])))
    refine(SimpleNamespace(run=tmp_path,source=tmp_path,index=0))
    result=json.loads((tmp_path/'0.refined.json').read_text())['slap'][0]
    assert result['events']['total']==0
    assert result['hydrogen_score_optimization']['optimal']
    assert result['mapping']==[[0,0],[1,2],[2,1]]
