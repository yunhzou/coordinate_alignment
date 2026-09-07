import importlib.util
from pathlib import Path

import numpy as np
import pytest

from rxn_core import AAMProblem
from rxn_core.domain import MolecularEndpoint
from golden_evaluation import rank_key

spec=importlib.util.spec_from_file_location('golden_fast_rerank',Path(__file__).parents[1]/'bench/golden_fast_rerank.py')
M=importlib.util.module_from_spec(spec);spec.loader.exec_module(M)


def endpoint(elements,bonds):
    matrix=np.zeros((len(elements),len(elements)))
    for a,b,w in bonds:matrix[a,b]=matrix[b,a]=w
    return MolecularEndpoint(elements,np.zeros((len(elements),3)),matrix)


def test_unchanged_mapping_has_no_cost():
    e=endpoint(('C','H'),[(0,1,1.)]);f=M.edit_features(AAMProblem(e,e),{0:0,1:1})
    assert f['events']==f['center_components']==f['center_atoms']==f['weighted_edits']==0


def test_explicit_H_edits_and_partial_mapping_match_baseline():
    r=endpoint(('C','O','H','Cl'),[(0,2,1.),(0,3,1.)])
    p=endpoint(('C','O','H','Br'),[(1,2,1.),(0,3,1.)]);problem=AAMProblem(r,p)
    mapping={0:0,1:1,2:2};f=M.edit_features(problem,mapping)
    assert f['events']==rank_key(mapping,problem)[0][2]==4
    assert f['center_components']==1
    assert f['weighted_edits']==pytest.approx((411+327+459+285)/346)


def test_unknown_energy_is_explicit_and_aromatic_convention_is_fixed():
    with pytest.raises(KeyError):M.bond_energy(('C','Hg'),1.)
    assert M.bond_energy(('C','C'),1.5)==474.
    r=endpoint(('C','Hg'),[(0,1,1.)]);p=endpoint(('C','Hg'),[])
    f=M.edit_features(AAMProblem(r,p),{0:0,1:1})
    assert f['weighted_edits'] is None and f['missing_energies']


def test_ranking_features_do_not_read_reference_labels():
    row=dict(rank=2,score=[-10,-20,4],reference_equivalent=True,
        features=dict(events=4,center_components=2,center_atoms=6,weighted_edits=3.5))
    for variant in M.VARIANTS:
        key=M.candidate_key(row,variant)
        assert key==M.candidate_key(dict(row,reference_equivalent=False),variant)
