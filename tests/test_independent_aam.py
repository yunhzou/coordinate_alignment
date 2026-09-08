import numpy as np
from rxn_core.domain import AAMProblem, MolecularEndpoint
from rxn_core.independent_aam import detect_independent, independent_paths, assemble_independent


def problem(elements, edges, target_edges):
    def endpoint(bonds):
        w = np.zeros((len(elements), len(elements)))
        for a,b,value in bonds:
            w[a,b] = w[b,a] = value
        return MolecularEndpoint(elements, np.zeros((len(elements),3)), w)
    return AAMProblem(endpoint(edges), endpoint(target_edges))


def test_all_requested_seeds_and_online_dedup():
    p = problem(('C','O','H'), [(0,1,1),(1,2,1)], [(0,1,1),(1,2,1)])
    d = detect_independent(p, seeds=(0,1))
    assert [x['seed'] for x in d.attempts] == [0,1]
    paths = tuple(independent_paths([d.graph]))
    assert len(paths) == 1
    answer = assemble_independent(p, paths)
    assert answer['status'] == 'covered'
    assert len(answer['mapping']) == 3


def test_no_invented_hydrogen_completion():
    p = problem(('C','H'), [], [])
    d = detect_independent(p, seeds=(0,))
    answer = assemble_independent(p, independent_paths([d.graph]))
    assert answer['status'] == 'not_covered'
    assert answer['missing_source_atoms'] == [1]


def test_compressed_symmetry_preserves_correlations():
    p = problem(('C',)*4, [(0,1,1),(1,2,1),(2,3,1)], [(0,1,1),(1,2,1),(2,3,1)])
    d = detect_independent(p, seeds=(0,1,2,3))
    paths = tuple(independent_paths([d.graph]))
    assert assemble_independent(p, paths, required_mapping={0:3,1:2,2:1,3:0})['status'] == 'covered'
    assert assemble_independent(p, paths, required_mapping={0:0,1:2,2:1,3:3})['status'] == 'not_covered'


def test_disjoint_fragments_assemble_and_do_not_share_targets():
    p = problem(('C','O','C','O'), [(0,1,1),(2,3,1)], [(0,1,1),(2,3,1)])
    d = detect_independent(p, seeds=(0,1,2,3))
    result = assemble_independent(p, independent_paths([d.graph]))
    assert result['status'] == 'covered'
    assert len(result['selected']) == 2
    assert len(set(dict(result['mapping']).values())) == 4


def test_unrelated_cuts_do_not_duplicate_identical_finalized_relations():
    p = problem(('C','O','N','F'), [(0,1,1),(2,3,1)], [(0,1,1),(2,3,1)])
    base = detect_independent(p, seeds=(0,))
    cut = detect_independent(p, seeds=(0,), cuts=((2,3),))
    assert len(tuple(independent_paths([base.graph, cut.graph]))) == 1


def test_compatible_overlap_is_not_rejected_as_a_partition_conflict():
    p = problem(('C','N','O'), [(0,1,1),(1,2,1)], [(0,1,1),(1,2,1)])
    ab = detect_independent(p, seeds=(0,), cuts=((1,2),))
    bc = detect_independent(p, seeds=(2,), cuts=((0,1),))
    result = assemble_independent(p, independent_paths([ab.graph, bc.graph]))
    assert result['status'] == 'covered'
    assert dict(result['mapping']) == {0:0,1:1,2:2}
