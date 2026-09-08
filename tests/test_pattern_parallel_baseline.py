from itertools import permutations
from benchmark_pattern_collection import baseline_keys
from rxn_core import AAMProblem
from rxn_core.pattern_collection import PatternEquivalence
from test_family_scoring import endpoint


def test_parallel_baseline_has_identical_keys_and_order():
    mol=endpoint('OOO',[(0,1,1)])
    eq=PatternEquivalence(AAMProblem(mol,mol))
    mappings=[dict(enumerate(p)) for p in permutations(range(3))]
    assert list(baseline_keys(eq,mappings,2))==list(baseline_keys(eq,mappings,1))
