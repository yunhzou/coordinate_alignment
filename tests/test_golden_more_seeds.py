from dataclasses import asdict
import networkx as nx

from golden_more_seeds import changed_config
from rxn_core import AAMSearchConfig
from rxn_core.aam import cut_seed
from rxn_core.alignment.branch import _generate_seed_orders


def test_seed_ablation_changes_no_other_parameter():
    old=AAMSearchConfig(seed_count=3);new=changed_config(old,10)
    assert {k for k,v in asdict(old).items() if asdict(new)[k]!=v}=={'seed_count'}


def test_more_seed_orders_preserve_original_prefix():
    graph=nx.path_graph(8)
    nx.set_node_attributes(graph,'C','element')
    for cut in ((),((1,2),)):
        seed=cut_seed(cut)
        first=list(_generate_seed_orders(graph,n_trials=3,rng_seed=seed))
        more=list(_generate_seed_orders(graph,n_trials=10,rng_seed=seed))
        assert more[:3]==first and len(more)==10


def test_high_cap_ablation_changes_only_requested_budgets():
    old=AAMSearchConfig(seed_count=3,branch_limit=100)
    new=changed_config(old,10,2000)
    assert new.seed_count==10 and new.branch_limit==2000
    assert {k for k,v in asdict(old).items() if asdict(new)[k]!=v}=={'seed_count','branch_limit'}
