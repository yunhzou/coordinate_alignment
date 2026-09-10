from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'bench'))
from aam_seed_ablation import claim,save,result_folder
from rxn_core.alignment.branch import _generate_seed_orders
from rxn_core.aam import cut_seed
import networkx as nx


class SeedAblationTests(unittest.TestCase):
    def test_reduced_orders_are_the_ten_order_prefix_for_every_cut(self):
        graph=nx.Graph()
        for i,e in enumerate(('C','C','O','H','H','H','H','H','H','Cl')):
            graph.add_node(i,element=e)
        graph.add_edges_from(((0,1),(1,2),(0,3),(0,4),(0,5),(1,6),(1,7),(2,8)))
        for cut in ((),*((edge,) for edge in graph.edges)):
            source=graph.copy();source.remove_edges_from(cut)
            original=_generate_seed_orders(source,10,rng_seed=cut_seed(cut))
            for count in (1,2,3):
                with self.subTest(cut=cut,seeds=count):
                    self.assertEqual(_generate_seed_orders(source,count,rng_seed=cut_seed(cut)),original[:count])

    def test_dynamic_queue_claims_each_task_exactly_once(self):
        with tempfile.TemporaryDirectory() as directory:
            run=Path(directory)
            save(run/'work_order.json',list(reversed(range(73))))
            save(run/'search_queue.json',dict(next=0))
            def consume(_):
                output=[]
                while (slot:=claim(run,'search')) is not None:output.append(slot)
                return output
            values=list(ThreadPoolExecutor(max_workers=8).map(consume,range(8)))
            self.assertEqual(sorted(v for group in values for v in group),list(range(73)))
            self.assertIsNone(claim(run,'search'))

    def test_variant_does_not_overwrite_baseline_paths(self):
        task=dict(dataset='golden',index=1,direction='R_to_P')
        self.assertNotEqual(result_folder(Path('/baseline'),task),result_folder(Path('/variant'),task))


if __name__=='__main__':unittest.main()
