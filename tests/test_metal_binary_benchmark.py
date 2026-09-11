import copy
from pathlib import Path
import sys
import unittest

import numpy as np
import z3

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'bench'))
from metal_binary_events import binary_metal_input,DeltaPatterns,delta_objective,query_pattern,scalar_events
from rxn_core import AAMProblem,AAMSearchConfig,search_aam
from rxn_core.domain import MolecularEndpoint
from rxn_core.family_query import compile_path


def raw_problem(elements,r_edges,p_edges):
    def endpoint(edges):
        w=np.zeros((len(elements),len(elements)))
        for a,b,v in edges:w[a,b]=w[b,a]=v
        return dict(elements=list(elements),coordinates=np.zeros((len(elements),3)).tolist(),wbo=w.tolist())
    return dict(name='unit',reactant=endpoint(r_edges),product=endpoint(p_edges))


class MetalBinaryTests(unittest.TestCase):
    def test_copy_preserves_nonmetal_weights_and_graph_threshold(self):
        raw=raw_problem(['V','O','C','H'],[(0,1,.19),(0,2,.2),(1,2,1.8),(2,3,1.01)],[])
        before=copy.deepcopy(raw);value=binary_metal_input(raw)
        self.assertEqual(raw,before)
        self.assertEqual(value['reactant']['wbo'][0][1],0)
        self.assertEqual(value['reactant']['wbo'][0][2],1)
        self.assertEqual(value['reactant']['wbo'][1][2],1.8)
        self.assertEqual(value['reactant']['wbo'][2][3],1.01)

    def test_raw_delta_thresholds_are_pair_specific_and_inclusive(self):
        raw=raw_problem(['V','O','C','H'],[(0,1,1.),(2,3,1.),(1,2,.1)],[(0,1,.65),(2,3,.65),(1,2,.6)])
        events=DeltaPatterns(raw).describe([0,1,2,3])
        self.assertEqual(events['counts'],dict(broken=1,formed=1))
        self.assertEqual(events['events']['broken'],((0,1),))
        self.assertEqual(events['events']['formed'],((1,2),))
        self.assertEqual(scalar_events(raw,dict(enumerate([0,1,2,3])))['total'],2)

    def test_equivalent_hydrogen_event_patterns_are_canonicalized(self):
        raw=raw_problem(['C','H','H'],[(0,1,1),(0,2,1)],[(0,1,1),(0,2,.1)])
        patterns=DeltaPatterns(raw)
        a,b=patterns.describe([0,1,2]),patterns.describe([0,2,1])
        self.assertNotEqual(a['events'],b['events'])
        self.assertEqual(a['id'],b['id'])

    def test_family_constraints_use_binary_graph_but_score_raw_wbo(self):
        raw=raw_problem(['V','O','O'],[(0,1,2.2),(0,2,.7)],[(0,1,.8),(0,2,2.1)])
        search=binary_metal_input(raw)
        problem=AAMProblem(*(MolecularEndpoint(**search[s]) for s in ('reactant','product')))
        result=search_aam(problem,AAMSearchConfig(seed_count=1,branch_limit=32,cut_floor=10,iso_tolerance=1.),
                          execution='reused_native',workers=1)
        canonical=DeltaPatterns(raw)
        # Tiny exhaustive oracle: the only two element-preserving mappings.
        targets=[canonical.describe(v) for v in ([0,1,2],[0,2,1])]
        self.assertEqual({p['total'] for p in targets},{0,2})
        found=set()
        for path in result.graph.paths():
            if len(path.mapping)!=3:continue
            compiled=compile_path(path,problem,{},source_atoms=(),complete_reference=False)
            objective,lower=delta_objective(compiled,canonical)
            for pattern in targets:
                status,witness=query_pattern(compiled,canonical,pattern,objective,5000)
                if status=='represented':
                    found.add(pattern['id']);self.assertEqual(witness['total'],pattern['total'])
            compiled.solver.push();compiled.solver.add(objective==1)
            self.assertEqual(compiled.solver.check(),z3.unsat);compiled.solver.pop()
        self.assertEqual(found,{p['id'] for p in targets})


if __name__=='__main__':unittest.main()
