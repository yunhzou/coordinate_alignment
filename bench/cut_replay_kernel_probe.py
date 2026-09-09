"""Same-input default-native regression against a saved pre-edit binary."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import time
import numpy as np

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--binary',type=Path,required=True)
p.add_argument('--output',type=Path,required=True)
a=p.parse_args()
spec=importlib.util.spec_from_file_location('rxn_core._engine',a.binary)
module=importlib.util.module_from_spec(spec)
sys.modules['rxn_core._engine']=module
spec.loader.exec_module(module)

from rxn_core.frag import build_graph
from rxn_core.matcher import _nauty_orbits
from rxn_core.growth import native

source=Path('/project/yunhengzou/coordinate_alignment/aam_benchmarks/elementary140_tol1_20260909/inputs')
rows=[]
for index in (0,4,25,28,76,114,135):
    raw=json.loads((source/str(index)/'input.json').read_text())
    r,p=[build_graph(raw[k]['elements'],np.asarray(raw[k]['wbo']),bond_cut=.2) for k in ('reactant','product')]
    rv=native.source_graph(r);pv=native.target_graph(p,_nauty_orbits(p,wbo_tol=1.))
    calls=[]
    for seed in (0,len(r)//2,len(r)-1):
        args=(rv.graph,pv.graph,seed,[-1]*len(r),.2,1.,1,100,None,(),False)
        cold=time.process_time();out=module.grow_island(*args);cold=time.process_time()-cold
        start=time.process_time()
        for _ in range(7):actual=module.grow_island(*args)
        cpu=time.process_time()-start
        assert actual==out
        calls.append(dict(seed=seed,cold_cpu=cold,warm_seven_cpu=cpu,
            digest=hashlib.sha256(json.dumps(out,sort_keys=True).encode()).hexdigest()))
    rows.append(dict(index=index,calls=calls))
a.output.write_text(json.dumps(dict(binary=str(a.binary),cases=rows),indent=2)+'\n')
print(sum(c['cold_cpu']+c['warm_seven_cpu'] for r in rows for c in r['calls']))
