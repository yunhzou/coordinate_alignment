"""Compare minimum-score bond-edit patterns in the saved 140-case outputs.

Events include explicit H and are jointly canonicalized on the reactant graph;
equal event counts or equal heavy-atom mappings are not treated as equal edits.
This is offline analysis of frozen outputs, never a new mapping search.
"""
import argparse
from collections import Counter, defaultdict
from functools import lru_cache
import hashlib
import json
import os
from pathlib import Path
import shlex
import shutil
import socket
import subprocess
import time

import numpy as np
import pynauty
from golden_evaluation import colored_graph
from compare_elementary_outputs import event_counts
from rxn_core import AAMProblem
from rxn_core.domain import MolecularEndpoint
from rxn_core.artifacts import read_aam_checkpoint

DATA = Path('/project/yunhengzou/coordinate_alignment/aam_benchmarks')
AAM = DATA/'holdout_cap1000_seed1_20260910'
SWEEP = DATA/'holdout_slap_xyz_sweep_20260910'
PYTHON = '/project/yunhengzou/coordinate_alignment/.venv/bin/python'
KINDS = ('broken','formed','strengthened','weakened')


def read(path):
    return json.loads(path.read_text())


def save(path,value):
    path.parent.mkdir(parents=True,exist_ok=True)
    tmp = path.with_suffix(path.suffix+'.tmp')
    tmp.write_text(json.dumps(value,indent=2)+'\n')
    tmp.replace(path)


def sha(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream,'sha256').hexdigest()


def signed_features(raw):
    result = []
    for side,other in [('reactant','product'),('product','reactant')]:
        w = np.asarray(raw[side]['wbo'])
        v = np.asarray(raw[other]['wbo'])
        values = np.unique(v[np.triu(v,1)>.2])
        bonds = []
        for a,b in zip(*np.where(np.triu(w,1)>.2)):
            delta = values-w[a,b]
            response = tuple(map(int,np.where(np.abs(delta)>.5,np.sign(delta),0)))
            bonds.append((int(a),int(b),response))
        result.append(dict(colors=[(e,) for e in raw[side]['elements']],bonds=bonds))
    return result


class EventPatterns:
    def __init__(self,raw):
        self.raw = raw
        self.r,self.p = [np.asarray(raw[k]['wbo']) for k in ('reactant','product')]
        self.n = len(self.r)
        self.tri = np.triu(np.ones_like(self.r,dtype=bool),1)
        self.features = signed_features(raw)
        self.base = colored_graph([self.features[0]])
        self.cache = {}

    def events(self,mapping):
        vector = np.asarray(mapping,dtype=int)
        p = self.p[np.ix_(vector,vector)]
        bonded_r,bonded_p = self.r>.2,p>.2
        masks = [bonded_r & ~bonded_p,~bonded_r & bonded_p,
                 bonded_r & bonded_p & (p-self.r>.5),bonded_r & bonded_p & (self.r-p>.5)]
        return tuple(tuple((int(a),int(b)) for a,b in zip(*np.where(mask & self.tri))) for mask in masks)

    def certificate(self,events):
        if events in self.cache:
            return self.cache[events]
        base = self.base
        adjacency = {a:set(ns) for a,ns in base.adjacency_dict.items()}
        colors = [set(c) for c in base.vertex_coloring]
        count = base.number_of_vertices
        markers = list(range(count,count+4))
        for marker in markers:
            adjacency[marker] = set()
            colors.append({marker})
        count += 4
        gadgets = set()
        for marker,edges in zip(markers,events,strict=True):
            for a,b in edges:
                adjacency[count] = {a,b,marker}
                for vertex in (a,b,marker):adjacency.setdefault(vertex,set()).add(count)
                gadgets.add(count)
                count += 1
        if gadgets:colors.append(gadgets)
        graph = pynauty.Graph(count,adjacency_dict={i:sorted(adjacency.get(i,())) for i in range(count)},
                              vertex_coloring=colors)
        key = hashlib.sha256(pynauty.certificate(graph)).hexdigest()
        self.cache[events] = key
        return key

    def describe(self,mapping):
        mapping = list(map(int,mapping))
        assert sorted(mapping) == list(range(self.n))
        assert all(self.raw['reactant']['elements'][a] == self.raw['product']['elements'][b]
                   for a,b in enumerate(mapping))
        events = self.events(mapping)
        counts = [len(edges) for edges in events]
        expected = event_counts(self.r,self.p,[mapping])[0].tolist()
        assert [counts[0],counts[1],counts[2]+counts[3]] == expected
        return dict(id=self.certificate(events),mapping=mapping,events=dict(zip(KINDS,events)),
                    counts=dict(zip(KINDS,counts)),total=sum(counts))


def compare_sets(left,right):
    left,right = set(left),set(right)
    return dict(left=len(left),right=len(right),shared=sorted(left & right),
        left_only=sorted(left-right),right_only=sorted(right-left),identical=left==right,
        jaccard=len(left & right)/len(left | right) if left or right else 1.)


def prepare(args):
    args.run.mkdir(parents=True,exist_ok=False)
    for folder in ('snapshots','families','status'):(args.run/folder).mkdir()
    shutil.copy2(__file__,args.run/'driver.py')
    save(args.run/'manifest.json',dict(cases=140,aam_source=str(AAM),slap_sweep_source=str(SWEEP),
        config=read(AAM/'manifest.json')['original_config'],reference_available=False,
        score_scope='The recorded minimum full-H event score of each saved benchmark run; '
                    'not a global minimum over every unmaterialized compressed assignment.',
        equivalence='Joint bond-edit patterns on original reactant atom indices, modulo automorphisms '
                    'of the full explicit-H reactant graph with signed exact score-response edge colors. '
                    'Broken, formed, strengthened and weakened bonds remain distinct. '
                    'Numerical WBO values within an event-response class, geometry and stereochemistry are not separate mechanisms.',
        sources_sha256={str(p):sha(p) for p in [AAM/'manifest.json',AAM/'case_metrics.json',SWEEP/'manifest.json',SWEEP/'case_metrics.json']},
        driver_sha256=sha(args.run/'driver.py')))
    print(args.run,flush=True)


def snapshot(args):
    started = time.perf_counter()
    raw = read(AAM/f'inputs/{args.index}/input.json')
    canonical = EventPatterns(raw)
    baseline = read(AAM/'case_metrics.json')[args.index]
    methods = {}
    minimum = baseline['aam']['best_events']
    found = {}
    scanned,at_minimum = 0,0
    for direction in ('R_to_P','P_to_R'):
        folder = AAM/f'results/holdout/{args.index}/{direction}/original'
        graph = read_aam_checkpoint(folder/'cuts/aam.pkl.gz').graph
        selected = [t for t in graph.terminals if len(graph.states[t].mapping)==canonical.n]
        for offset in range(0,len(selected),256):
            terminals = selected[offset:offset+256]
            vectors = [[dict(graph.states[t].mapping)[a] for a in range(canonical.n)] for t in terminals]
            if direction=='P_to_R':vectors = np.argsort(vectors,axis=1).tolist()
            scores = event_counts(canonical.r,canonical.p,vectors).sum(axis=1)
            scanned += len(vectors)
            assert min(scores,default=minimum)>=minimum
            for t,vector,score in zip(terminals,vectors,scores,strict=True):
                if score != minimum:continue
                at_minimum += 1
                value = canonical.describe(vector)
                if value['id'] not in found:
                    found[value['id']] = dict(value,direction=direction,terminal=t)
    assert found
    methods['aam'] = dict(minimum=minimum,patterns=found,full_terminals_scanned=scanned,
        terminals_at_minimum=at_minimum)
    for name,root in [('native_slap',AAM),('slap_sweep',SWEEP)]:
        scored = read(root/f'slap_scoring/{args.index}.refined.json')
        minimum = scored['slap_min_representative_events']
        patterns = {}
        for candidate in scored['slap']:
            if candidate['events']['total']!=minimum:continue
            m = dict(candidate['mapping'])
            value = canonical.describe([m[a] for a in range(canonical.n)])
            assert value['total']==minimum
            if value['id'] not in patterns:patterns[value['id']] = dict(value,candidate=candidate['candidate'])
        methods[name] = dict(minimum=minimum,patterns=patterns,
            candidates_at_minimum=sum(c['events']['total']==minimum for c in scored['slap']))
    result = dict(index=args.index,name=raw['name'],methods=methods,
        comparisons={name:compare_sets(methods['aam']['patterns'],methods[name]['patterns'])
                     for name in ('native_slap','slap_sweep')},
        score_tie={name:methods['aam']['minimum']==methods[name]['minimum'] for name in ('native_slap','slap_sweep')},
        seconds=time.perf_counter()-started,scope='All saved full AAM terminal witnesses and each optimized SLAP label representative; '
            'missing representative patterns require compressed-family checks before claiming absence.')
    save(args.run/f'snapshots/{args.index}.json',result)
    print(json.dumps(dict(index=args.index,counts={m:len(v['patterns']) for m,v in methods.items()},seconds=result['seconds'])),flush=True)


def worker(args):
    output = args.run/f'status/{args.phase}_{args.index}.json'
    assert not output.exists()
    row = dict(index=args.index,phase=args.phase,started=time.time(),host=socket.gethostname())
    save(output,row)
    with output.with_suffix('.log').open('w') as log:
        command = [PYTHON,str(args.run/'driver.py'),args.phase,'--run',str(args.run),'--index',str(args.index)]
        code = subprocess.run(['timeout','--kill-after=5s',str(args.seconds),*command],stdout=log,stderr=subprocess.STDOUT).returncode
    save(output,dict(row,exit=code,finished=time.time()))


def submit(args):
    path = args.run/f'{args.phase}_submission.json'
    assert not path.exists()
    env = ['env',f'PYTHONPATH={AAM}/original/src:{AAM}/engine/bench','PYTHONHASHSEED=0',
           'PYTHONDONTWRITEBYTECODE=1','OPENBLAS_NUM_THREADS=1','OMP_NUM_THREADS=1','MKL_NUM_THREADS=1']
    command = ['sbatch','--parsable','--partition=cpunodes_nia','--exclude=bosque49,bosque56',
        '--nodes=1','--cpus-per-task=1','--mem=8G','--time=00:20:00','--no-requeue',
        '--array=0-139%32',f'--job-name=min_events_{args.phase}',f'--output={args.run}/status/{args.phase}_%A_%a.out',
        '--wrap',shlex.join([*env,PYTHON,str(args.run/'driver.py'),'worker','--run',str(args.run),
                            '--phase',args.phase,'--seconds',str(args.seconds),'--index'])+' "$SLURM_ARRAY_TASK_ID"']
    job = subprocess.check_output(command,text=True).strip()
    save(path,dict(job=job,command=command))
    print(job,flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command',choices=('prepare','snapshot','worker','submit'))
    parser.add_argument('--run',type=Path,required=True)
    parser.add_argument('--index',type=int)
    parser.add_argument('--phase',choices=('snapshot',),default='snapshot')
    parser.add_argument('--seconds',type=int,default=900)
    args = parser.parse_args()
    globals()[args.command](args)
