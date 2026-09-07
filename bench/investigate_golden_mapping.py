"""Read-only diagnosis of saved Golden misses; never rerun or alter AAM."""
import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
import csv
import gc
import gzip
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import time

import pynauty
from rdkit import Chem

from golden_evaluation import colored_graph, exact_action, project
from rxn_core.artifacts import read_aam, read_aam_checkpoint


def save(path, value):
    temporary=path.with_suffix('.tmp')
    temporary.write_text(json.dumps(value,indent=2)+'\n')
    temporary.replace(path)


def original_reference_certificate(smiles):
    """Build directly in RDF atom order, without the input canonical reorder."""
    features=[];labels=[]
    for side in smiles.split('>>'):
        mol=Chem.MolFromSmiles(side)
        map_labels={a.GetIdx():a.GetAtomMapNum() for a in mol.GetAtoms()}
        for atom in mol.GetAtoms():atom.SetAtomMapNum(0)
        mol=Chem.AddHs(mol)
        Chem.AssignStereochemistry(mol,cleanIt=True,force=True)
        heavy=[a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum()!=1]
        index={a:i for i,a in enumerate(heavy)}
        labels.append({label:index[i] for i,label in map_labels.items() if label and i in index})
        colors=[]
        for i in heavy:
            a=mol.GetAtomWithIdx(i)
            colors.append((a.GetSymbol(),a.GetFormalCharge(),a.GetIsotope(),
                a.GetTotalNumHs(includeNeighbors=True),a.GetProp('_CIPCode') if a.HasProp('_CIPCode') else ''))
        bonds=[(index[b.GetBeginAtomIdx()],index[b.GetEndAtomIdx()],
                (b.GetBondTypeAsDouble(),str(b.GetStereo()))) for b in mol.GetBonds()
               if b.GetBeginAtomIdx() in index and b.GetEndAtomIdx() in index]
        features.append(dict(heavy=heavy,colors=colors,bonds=bonds))
    reference={labels[0][k]:labels[1][k] for k in labels[0].keys()&labels[1].keys()}
    return pynauty.certificate(colored_graph(features,reference))


def initialize(args):
    args.run.mkdir(parents=True,exist_ok=False)
    rows=list(csv.DictReader(args.report.open()))
    selected=[r for r in rows if r['reference_complete']=='True' and
              r['reference_recovery']=='not_recovered' and int(r['candidate_terminals'])>0]
    assert len(selected)==137
    raw={r['index']:r for r in map(json.loads,args.audit.read_text().splitlines())}
    save(args.run/'manifest.json',dict(records=selected,audit=str(args.audit),
         source_report=str(args.report),purpose='137 nonempty ground-truth misses'))
    for row in selected:
        directory=args.run/row['index'];directory.mkdir()
        save(directory/'case.json',dict(row=row,mapped_reaction=raw[int(row['index'])]['mapped_reaction']))


def investigate(args):
    started=time.perf_counter();directory=args.run/str(args.index)
    item=json.loads((directory/'case.json').read_text());row=item['row']
    reference=json.loads((Path(row['evaluation']).parent/'reference.json').read_text())
    features=reference['features'];expected=project(reference['mapping'],features)
    original_ok=(original_reference_certificate(item['mapped_reaction'])==
                 pynauty.certificate(colored_graph(features,expected)))
    save(directory/'progress.json',dict(stage='load',original_reference_verified=original_ok))
    path=Path(row['archive'])
    if path.name.endswith('.pkl.gz'):result=read_aam_checkpoint(path)
    else:
        with gzip.open(path,'rt') as stream:result=read_aam(stream)
    gc.disable()
    graph=result.graph;terminals=graph.terminals
    loaded=time.perf_counter()
    save(directory/'progress.json',dict(stage='inspect',load_seconds=loaded-started))
    selected=graph.ancestor_transitions(terminals)
    target_index={a:i for i,a in enumerate(features[1]['heavy'])}
    actions=set();unfinalized=0
    for edge in graph.transitions:
        if edge.id not in selected or edge.match is None:continue
        symmetry=edge.match['symmetry']
        if symmetry.get('automorph_group_source')!='conditioned_search_transition':
            unfinalized+=1;continue
        for g in symmetry['automorph_generators']:
            actions.add(tuple(target_index[g[a]] for a in features[1]['heavy']))
    nonexact=[g for g in actions if not exact_action(g,features[1])]
    orbits=[pynauty.autgrp(colored_graph([f]))[3][:len(f['heavy'])] for f in features]
    def source_key(mapping):return tuple(sorted(Counter(orbits[0][r] for r in mapping).items()))
    def pair_key(mapping):return tuple(sorted(Counter((orbits[0][r],orbits[1][p]) for r,p in mapping.items()).items()))
    expected_source=source_key(expected);expected_pairs=pair_key(expected)
    exact_source=exact_pairs=0;best=None;best_distance=None;keys=set()
    expected_counts=Counter(dict(expected_pairs))
    for terminal in terminals:
        mapping=project(graph.states[terminal].mapping,features)
        key=pair_key(mapping);keys.add(key)
        exact_source+=source_key(mapping)==expected_source
        exact_pairs+=key==expected_pairs
        actual=Counter(dict(key))
        distance=sum((actual-expected_counts).values())+sum((expected_counts-actual).values())
        if best_distance is None or distance<best_distance:best_distance=distance;best=terminal
    save(directory/'diagnosis.json',dict(index=args.index,
        original_reference_verified=original_ok,terminals=len(terminals),
        contexts=len(graph.contexts),states=len(graph.states),transitions=len(graph.transitions),
        cap_stages=dict(Counter(s.stage for s in graph.stops if s.reason=='capped')),
        stop_reasons=dict(Counter(s.reason for s in graph.stops)),
        unique_projected_actions=len(actions),non_endpoint_automorphisms=len(nonexact),
        unfinalized_terminal_edges=unfinalized,
        reference_source_orbit_selection_matches=exact_source,
        reference_atom_orbit_pairing_matches=exact_pairs,
        unique_orbit_pairings=len(keys),closest_orbit_distance=best_distance,
        closest_terminal=best,load_seconds=loaded-started,total_seconds=time.perf_counter()-started,
        slurm_job=os.environ.get('SLURM_JOB_ID'),hostname=os.uname().nodename))


def pool(args):
    rows=json.loads((args.run/'manifest.json').read_text())['records'][args.shard::args.shards]
    def one(row):
        directory=args.run/row['index']
        if (directory/'diagnosis.json').exists():return
        with (directory/'diagnose.log').open('a') as log:
            child=subprocess.Popen([sys.executable,__file__,'case','--run',str(args.run),
                '--index',row['index']],stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
            try:code=child.wait(timeout=300)
            except subprocess.TimeoutExpired:
                os.killpg(child.pid,signal.SIGKILL);child.wait();code='timeout'
        save(directory/'status.json',dict(exit_code=code))
    with ThreadPoolExecutor(max_workers=args.workers) as executor:list(executor.map(one,rows))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode',choices=['init','case','pool'])
    parser.add_argument('--run',type=Path,required=True)
    parser.add_argument('--report',type=Path)
    parser.add_argument('--audit',type=Path)
    parser.add_argument('--index',type=int)
    parser.add_argument('--shard',type=int,default=0)
    parser.add_argument('--shards',type=int,default=1)
    parser.add_argument('--workers',type=int,default=2)
    args=parser.parse_args();dict(init=initialize,case=investigate,pool=pool)[args.mode](args)
