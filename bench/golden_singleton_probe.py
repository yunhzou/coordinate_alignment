"""Positive-only check of saved trailing singleton domains, not a complete family verifier."""
import argparse
import json
from pathlib import Path
import time
from itertools import chain

import pynauty
from golden_policy_campaign import load_case,save
from golden_evaluation import project,colored_graph,symbolic_path_query
from rxn_core.artifacts import read_aam_checkpoint


def main(args):
    source,plan=load_case(args.source,args.index)
    out=args.output/str(args.index);out.mkdir(parents=True,exist_ok=True)
    archive=source/'cuts/aam.pkl.gz'
    if not archive.exists():
        save(out/'result.json',dict(index=args.index,reference_recovery='unknown',reason='no full checkpoint'))
        return
    start=time.perf_counter();aam=read_aam_checkpoint(archive)
    ref=json.loads((source/'reference.json').read_text())
    features=list(reversed(ref['features'])) if plan.reversed else ref['features']
    expected=project(plan.to_search_mapping(ref['mapping']),features)
    certificate=pynauty.certificate(colored_graph(features,expected))
    ri={a:i for i,a in enumerate(features[0]['heavy'])};pi={a:i for i,a in enumerate(features[1]['heavy'])}
    seen=set();queries=0;unknown=0
    top=json.loads((source/'evaluation.json').read_text()).get('top_terminal')
    paths=chain(aam.graph.paths(top),aam.graph.paths()) if top is not None else aam.graph.paths()
    for path in paths:
        if time.perf_counter()-start>150:break
        domains={}
        for edge_id in reversed(path.transitions):
            edge=aam.graph.transitions[edge_id]
            if edge.match is None:continue
            fragment=edge.match['fragment']
            if len(fragment)!=1:break
            atom=fragment[0]
            if atom not in ri:continue
            for block in edge.match['symmetry']['blocks']:
                if block['r_atoms']==[atom] and atom not in edge.match['symmetry']['exact_fixed']:
                    domains[ri[atom]]=tuple(pi[p] for p in block['p_atoms'])
        if not domains:continue
        mapping=project(path.mapping,features)
        key=(tuple(sorted(mapping.items())),tuple(sorted(domains.items())))
        if key in seen:continue
        seen.add(key);queries+=1
        # Linear, one-domain counterexample probe, not a Cartesian expansion
        # of assignments. Each positive is independently certified below.
        witness=None
        for atom,domain in domains.items():
            for target in domain:
                trial={**mapping,atom:target}
                if len(set(trial.values()))!=len(trial):continue
                if pynauty.certificate(colored_graph(features,trial))==certificate:
                    witness=dict(heavy_mapping=sorted(trial.items()),singleton_assignments={atom:target},
                                 method='one saved-domain substitution, exact certificate')
                    break
            if witness is not None:break
        if witness is not None:status='recovered'
        else:
            status,witness=symbolic_path_query(mapping,[],features,expected,2000,
                finite_domain=True,independent_singleton_domains=domains)
        if status=='unknown':unknown+=1
        if status!='recovered':continue
        heavy=dict(witness['heavy_mapping'])
        # Certificate is checked independently of SAT. Keep every non-singleton
        # assignment and explicit-H assignment exactly as archived.
        assert pynauty.certificate(colored_graph(features,heavy))==certificate
        full=dict(path.mapping)
        for r,p in heavy.items():full[features[0]['heavy'][r]]=features[1]['heavy'][p]
        assert len(set(full.values()))==len(full)
        for atom in full:
            assert plan.problem.reactant.elements[atom]==plan.problem.product.elements[full[atom]]
        for edge_id in path.transitions:
            edge=aam.graph.transitions[edge_id]
            if edge.match is None:continue
            atoms=edge.match['fragment']
            if len(atoms)>1:assert all(full[r]==path.mapping[r] for r in atoms)
        save(out/'result.json',dict(index=args.index,reference_recovery='recovered',
            archive=str(archive),terminal=path.terminal,transitions=path.transitions,
            domains=domains,witness=witness,input_orientation_witness=sorted(plan.to_input_mapping(full).items()),
            seconds=time.perf_counter()-start,queries=queries,
            scope='Verified saved trailing singleton-domain witness; no search rerun; no top-1 claim'))
        return
    save(out/'result.json',dict(index=args.index,reference_recovery='unknown',
        queries=queries,unknown_queries=unknown,seconds=time.perf_counter()-start,
        scope='Positive-only trailing singleton probe; absence is not a family miss'))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',type=Path,required=True);p.add_argument('--index',type=int,required=True)
    p.add_argument('--output',type=Path,required=True);main(p.parse_args())
