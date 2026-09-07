"""Exact necessary source-orbit check before saved-family verification.

This is a separate diagnostic, not a replacement top-1 benchmark ranking.
Target actions cannot change the set of mapped source atoms. Endpoint source
automorphisms cannot move atoms between their exact chemical orbits.
"""
import argparse
from collections import Counter
from dataclasses import replace
import gc
import gzip
import json
from pathlib import Path
import time
from functools import partial

import pynauty

from golden_evaluation import colored_graph,project,evaluate,symbolic_path_query,endpoint_generators
from investigate_golden_mapping import save
from rxn_core.artifacts import read_aam,read_aam_checkpoint


def source_orbit_key(mapping, orbits):
    return tuple(sorted(Counter(orbits[r] for r in mapping).items()))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--archive',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--finite',action='store_true')
    p.add_argument('--query-seconds',type=float,default=10)
    p.add_argument('--target-orbits',action='store_true')
    args=p.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    started=time.perf_counter()
    if args.archive.name.endswith('.pkl.gz'):aam=read_aam_checkpoint(args.archive)
    else:
        with gzip.open(args.archive,'rt') as f:aam=read_aam(f)
    gc.disable()
    reference=json.loads((args.archive.parent/'reference.json').read_text())
    features=reference['features'];expected=project(reference['mapping'],features)
    assert len(set(expected.values()))==len(features[1]['heavy'])
    orbits=pynauty.autgrp(colored_graph([features[0]]))[3][:len(features[0]['heavy'])]
    expected_key=source_orbit_key(expected,orbits)
    terminals=aam.graph.terminals
    eligible={t for t in terminals if source_orbit_key(project(aam.graph.states[t].mapping,features),orbits)==expected_key}
    source_eligible=len(eligible)
    if args.target_orbits and eligible:
        # Overapproximate ALL returned target actions by their union orbits.
        # Failure is conclusive; passing is not an acceptance criterion.
        degree=len(features[1]['heavy']);parents=list(range(degree))
        def root(a):
            while parents[a]!=a:
                parents[a]=parents[parents[a]];a=parents[a]
            return a
        def action(g):
            for a,b in enumerate(g):parents[root(a)]=root(b)
        for g in endpoint_generators(features[1]):action(g)
        index={a:i for i,a in enumerate(features[1]['heavy'])}
        ancestors=aam.graph.ancestor_transitions(eligible)
        for edge in aam.graph.transitions:
            if edge.id not in ancestors or edge.match is None:continue
            for g in edge.match['symmetry']['automorph_generators']:
                action(tuple(index[g[a]] for a in features[1]['heavy']))
        target_orbits=[root(i) for i in range(degree)]
        def pair_key(mapping):
            return tuple(sorted(Counter((orbits[r],target_orbits[p]) for r,p in mapping.items()).items()))
        required=pair_key(expected)
        eligible={t for t in eligible if pair_key(project(aam.graph.states[t].mapping,features))==required}
    save(args.output/'filter.json',dict(original_terminals=len(terminals),eligible_terminals=len(eligible),
        source_eligible_terminals=source_eligible,target_orbits_checked=args.target_orbits,
        rejected_by_necessary_condition=len(terminals)-len(eligible),seconds=time.perf_counter()-started))
    view=replace(aam,graph=replace(aam.graph,stops=tuple(s for s in aam.graph.stops
        if s.reason not in {'objective_met','stalled'} or s.state in eligible)))
    checked=evaluate(view,features,reference['mapping'],seconds=120,
        query=partial(symbolic_path_query,finite_domain=True) if args.finite else None,
        query_timeout_ms=1000*args.query_seconds)
    save(args.output/'verification.json',dict(reference_recovery=checked['reference_recovery'],
        witness_terminal=checked.get('witness_terminal'),witness_path=checked.get('witness_path'),
        witness_group_edges=checked.get('witness_group_edges'),witness_actions=checked.get('witness_actions'),
        symbolic_queries=checked['symbolic_queries'],unknown_queries=checked['unknown_queries'],
        evaluation_seconds=checked['evaluation_seconds'],total_seconds=time.perf_counter()-started,
        original_archive=str(args.archive),top1_not_recomputed=True,finite_domain=args.finite,
        query_seconds=args.query_seconds))


if __name__=='__main__':main()
