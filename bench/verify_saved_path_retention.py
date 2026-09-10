"""Verify an original successful compressed path inside saved candidate output.

Diagnostic only: no mapper runs and no original mapping is inserted into the
candidate. Match the complete fragment/constraint records, then independently
query the candidate path. Keep this certificate separate from timed evaluation.
"""
import argparse
from collections import defaultdict
from pathlib import Path
from types import SimpleNamespace

import pynauty

from adaptive_full_benchmark import read, save, problem_plan
from golden_evaluation import colored_graph, endpoint_generators, project
from rxn_core.alignment.sweep import cut_sweep_items
from rxn_core.artifacts import read_aam_checkpoint, read_graph_checkpoint
from rxn_core.family_query import query_path
from rxn_core.search_graph import SearchPath


def state_key(state):
    return state.mapping,state.islands,state.deferred_edges


def identical_path(path,graph):
    """Find a connected path, not independent atom/fragment membership."""
    outgoing=defaultdict(list)
    for edge in graph.transitions:
        outgoing[edge.source].append(edge)
    old=path.graph
    steps=[old.transitions[i] for i in path.transitions if old.transitions[i].match is not None]
    initial=state_key(old.states[old.transitions[path.transitions[0]].source])
    frontier={root:() for root in graph.roots if state_key(graph.states[root])==initial
              and graph.contexts[graph.states[root].context].cuts==path.context.cuts}
    for expected in steps:
        next_frontier={}
        for node,prefix in frontier.items():
            for edge in outgoing[node]:
                if (edge.match==expected.match and edge.preserved_bonds==expected.preserved_bonds
                    and state_key(graph.states[edge.target])==state_key(old.states[expected.target])):
                    next_frontier.setdefault(edge.target,prefix+(edge.id,))
        frontier=next_frontier
    terminals=set(graph.terminals)
    for node,edges in frontier.items():
        if node in terminals and state_key(graph.states[node])==state_key(old.states[path.terminal]):
            return SearchPath(graph,node,edges)
    return None


def verify(args):
    spec=dict(dataset='golden',index=args.index,direction=args.direction)
    _,plan=problem_plan(SimpleNamespace(run=args.run,method='adaptive'),spec)
    folder=args.run/f'results/golden/{args.index}/{args.direction}'
    evaluation=read(folder/'original/full_sweep_evaluation.json')
    original=read_aam_checkpoint(folder/'original/cuts/aam.pkl.gz')
    terminal=evaluation['witness_terminal']
    path=(SearchPath(original.graph,terminal,tuple(evaluation['witness_path']))
          if 'witness_path' in evaluation else next(original.graph.paths(terminal)))
    cuts=cut_sweep_items(plan.problem.reactant.wbo,plan.config.cut_floor)
    index=next(i for i,cut in enumerate(cuts) if tuple(cut)==path.context.cuts)
    archive=folder/f'adaptive/cuts/cut_{index:05d}.finalized.pkl.gz'
    graph=read_graph_checkpoint(archive)
    found=identical_path(path,graph)
    result=dict(index=args.index,direction=args.direction,original_terminal=terminal,
        candidate_archive=str(archive),identical_compressed_path=found is not None,
        scope='Post-hoc output containment certificate, separate from fixed-time reference evaluation; no search run')
    if found is not None:
        reference=read(args.run/f'inputs/golden/{args.index}/reference.json')
        features=list(reversed(reference['features'])) if plan.reversed else reference['features']
        def lifted(feature,degree):
            output=[]
            for action in endpoint_generators(feature):
                images=list(range(degree))
                for i,atom in enumerate(feature['heavy']):images[atom]=feature['heavy'][action[i]]
                output.append(tuple(images))
            return tuple(output)
        wanted=plan.to_search_mapping(reference['mapping'])
        status,witness=query_path(found,plan.problem,wanted,source_atoms=features[0]['heavy'],
            source_generators=lifted(features[0],plan.problem.source_atom_count),
            target_generators=lifted(features[1],plan.problem.target_atom_count),
            complete_reference=True,timeout_ms=5000)
        if status=='recovered':
            assert pynauty.certificate(colored_graph(features,project(witness['mapping'],features)))==\
                   pynauty.certificate(colored_graph(features,project(wanted,features)))
        result.update(reference_recovery=status,candidate_terminal=found.terminal,
                      candidate_path=list(found.transitions),witness=witness)
    save(args.run/f'diagnosis/path_retention/{args.index}_{args.direction}.json',result)
    print({k:v for k,v in result.items() if k!='witness'},flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run',type=Path,required=True)
    parser.add_argument('--index',type=int,required=True)
    parser.add_argument('--direction',choices=['R_to_P','P_to_R'],required=True)
    verify(parser.parse_args())
