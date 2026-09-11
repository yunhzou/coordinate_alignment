"""Necessary-condition certificate for missing forward bond-edit patterns.

If every recorded family action preserves the product's binary heavy-atom graph,
it preserves the literal broken/formed heavy bonds of its terminal witness.
A target whose restricted edit pattern matches no such terminal is excluded.
Bond-order and hydrogen edits are deliberately ignored: matching these remaining
edits is necessary, but not sufficient, for matching the complete target pattern.
"""
import argparse
from collections import Counter
from functools import lru_cache
from pathlib import Path
import time

from holdout_missing_pattern_seeds import CASES,SEEDS
from holdout_minimum_events import EventPatterns,read,save,sha,read_aam_checkpoint,KINDS
from golden_evaluation import exact_action


def certify(run,variant,index):
    source=Path(read(run/'variants.json')[variant]['run'])
    raw=read(source/f'inputs/holdout/{index}/input.json')
    c=EventPatterns(raw);target=read(run/'targets.json')[str(index)]['pattern']
    archive=source/f'results/holdout/{index}/R_to_P/original/cuts/aam.pkl.gz'
    started,cpu=time.perf_counter(),time.process_time()
    graph=read_aam_checkpoint(archive).graph
    rheavy={i for i,e in enumerate(raw['reactant']['elements']) if e!='H'}
    pheavy={i for i,e in enumerate(raw['product']['elements']) if e!='H'}
    feature=dict(c.features[1],bonds=[(a,b,('bond',)) for a,b,_ in c.features[1]['bonds'] if a in pheavy and b in pheavy])
    actions=[]
    @lru_cache(None)
    def invariant_action(images):
        result=exact_action(images,feature)
        if result and len(actions)<64:actions.append(images)
        return result
    @lru_cache(None)
    def invariant_edge(edge):
        placement=graph.fragment_placement(edge)
        if placement is None:return True
        assert placement.target_generators is not None
        if any(not invariant_action(tuple(g.images)) for g in placement.target_generators):return False
        for domain in placement.symmetry_domains:
            if domain.source=='exact_automorph_group':continue
            for a,b in zip(domain.p_atoms,domain.p_atoms[1:]):
                images=list(range(c.n));images[a],images[b]=b,a
                if not invariant_action(tuple(images)):return False
        return True
    incoming=[[] for _ in graph.states]
    for edge in graph.transitions:incoming[edge.target].append(edge)
    @lru_cache(None)
    def can_change(state):return any(can_change(e.source) or not invariant_edge(e.id) for e in incoming[state])
    def restricted(events):
        return tuple(tuple(tuple(pair) for pair in events[k] if all(a in rheavy for a in pair))
                     if k in ('broken','formed') else () for k in KINDS)
    target_key=c.certificate(restricted(target['events']))
    signatures=Counter();noninvariant=[];full=[];samples=[];matching=[]
    for terminal in graph.terminals:
        mapping=dict(graph.states[terminal].mapping)
        if len(mapping)!=c.n:continue
        full.append(terminal)
        if can_change(terminal):noninvariant.append(terminal);continue
        vector=[mapping[a] for a in range(c.n)]
        events=dict(zip(KINDS,c.events(vector)));signature=c.certificate(restricted(events))
        signatures[signature]+=1
        if signature==target_key:matching.append(terminal)
        if len(samples)<3:samples.append(vector)
    checks=0
    for vector in samples:
        original=restricted(dict(zip(KINDS,c.events(vector))))
        for images in actions:
            moved=[images[b] for b in vector]
            assert restricted(dict(zip(KINDS,c.events(moved))))==original
            checks+=1
    excluded=bool(full) and not noninvariant and not matching
    initial=read(run/f'analysis/{variant}/{index}.json')
    if initial['target_result']['status']=='represented':assert not excluded
    if excluded:assert initial['target_result']['status']!='represented'
    # A known saved minimum witness must remain possible under this relaxation.
    minimum=initial['minimum_witness']
    assert noninvariant or c.certificate(restricted(minimum['events'])) in signatures
    result=dict(index=index,variant=variant,archive=str(archive),archive_sha256=sha(archive),target_pattern=target['id'],
        restricted_target=target_key,full_terminals=len(full),noninvariant_terminals=noninvariant,
        matching_invariant_terminals=matching,restricted_signature_counts=dict(signatures),complete_exclusion=excluded,
        action_checks=invariant_action.cache_info().misses,literal_invariance_checks=checks,
        initial_target_status=initial['target_result']['status'],cpu=time.process_time()-cpu,wall=time.perf_counter()-started,
        proof='Every path action in each invariant terminal family preserves the full product binary heavy-atom graph. '
              'Thus broken/formed heavy bonds cannot vary within that family. No target restricted pattern '
              'matches any invariant terminal modulo the same reactant symmetry; noninvariant terminals prevent exclusion.')
    save(run/f'coarse_checks/{variant}/{index}.json',result)
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run',type=Path,required=True)
    args=parser.parse_args()
    results=[certify(args.run,f'seed_{seed}',index) for seed in SEEDS for index in CASES]
    save(args.run/'coarse_validation.json',dict(cases=list(CASES),seeds=list(SEEDS),
        certificates=len(results),literal_invariance_checks=sum(r['literal_invariance_checks'] for r in results),
        known_positive_controls=sum(r['initial_target_status']=='represented' for r in results),
        false_exclusions=0,results=[dict(index=r['index'],variant=r['variant'],excluded=r['complete_exclusion'],
            full_terminals=r['full_terminals'],noninvariant=len(r['noninvariant_terminals']),
            matching=len(r['matching_invariant_terminals'])) for r in results]))
    print([(r['index'],r['variant'],r['complete_exclusion']) for r in results],flush=True)
