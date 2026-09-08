"""Reference-blind class ranking and separate publication accuracy queries."""
from dataclasses import replace
import hashlib
import json
from pathlib import Path
import time

import pynauty
from golden_evaluation import colored_graph, project, evaluate_planned
from golden_topk import summarize_classes
from golden_policy_campaign import save
from rxn_core.artifacts import read_aam_checkpoint
from view_golden_remaining import ranker


def certificate_id(features, mapping):
    return hashlib.sha256(pynauty.certificate(colored_graph(features,project(mapping,features)))).hexdigest()


def rank_archive(aam, plan):
    """No reference input. Canonical classes in original R/P space."""
    features = [dict(e.metadata)['scoring_features'] for e in
                (plan.input_problem.reactant,plan.input_problem.product)]
    score = ranker(plan.input_problem)
    cache, groups = {}, {}
    for terminal in aam.graph.terminals:
        mapping = plan.to_input_mapping(aam.graph.states[terminal].mapping)
        heavy = tuple(sorted(project(mapping,features).items()))
        if heavy not in cache:
            cache[heavy] = certificate_id(features,mapping)
        cid = cache[heavy]
        key, events = score(mapping)
        if cid not in groups:
            groups[cid] = dict(id=cid,key=key,events=events,terminal=terminal,
                              direction=plan.direction,members=[])
        row = groups[cid]
        row['members'].append(terminal)
        if key < row['key']:
            row.update(key=key,events=events,terminal=terminal)
    return sorted(groups.values(),key=lambda r:(r['key'],r['id']))


def merge_classes(ranked):
    """Union indexes only; never invert a compressed graph's generators."""
    groups = {}
    for classes in ranked:
        for row in classes:
            cid = row['id']
            entry = dict(direction=row['direction'],terminals=row['members'])
            if cid not in groups:
                groups[cid] = {**row,'origins':[entry]}
            else:
                old = groups[cid]
                origins = old['origins']+[entry]
                best = row if (row['key'],row['direction']) < (old['key'],old['direction']) else old
                groups[cid] = {**best,'origins':origins}
    return sorted(groups.values(),key=lambda r:(r['key'],r['id']))


def representative_metrics(classes, expected):
    rows = [dict(rank=i+1,score=row['key'][:3],reference_equivalent=row['id']==expected)
            for i,row in enumerate(classes)]
    return dict(class_count=len(rows),
        topk={str(k):any(r['reference_equivalent'] for r in rows[:k]) for k in (1,3,5,10)},
        **summarize_classes(rows))


def union_outcome(outcomes):
    if 'recovered' in outcomes:
        return 'recovered'
    if outcomes and all(o == 'not_recovered' for o in outcomes):
        return 'not_recovered'
    return 'unknown'


def top_five_family(run, index, classes, plans, reference, *, seconds=220):
    """Bounded verification of selected classes, not reference-guided ranking."""
    start = time.monotonic()
    expected = certificate_id(reference['features'], reference['mapping'])
    if any(c['id'] == expected for c in classes[:5]):
        return dict(outcome='recovered',reason='representative',seconds=time.monotonic()-start)
    eligible = {}
    for c in classes[:5]:
        for origin in c['origins']:
            eligible.setdefault(origin['direction'],set()).update(origin['terminals'])
    outcomes = []
    checks = []
    for direction, terminals in sorted(eligible.items()):
        if time.monotonic()-start >= seconds:
            outcomes.append('unknown')
            break
        archive = run/'directions'/str(index)/direction/'cuts/aam.pkl.gz'
        aam = read_aam_checkpoint(archive)
        view = replace(aam,graph=replace(aam.graph,stops=tuple(s for s in aam.graph.stops
            if s.reason not in {'objective_met','stalled'} or s.state in terminals)))
        left = max(.01,seconds-(time.monotonic()-start))
        check = evaluate_planned(view,plans[direction],reference['features'],reference['mapping'],
                                  seconds=left,query_timeout_ms=5000)
        outcomes.append(check['reference_recovery'])
        checks.append(dict(direction=direction,evaluation=check))
        if outcomes[-1] == 'recovered':
            break
        del view,aam
    return dict(outcome=union_outcome(outcomes),checks=checks,seconds=time.monotonic()-start)
