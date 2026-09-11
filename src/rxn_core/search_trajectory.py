"""Reusable, verified replay of fragment calls along saved AAM search paths.

Run capture in an isolated process: diagnostic observers temporarily wrap the
matcher while preserving its decisions. No sweep or seed-order search is run.
"""
import hashlib
import importlib
import gzip
import json
from pathlib import Path

import numpy as np
from .artifacts import read_aam_checkpoint, aam_from_record
from .frag import build_graph
from .matcher import _nauty_orbits, _cand_possible_p_atoms, _support_witness_for_value, _sym_block_indexes
from .matcher.support import _refine_sym_assignments
from .growth.trace import cands_pattern_sample

growth = importlib.import_module(__package__ + '.growth.island')
extension = importlib.import_module(__package__ + '.matcher.extend')


def snapshot(cands):
    result = cands_pattern_sample(cands, len(cands))
    for row,cand in zip(result,cands):
        row['automorph_blocks'] = [dict(r_atoms=list(b.r_atoms),p_atoms=list(b.p_atoms))
                                  for b in cand.automorph_blocks]
    return result


def replay(source, target, state, edge, po, ro, context):
    calls, events, quotient = [], [], []
    original_extend, original_dedupe = growth._extend_sym_cands, extension._dedupe_children

    def dedupe(children, ctx):
        result = original_dedupe(children, ctx)
        quotient.append(dict(before=len(children), after=len(result),
                             children=snapshot(children)))
        return result

    def extend(cands, fragment, n, r, p, mapping, tol, islands, **kwargs):
        before = snapshot(cands)
        trial_groups = []
        if n not in mapping:
            for ci, cand in enumerate(cands):
                bonded=sorted(a for a in fragment if r.has_edge(a,n))
                weights=[(a,float(r.graph['wbo_matrix'][a,n])) for a in bonded]
                anchor=kwargs.get('anchor_u')
                strict={anchor:kwargs['anchor_wbo']} if anchor in fragment else {}
                _,p_to_block=_sym_block_indexes(cand)
                occupied = {v:k for k,v in mapping.items()}
                local = {v:k for k,v in cand.items()}
                possible = _cand_possible_p_atoms(cand)
                rows = []
                for v in p:
                    if p.nodes[v]['element'] != r.nodes[n]['element']:
                        continue
                    row = dict(target=v)
                    if v in occupied:
                        row.update(status='occupied', reason=f'Already committed to R{occupied[v]}')
                    elif v in local:
                        row.update(status='occupied', reason=f'Already used by R{local[v]} in this candidate')
                    elif v in possible:
                        row.update(status='symmetry', reason='Reserved within this candidate symmetry block')
                    else:
                        support=_support_witness_for_value(cand,n,v,bonded,weights,p,tol,
                            join_block_idx=p_to_block.get(v),strict_r_wbos=strict)
                        refined=_refine_sym_assignments(cand,support) if support is not None else None
                        if refined is None:support=None
                        row.update(status='rejected' if support is None else 'compatible',
                            reason='No supported assignment satisfies the active bond checks.' if support is None else 'Active bond checks pass.')
                        if support is not None:row['witness']={**dict(refined.items()),n:v}
                    witness=row.get('witness',cand)
                    row['bonds'] = [dict(r=[a,n], p=[witness[a],v],
                        wbo=[float(r.graph['wbo_matrix'][a,n]),float(p.graph['wbo_matrix'][witness[a],v])],
                        delta=abs(float(r.graph['wbo_matrix'][a,n])-float(p.graph['wbo_matrix'][witness[a],v])))
                        for a in sorted(fragment) if r.has_edge(a,n)]
                    if row['status']=='compatible':
                        merged={**mapping,**row['witness']}
                        if len(set(merged.values()))!=len(merged) or any(
                            bond['wbo'][1]<context.graph_floor or bond['delta']>tol for bond in row['bonds']):
                            raise ValueError('Invalid compatible-trial witness')
                    rows.append(row)
                trial_groups.append(rows)
        quotient.clear()
        result = original_extend(cands, fragment, n, r, p, mapping, tol, islands, **kwargs)
        calls.append(dict(fragment=sorted(fragment), added=n,
            edge=[kwargs['anchor_u'],n], before=before, after=snapshot(result),
            trials=trial_groups, dedupe=list(quotient)))
        return result

    growth._extend_sym_cands, extension._dedupe_children = extend, dedupe
    try:
        isos = growth.grow_island(source,target,edge.seed,dict(state.mapping),
            graph_floor=context.graph_floor,iso_tol=context.iso_tolerance,max_branches=context.branch_limit,events=events,
            islands_R=dict(state.islands),p_orbits=po,r_orbits=ro,
            prior_deferred_edges=state.deferred_edges)
    finally:
        growth._extend_sym_cands, extension._dedupe_children = original_extend, original_dedupe
    expected = dict(edge.match['symmetry']['witness'])
    matching = [i for i in isos if dict(i)==expected and sorted(i.fragment)==sorted(edge.match['fragment'])]
    if not matching or sorted(matching[0].deferred_edges)!=sorted(map(tuple,edge.match['deferred_edges'])):
        raise ValueError(f'Replay differs from archived transition {edge.id}; use its producing engine')
    return dict(events=events,calls=calls,isos=[dict(i) for i in isos],verified=True)


def bond_events(raw, mapping, floor=.2, tolerance=.5):
    r=np.array(raw['reactant']['wbo']);p=np.array(raw['product']['wbo'])
    out=[]
    for a in range(len(r)):
        for b in range(a+1,len(r)):
            x,y=float(r[a,b]),float(p[mapping[a],mapping[b]])
            kind=('broken' if x>floor and y<=floor else 'formed' if x<=floor and y>floor else
                  'weakened' if x>floor and y>floor and x-y>tolerance else
                  'strengthened' if x>floor and y>floor and y-x>tolerance else None)
            if kind:out.append(dict(kind=kind,r=[a,b],p=[mapping[a],mapping[b]],wbo=[x,y]))
    return out


def frames_for(graph, terminal, traces):
    frames=[];stage=0
    path=next(graph.paths(terminal))
    atom_count=len(path.context.source_atoms)
    for stage,eid in enumerate(path.transitions,1):
        edge=graph.transitions[eid];state=graph.states[edge.source]
        if edge.match is None:
            after=graph.states[edge.target]
            frames.append(dict(kind='locked',stage=stage,title='Reuse recorded state',locked=dict(after.mapping),
                active=[],deferred=list(map(list,after.deferred_edges)),candidates=[],preferred=0))
            continue
        trace=traces[eid];calls=iter(trace['calls']);call=None
        locked=dict(state.mapping);deferred=list(map(list,state.deferred_edges))
        expected=dict(edge.match['symmetry']['witness'])
        patterns=[]
        for event in trace['events']:
            kind=event['type']
            if kind=='pop_skip':continue
            data=dict(stage=stage,seed=edge.seed,transition=eid,locked=locked,
                      deferred=[*deferred],kind=kind,expected=expected)
            if kind=='seed_start':
                patterns=trace['calls'][0]['before'] if trace['calls'] else event.get('cand_patterns',[])
                data.update(title=f'Start fragment {stage} at R{edge.seed}',active=event['fragment'],candidates=patterns)
            elif kind=='pop':
                call=next(calls);assert call['added']==event['edge']['ext_atom']
                patterns=call['before']
                data.update(title=f'Test R{call["edge"][0]}–R{call["edge"][1]}',
                    active=call['fragment'],candidates=patterns,trial=call,edge=call['edge'])
            elif kind in ('commit','consumed'):
                assert call is not None
                patterns=call['after'] if kind=='commit' else call['before']
                if kind=='consumed':deferred.append(call['edge'])
                data.update(title=(f'Grow to R{call["added"]}' if kind=='commit' else
                    f'Defer R{call["edge"][0]}–R{call["edge"][1]}'),
                    active=event['fragment'],candidates=patterns,edge=call['edge'],
                    refinement=dict(before=len(call['before']),after=len(call['after']),dedupe=call['dedupe']),
                    deferred=[*deferred])
            elif kind=='seed_end':
                patterns=[dict(witness=i,blocks=[]) for i in trace['isos']]
                data.update(title=f'Fragment {stage} saturated · {len(patterns)} retained placement(s)',
                    active=event['fragment'],candidates=patterns)
            else:continue
            # All previews are actual live representatives, never a made-up completion.
            data['preferred']=max(range(len(patterns)),key=lambda i:sum(expected.get(a)==b for a,b in patterns[i]['witness'].items())) if patterns else 0
            frames.append(data)
        after=graph.states[edge.target]
        frames.append(dict(kind='locked',stage=stage,seed=edge.seed,transition=eid,
            title=f'Commit fragment {stage} · {len(after.mapping)} / {atom_count} atoms assigned',
            locked=dict(after.mapping),active=[],deferred=list(map(list,after.deferred_edges)),
            candidates=[],preferred=0,expected=expected))
    final=dict(graph.states[terminal].mapping)
    frames.append(dict(kind='terminal',stage=stage,title='Saved terminal mapping · final bond events',
        locked=final,active=[],deferred=list(map(list,graph.states[terminal].deferred_edges)),
        candidates=[],preferred=0))
    return frames


def endpoint_record(endpoint):
    if endpoint.coordinates is None:
        raise ValueError('A 3D search trajectory requires endpoint coordinates in the archive')
    return dict(elements=list(endpoint.elements),coordinates=endpoint.coordinates.tolist(),wbo=endpoint.wbo.tolist())


def read_archive(path):
    if path.name.endswith(('.json','.json.gz')):
        opener=gzip.open if path.suffix=='.gz' else open
        with opener(path,'rt') as stream:return aam_from_record(json.load(stream))
    return read_aam_checkpoint(path)


def build_trajectory(selections, *, title=None, key_atoms=None, watch_targets=None, event_tolerance=.5):
    """Capture selected contexts from compatible archives of the same endpoints.

    Each selection specifies archive + context, with optional terminals, label,
    path_labels, default_terminal and focus={source_edge:[a,b], target:p}.
    One recorded history is followed per terminal; compressed permutations are
    represented, not expanded. A replay mismatch raises instead of publishing.
    """
    runs,checks,raw=[],[],None
    for selection in selections:
        archive=Path(selection['archive'])
        result=read_archive(archive);g=result.graph
        context_id=int(selection.get('context',0));c=g.contexts[context_id]
        current=dict(name=result.problem.name,
                     reactant=endpoint_record(result.problem.reactant),product=endpoint_record(result.problem.product))
        if raw is None:raw=current
        elif any(raw[k]!=current[k] for k in ('reactant','product')):
            raise ValueError('Comparison archives must have identical endpoint atom order, coordinates and WBOs')
        source,target=[build_graph(raw[k]['elements'],np.array(raw[k]['wbo']),bond_cut=c.graph_floor)
                       for k in ('reactant','product')]
        source.remove_edges_from(c.cuts)
        po=_nauty_orbits(target,wbo_tol=c.iso_tolerance);ro=_nauty_orbits(source,wbo_tol=c.iso_tolerance)
        terminals=selection.get('terminals',[t for t in g.terminals if g.states[t].context==context_id])
        if not terminals:raise ValueError(f'Context {context_id} has no selected terminal paths')
        if any(g.states[t].context!=context_id or t not in g.terminals for t in terminals):
            raise ValueError('Selected terminal does not belong to the selected context')
        ids={eid for t in terminals for eid in next(g.paths(t)).transitions}
        traces={eid:replay(source,target,g.states[g.transitions[eid].source],g.transitions[eid],po,ro,c)
                for eid in sorted(ids) if g.transitions[eid].match is not None}
        paths=[]
        for t in terminals:
            frames=frames_for(g,t,traces);mapping=dict(g.states[t].mapping)
            complete=len(mapping)==len(raw['reactant']['elements'])
            events=bond_events(raw,mapping,c.graph_floor,event_tolerance) if complete else []
            for frame in frames:
                for candidate in frame['candidates']:
                    preview={**frame['locked'],**candidate['witness']}
                    if len(set(preview.values()))!=len(preview) or any(
                        raw['reactant']['elements'][a]!=raw['product']['elements'][b] for a,b in preview.items()):
                        raise ValueError('Invalid replayed candidate mapping')
            paths.append(dict(terminal=t,label=selection.get('path_labels',{}).get(str(t),f'Terminal {t}'),
                              frames=frames,mapping=mapping,events=events,complete=complete))
            checks.append(dict(archive=str(archive),context=context_id,terminal=t,frames=len(frames),
                               archived_fragment_matches='exact',complete=complete,events=len(events)))
        nodes=[s for s in g.states if s.context==context_id]
        runs.append(dict(label=selection.get('label',f'Context {context_id} · R{c.seed_order[0]} first'),
            context=context_id,seed_order=list(c.seed_order),cuts=list(c.cuts),paths=paths,
            default_terminal=selection.get('default_terminal',terminals[0]),focus=selection.get('focus'),
            config=dict(graph_floor=c.graph_floor,iso_tolerance=c.iso_tolerance,branch_limit=c.branch_limit,
                        event_tolerance=event_tolerance),
            graph=dict(states=[dict(id=s.id,assigned=len(s.mapping)) for s in nodes],
                transitions=[dict(id=eid,source=g.transitions[eid].source,target=g.transitions[eid].target,
                    seed=g.transitions[eid].seed,size=len(g.transitions[eid].match['fragment']) if g.transitions[eid].match else 0)
                    for eid in sorted(ids)]),
            provenance=dict(archive=str(archive),sha256=hashlib.sha256(archive.read_bytes()).hexdigest()),traces=traces))
    if raw is None:raise ValueError('At least one archive selection is required')
    default_keys=sorted({a for run in runs for path in run['paths'] for event in path['events'] for a in event['r']})
    return dict(schema='rxn_core.search_trajectory/v1',name=title or raw['name'] or 'AAM search trajectory',
        input=raw,runs=runs,key_atoms=default_keys if key_atoms is None else key_atoms,
        watch_targets=[] if watch_targets is None else watch_targets,
        capture_source={str(p):hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(set(Path(growth.__file__).parents[1].rglob('*.py'))|{Path(__file__)})},
        scope='Diagnostic replay of selected saved fragment calls; one recorded history per terminal. '
              'All live compressed candidates are included, each with an actual representative and symmetry blocks. '
              'Trial rows inspect support before deduplication. Every fragment placement is checked against its archive. '
              'The endpoint coordinates stay fixed throughout playback.',checks=checks)
