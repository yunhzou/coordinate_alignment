"""Reactant-to-product scores, timings and alternatives from saved holdout runs.

No mapper is rerun. Candidate provenance filters SLAP before scoring; only the
forward AAM archive supplies witnesses. Reverse outputs can provide an exclusion
bound, but never a forward positive witness.
"""
import argparse
from collections import Counter
from functools import lru_cache
from pathlib import Path
import shlex
import shutil
import socket
import subprocess
import time
from types import SimpleNamespace

from holdout_minimum_events import *
from holdout_minimum_event_families import aam_membership, exact_action
from enumerate_slap_minimum_events import enumerate_family
from enumerate_aam_minimum_events import case as enumerate_aam

PRIOR = DATA/'holdout_minimum_event_patterns_20260910'


def forward_candidates(index):
    origins = read(SWEEP/f'evaluations/{index}_sources.json')['sources']
    native = read(SWEEP/f'slap_xyz/{index}.json')['candidates']
    scored = read(SWEEP/f'slap_scoring/{index}.refined.json')['slap']
    assert len(origins)==len(native)==len(scored)
    return [(j,native[j],scored[j],sources) for j,sources in enumerate(origins)
            if any(source['direction']=='R_to_P' for source in sources)]


def directional_timing(index,baseline):
    folder = SWEEP/f'outputs/{index}/R_to_P'
    header = read(folder/'input.json')
    journal = [json.loads(line) for line in (folder/'records.jsonl').read_text().splitlines()]
    assert read(folder/'complete.json')['complete']
    assert len(journal)==header['expected_calls']==len(header['edges'])+1
    assert [row['ordinal'] for row in journal]==list(range(len(journal)))
    assert [row['cut'] for row in journal]==[None,*header['edges']]
    assert all(row['status']=='mapped' and row['candidate_count']>0 for row in journal)
    assert read(SWEEP/f'status/{index}.json')['host']==baseline['environment']['host']
    d = baseline['aam']['directions']['R_to_P']
    assert d['search_complete'] and d['evaluation_complete']
    search = read(AAM/f'results/holdout/{index}/R_to_P/original/search.json')['rows'][-1]
    assert d['cpu_seconds']==search['compute_cpu_excluding_persistence_and_loading_seconds']
    return dict(aam_cpu=d['cpu_seconds'],native_slap_cpu=baseline['slap']['cpu_seconds'],
        slap_sweep_cpu=header['preparation_cpu']+sum(row['graph_cpu']+row['mapping_cpu']+row['export_cpu'] for row in journal),
        aam_wall=d['search_wall_including_io'],native_slap_wall=baseline['slap']['search_wall_including_io'],
        slap_sweep_wall=header['preparation_wall']+sum(row['graph_wall']+row['mapping_wall']+row['export_wall'] for row in journal),
        aam_scoring_cpu=d['scoring_cpu_seconds'],native_slap_scoring_cpu=baseline['slap']['scoring_cpu_seconds'],
        slap_sweep_scoring_cpu=None,calls=len(journal),native_outputs=sum(row['candidate_count'] for row in journal),
        host=baseline['environment']['host'],aam_workers=8,slap_workers=1,
        scope='Measured R_to_P mapping CPU only; AAM includes parent and workers and excludes measured persistence/loading. '
              'SLAP sweep includes XYZ preparation, graph construction, mapping and export. '
              'Wall scopes differ: AAM includes archive I/O, SLAP sweep excludes encoding/persistence. '
              'The historical SLAP-sweep scoring timer covers the bidirectional union and cannot be divided by direction.')


def prepare(args):
    args.run.mkdir(parents=True,exist_ok=False)
    for folder in ('status','snapshots','families','enumeration','aam_enumeration','timings'):
        (args.run/folder).mkdir()
    files = ['holdout_forward_analysis.py','holdout_minimum_events.py',
             'holdout_minimum_event_families.py','enumerate_slap_minimum_events.py',
             'enumerate_aam_minimum_events.py']
    for name in files:shutil.copy2(Path(__file__).with_name(name),args.run/name)
    save(args.run/'manifest.json',dict(direction='R_to_P',cases=140,
        config=read(AAM/'manifest.json')['original_config'],aam_source=str(AAM),slap_sweep_source=str(SWEEP),
        prior_alternative_analysis=str(PRIOR),source_commit='6268b0d8bce06d58a5cc82b7d375a2c8af8d270a',
        original_engine_commit=read(AAM/'manifest.json')['original_commit'],
        score_scope='Recorded best forward full-H event score. Alternative catalogues hold this score fixed.',
        timing_scope='Use measured forward-direction CPU and wall records, never half a bidirectional total.',
        scoring_note='Forward SLAP-sweep scoring CPU is not separately available in historical records.',
        sources_sha256={str(p):sha(p) for p in [AAM/'manifest.json',AAM/'case_metrics.json',
            SWEEP/'manifest.json',SWEEP/'case_metrics.json',PRIOR/'manifest.json']},
        frozen_sha256={name:sha(args.run/name) for name in files}))
    print(args.run,flush=True)


def case(args):
    status_path = args.run/f'status/forward_{args.index}.json'
    assert not status_path.exists()
    status = dict(index=args.index,host=socket.gethostname(),started=time.time())
    save(status_path,status)
    started,started_cpu = time.perf_counter(),time.process_time()
    raw = read(AAM/f'inputs/{args.index}/input.json');canonical = EventPatterns(raw)
    baseline = read(AAM/'case_metrics.json')[args.index]
    timing = directional_timing(args.index,baseline)
    minimum = baseline['aam']['directions']['R_to_P']['best_events']
    assert minimum==baseline['aam']['best_events'], 'Do not reuse a bound at another score.'
    tick,cpu = time.perf_counter(),time.process_time()
    graph = read_aam_checkpoint(AAM/f'results/holdout/{args.index}/R_to_P/original/cuts/aam.pkl.gz').graph
    terminals = [t for t in graph.terminals if len(graph.states[t].mapping)==canonical.n]
    patterns={};at_minimum=0
    for offset in range(0,len(terminals),256):
        selected = terminals[offset:offset+256]
        vectors = [[dict(graph.states[t].mapping)[a] for a in range(canonical.n)] for t in selected]
        scores = event_counts(canonical.r,canonical.p,vectors).sum(axis=1)
        assert min(scores,default=minimum)>=minimum
        for terminal,vector,score in zip(selected,vectors,scores,strict=True):
            if score!=minimum:continue
            pattern=canonical.describe(vector);at_minimum+=1
            patterns.setdefault(pattern['id'],dict(pattern,direction='R_to_P',terminal=terminal))
    assert patterns
    timing['forward_terminal_scan'] = dict(cpu=time.process_time()-cpu,wall=time.perf_counter()-tick)
    methods={'aam':dict(minimum=minimum,patterns=patterns,full_terminals_scanned=len(terminals),terminals_at_minimum=at_minimum)}
    del graph
    snapshot=dict(index=args.index,name=raw['name'],methods=methods,directions=['R_to_P'])

    # Native SLAP is already a forward-only call. Reuse its complete catalogue.
    old_enum=read(PRIOR/f'enumeration/{args.index}.json')
    extended=PRIOR/f'enumeration_extended/{args.index}.json'
    if extended.exists():old_enum=read(extended)
    native=old_enum['methods']['native_slap'];assert native['complete']
    methods['native_slap']=dict(minimum=native['minimum'],patterns=native['patterns'])
    enum_methods={'native_slap':native}
    tick,cpu=time.perf_counter(),time.process_time()
    selected=forward_candidates(args.index)
    slap_minimum=min(score['events']['total'] for _,_,score,_ in selected)
    slab_patterns={};families=[];deadline=time.perf_counter()+180
    @lru_cache(None)
    def twin(side,a,b):
        images=list(range(canonical.n));images[a],images[b]=b,a
        return exact_action(images,canonical.features[side])
    for ordinal,candidate,score,origins in selected:
        assert score['hydrogen_score_optimization']['optimal']
        if score['events']['total']!=slap_minimum:continue
        result=enumerate_family(raw,canonical,candidate,score,twin,deadline,check_ms=30000)
        for key,value in result['patterns'].items():
            slab_patterns.setdefault(key,dict(value,candidate=ordinal,sources=[o for o in origins if o['direction']=='R_to_P']))
        families.append(dict(candidate=ordinal,**{k:v for k,v in result.items() if k!='patterns'},pattern_ids=sorted(result['patterns'])))
    enum_methods['slap_sweep']=dict(minimum=slap_minimum,patterns=slab_patterns,families=families,
        complete=all(f['complete'] for f in families),selected_candidates=[j for j,_,_,_ in selected])
    methods['slap_sweep']=dict(minimum=slap_minimum,patterns=slab_patterns)
    timing['forward_slap_event_enumeration']=dict(cpu=time.process_time()-cpu,wall=time.perf_counter()-tick)
    save(args.run/f'enumeration/{args.index}.json',dict(index=args.index,methods=enum_methods,directions=['R_to_P']))
    save(args.run/f'snapshots/{args.index}.json',snapshot)

    # A complete bidirectional catalogue bounds the forward subset at this same K.
    old_aam=read(PRIOR/f'aam_enumeration/{args.index}.json')
    old_family_path=PRIOR/f'families_enriched/{args.index}.json'
    if not old_family_path.exists():old_family_path=PRIOR/f'families/{args.index}.json'
    old_family=read(old_family_path)
    union=dict(old_aam['patterns'])
    for m in methods.values():union.update(m['patterns'])
    known={}
    for key,pattern in union.items():
        if pattern['total']!=minimum:continue
        negative=old_family['methods']['aam'].get(key,{}).get('status')=='excluded_from_saved_families'
        negative|=old_aam['complete'] and key not in old_aam['patterns']
        if negative:
            assert key not in patterns
            known[key]=dict(status='excluded_from_saved_families',method='bidirectional_superset_exclusion')
    generators=tuple(tuple(g[:canonical.n]) for g in pynauty.autgrp(colored_graph([canonical.features[0]]))[0])
    tick,cpu=time.perf_counter(),time.process_time()
    found,metrics=aam_membership(raw,canonical,union,snapshot,generators,time.perf_counter()+args.seconds,
        directions=('R_to_P',),known_results=known)
    timing['forward_aam_membership']=dict(cpu=time.process_time()-cpu,wall=time.perf_counter()-tick)
    family=dict(index=args.index,name=raw['name'],patterns=union,methods={'aam':found},query_metrics={'aam':metrics},directions=['R_to_P'])
    save(args.run/f'families/{args.index}.json',family)
    complete=old_aam['complete'] and all(found[k]['status']!='unresolved' for k in old_aam['patterns'])
    if complete:
        save(args.run/f'aam_enumeration/{args.index}.json',dict(index=args.index,minimum=minimum,
            patterns={k:v['witness'] for k,v in found.items() if v['status']=='represented'},complete=True,
            metrics=dict(proof='Complete bidirectional upper catalogue plus resolved forward membership for every class',
                         upper_patterns=len(old_aam['patterns'])),directions=['R_to_P']))
        save(args.run/f'status/aam_enumeration_{args.index}.json',dict(index=args.index,finished=time.time(),exit=0,proof='subset'))
    else:
        tick,cpu=time.perf_counter(),time.process_time()
        enumerate_aam(SimpleNamespace(run=args.run,index=args.index,seconds=args.seconds,directions=('R_to_P',)))
        timing['forward_aam_event_enumeration']=dict(cpu=time.process_time()-cpu,wall=time.perf_counter()-tick)
    timing['total_new_analysis']=dict(cpu=time.process_time()-started_cpu,wall=time.perf_counter()-started)
    save(args.run/f'timings/{args.index}.json',timing)
    save(status_path,dict(status,finished=time.time(),exit=0))
    print(dict(index=args.index,aam_minimum=minimum,slap_minimum=slap_minimum,seconds=time.perf_counter()-started),flush=True)


def submit(args):
    assert not (args.run/'submission.json').exists()
    command=['sbatch','--parsable','--partition=cpunodes_nia',
        '--exclude=bosque49,bosque50,bosque56,bosque59,bosque60,bosque70','--nodes=1','--cpus-per-task=1','--mem=8G',
        '--time=00:15:00','--no-requeue','--array=0-139%32','--job-name=forward_events',
        f'--output={args.run}/status/forward_%A_%a.out','--wrap',
        shlex.join(['env',f'PYTHONPATH={AAM}/original/src:{AAM}/engine/bench','PYTHONHASHSEED=0','PYTHONDONTWRITEBYTECODE=1',
            'OPENBLAS_NUM_THREADS=1','OMP_NUM_THREADS=1','MKL_NUM_THREADS=1','timeout','--kill-after=5s','800',PYTHON,
            str(args.run/'holdout_forward_analysis.py'),'case','--run',str(args.run),'--seconds',str(args.seconds),'--index'])+' "$SLURM_ARRAY_TASK_ID"']
    job=subprocess.check_output(command,text=True).strip()
    save(args.run/'submission.json',dict(job=job,command=command));print(job,flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command',choices=('prepare','case','submit'))
    parser.add_argument('--run',type=Path,required=True)
    parser.add_argument('--index',type=int)
    parser.add_argument('--seconds',type=int,default=240)
    args=parser.parse_args();globals()[args.command](args)
