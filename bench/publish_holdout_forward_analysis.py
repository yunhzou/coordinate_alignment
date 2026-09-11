"""Publish the audited forward-only elementary-step comparison."""
import argparse
from collections import Counter
import csv
import gzip
import json
import os
from pathlib import Path
import shutil
import statistics
import tarfile

from holdout_minimum_events import read,save,sha,AAM,SWEEP
from holdout_forward_analysis import PRIOR
from publish_minimum_event_analysis import collect,METHODS


def publish(args):
    rows,alternatives,hashes=collect(args.run)
    baseline=read(AAM/'case_metrics.json')
    timing_rows=[]
    for row in rows:
        index=row['index'];status=read(args.run/f'status/forward_{index}.json')
        assert status['exit']==0 and 'finished' in status
        timing=read(args.run/f'timings/{index}.json');timing_rows.append(timing)
        row['timing']=timing;row['directions']=['R_to_P']
        assert row['minima']['aam']==baseline[index]['aam']['best_events']
        snapshot=read(args.run/f'snapshots/{index}.json')
        assert all(p.get('direction')=='R_to_P' for p in snapshot['methods']['aam']['patterns'].values())
        families=read(args.run/f'families/{index}.json')['methods']['aam']
        for value in families.values():
            if value['status']=='represented':
                assert value.get('direction',value['witness'].get('direction'))=='R_to_P'
        aam_enum=read(args.run/f'aam_enumeration/{index}.json')
        for key,value in aam_enum['patterns'].items():
            if families.get(key,{}).get('status')!='represented':assert value['direction']=='R_to_P'
        origins=read(SWEEP/f'evaluations/{index}_sources.json')['sources']
        journal=[json.loads(line) for line in (SWEEP/f'outputs/{index}/R_to_P/records.jsonl').read_text().splitlines()]
        actual=[(o['call_ordinal'],o['native_candidate']) for group in origins for o in group if o['direction']=='R_to_P']
        expected={(j['ordinal'],c) for j in journal for c in range(j['candidate_count'])}
        assert len(actual)==len(set(actual)) and set(actual)==expected
        candidates=read(SWEEP/f'slap_xyz/{index}.json')['candidates']
        enum=read(args.run/f'enumeration/{index}.json')['methods']['slap_sweep']
        for pattern in enum['patterns'].values():
            assert pattern['sources'] and all(s['direction']=='R_to_P' for s in pattern['sources'])
            graphs=candidates[pattern['candidate']]['graphs']
            assert all(graphs[0]['labels'][a]==graphs[1]['labels'][b] for a,b in enumerate(pattern['mapping']))
        for p in [args.run/f'timings/{index}.json',args.run/f'status/forward_{index}.json',
                  SWEEP/f'evaluations/{index}_sources.json',SWEEP/f'outputs/{index}/R_to_P/records.jsonl',
                  AAM/f'results/holdout/{index}/R_to_P/original/search.json']:
            hashes[str(p)]=sha(p)
    scores={m:dict(Counter('aam_lower' if r['minima']['aam']<r['minima'][m] else
        'slap_lower' if r['minima']['aam']>r['minima'][m] else 'tie' for r in rows)) for m in METHODS[1:]}
    timings={}
    for method in METHODS:
        result={}
        for metric in ('cpu','wall'):
            values=[r[f'{method}_{metric}'] for r in timing_rows]
            result[metric]=dict(total=sum(values),mean=statistics.mean(values),median=statistics.median(values),maximum=max(values))
        result['historical_scoring_cpu']=sum(r[f'{method}_scoring_cpu'] for r in timing_rows) if method!='slap_sweep' else None
        timings[method]=result
    phases=['forward_terminal_scan','forward_slap_event_enumeration','forward_aam_membership',
            'forward_aam_event_enumeration','total_new_analysis']
    audit_timing={phase:{unit:sum(r.get(phase,{}).get(unit,0) for r in timing_rows) for unit in ('cpu','wall')} for phase in phases}
    summary=dict(direction='R_to_P',cases=140,full_mappings={m:140 for m in METHODS},score_comparisons=scores,
        timing=timings,aam_cpu_over_slap_sweep=timings['aam']['cpu']['total']/timings['slap_sweep']['cpu']['total'],
        calls=sum(r['calls'] for r in timing_rows),native_outputs=sum(r['native_outputs'] for r in timing_rows),
        alternatives=alternatives,new_analysis_timing=audit_timing,
        chemical_accuracy=None,scope='Forward-only, one seed, cap 1000; recorded minimum full-H event scores and mapping CPU. '
        'Alternative counts merge joint reactant symmetries; incomplete catalogues remain lower bounds. '
        'Elapsed times use eight AAM workers versus one SLAP worker and different persistence scopes.')
    destination=args.destination;destination.mkdir(parents=True,exist_ok=False)
    for name,value in [('summary.json',summary),('per_case.json',rows),('analysis_hashes.json',hashes)]:save(destination/name,value)
    with (destination/'per_case.csv').open('w',newline='') as stream:
        writer=csv.writer(stream,lineterminator='\n')
        writer.writerow(['index','name','aam_events','native_slap_events','slap_sweep_events',
            'aam_cpu','native_slap_cpu','slap_sweep_cpu','aam_wall','native_slap_wall','slap_sweep_wall',
            'aam_patterns','native_slap_patterns','slap_sweep_patterns','aam_complete','sweep_score_tied',
            'shared_patterns','aam_only_patterns','slap_only_proven_patterns','unresolved_patterns','host'])
        for r in rows:
            c=r['comparisons']['slap_sweep'];t=r['timing']
            writer.writerow([r['index'],r['name'],*(r['minima'][m] for m in METHODS),
                *(t[f'{m}_{unit}'] for unit in ('cpu','wall') for m in METHODS),
                *(len(r['pattern_ids'][m]) for m in METHODS),r['catalogue_complete']['aam'],c['comparable'],
                *(len(c[k]) for k in ('shared','aam_only','slap_only_proven','unresolved')),t['host']])
    for folder in ('snapshots','families','enumeration','aam_enumeration','timings','status'):
        with gzip.open(destination/f'{folder}.json.gz','wt') as stream:
            json.dump({p.stem:read(p) for p in sorted((args.run/folder).glob('*.json'))},stream)
    for name in ['manifest.json','submission.json','slurm_accounting.tsv']:
        shutil.copy2(args.run/name,destination/name)
    shutil.copy2(args.run/'validation/direct_checks.json',destination/'directional_checks.json')
    for p in args.run.glob('*relocation*.json'):shutil.copy2(p,destination/p.name)
    provenance=dict(frozen_sha256={p.name:sha(p) for p in args.run.glob('*.py')},
        published_sha256={p.name:sha(p) for p in [Path(__file__),Path(__file__).with_name('publish_minimum_event_analysis.py'),
            Path(__file__).with_name('minimum_events_viewer.html')]},
        prior_catalogue_hashes={str(p):sha(p) for p in PRIOR.glob('aam_enumeration/*.json')},
        prior_checks=str(PRIOR/'validation/finite_checks.json'))
    save(destination/'provenance.json',provenance)
    with tarfile.open(destination/'analysis_sources.tar.gz','w:gz') as archive:
        for p in args.run.glob('*.py'):archive.add(p,arcname=f'frozen/{p.name}')
        for p in [Path(__file__),Path(__file__).with_name('publish_minimum_event_analysis.py'),Path(__file__).with_name('minimum_events_viewer.html')]:
            archive.add(p,arcname=f'published/{p.name}')
        archive.add(PRIOR/'validation/finite_checks.json',arcname='validation/finite_checks.json')
        archive.add(args.run/'validation/direct_checks.json',arcname='validation/direct_checks.json')
    report(rows,summary,destination,args.run)
    figure(rows,summary,destination)
    template=Path(__file__).with_name('minimum_events_viewer.html').read_text()
    template=template.replace('The same total can describe different bond changes.','Reactant → product only; one seed, cap 1000. The same total can describe different bond changes.')
    template=template.replace('every individual edge deletion in both directions','every individual reactant-edge deletion in the reactant → product direction')
    template=template.replace('The analysis scans all saved full AAM terminal witnesses','The analysis scans forward AAM terminal witnesses')
    payload=json.dumps(dict(cases=rows,summary=alternatives),separators=(',',':')).replace('</','<\\/')
    (destination/'viewer.html').write_text(template.replace('__PAYLOAD__',payload))
    print(json.dumps(summary,indent=2),flush=True)


def report(rows,s,destination,run):
    alt=s['alternatives'];a=alt['comparisons']['slap_sweep'];n=alt['comparisons']['native_slap']
    names={'aam':'AAM + sweep','native_slap':'Native SLAP','slap_sweep':'SLAP + sweep'}
    timing='\n'.join(f"| {names[m]} | {s['timing'][m]['cpu']['total']:.3f} | {s['timing'][m]['cpu']['mean']:.3f} | {s['timing'][m]['cpu']['median']:.3f} | {s['timing'][m]['wall']['mean']:.3f} |" for m in METHODS)
    differences='\n'.join(f"| {r['index']} | {r['name']} | {r['minima']['aam']} | {r['minima']['native_slap']} | {r['minima']['slap_sweep']} |" for r in rows if len(set(r['minima'].values()))>1)
    examples=sorted(set([11,64,101,129]+a['difference_cases']))
    patterns='\n'.join(f"| {i} | {rows[i]['minima']['aam']} / {rows[i]['minima']['slap_sweep']} | {'≥' if not rows[i]['catalogue_complete']['aam'] else ''}{len(rows[i]['pattern_ids']['aam'])} | {len(rows[i]['pattern_ids']['slap_sweep'])} | {len(rows[i]['comparisons']['slap_sweep']['shared'])} | {len(rows[i]['comparisons']['slap_sweep']['aam_only'])} / {len(rows[i]['comparisons']['slap_sweep']['slap_only_proven'])} |" for i in examples)
    (destination/'README.md').write_text(f"""# Forward-only comparison on 140 elementary-step cases

**Reactant → product only.** AAM uses one seed order per cut, branch cap 1000,
tolerance 1.0 and the uncut/single-edge sweep. SLAP is shown as the native XYZ
call and as the uncut plus every individual reactant-edge deletion, including H.
The frozen engine is `{read(run/'manifest.json')['original_engine_commit']}`.

All three protocols provide full mappings for all 140 cases. AAM has
{s['score_comparisons']['native_slap'].get('aam_lower',0)} lower scores and
{s['score_comparisons']['native_slap'].get('tie',0)} ties against native SLAP;
against SLAP with the sweep it has
**{s['score_comparisons']['slap_sweep'].get('aam_lower',0)} lower scores,
{s['score_comparisons']['slap_sweep'].get('tie',0)} ties and
{s['score_comparisons']['slap_sweep'].get('slap_lower',0)} higher scores**.
No reference mappings are available, so these are minimum-change scores rather
than mapping-accuracy percentages.

## Mapping time

| Method | Total CPU s | Mean CPU s/case | Median CPU s/case | Mean observed wall s/case |
| --- | ---: | ---: | ---: | ---: |
{timing}

AAM uses **{s['aam_cpu_over_slap_sweep']:.3f}×** the mapping CPU of SLAP+sweep.
These are actual saved forward-direction measurements, never halved totals.
Each comparison uses the same host per case. AAM CPU sums the parent and eight
workers and excludes measured checkpoint persistence/loading. SLAP sweep CPU
includes XYZ preparation, graph construction, native mapping and export; it
excludes frame encoding, persistence and subsequent candidate aggregation.
The forward sweep contains {s['calls']:,} native calls and {s['native_outputs']:,} native outputs.

**Wall time is an observed execution comparison with different core counts:**
AAM used eight workers; SLAP used one. AAM wall includes archive I/O; SLAP-sweep
wall is the sum of measured preparation/construction/mapping/export durations.
It is not an equal-core speedup or full end-to-end latency comparison.

Ranking/scoring is separate from mapping: the saved forward AAM scorer used
{s['timing']['aam']['historical_scoring_cpu']:.3f} CPU s and native SLAP scoring used
{s['timing']['native_slap']['historical_scoring_cpu']:.3f} CPU s. The saved sweep
scoring timer covers the deduplicated bidirectional candidate union. It cannot
be assigned to the forward subset, so forward SLAP-sweep scoring CPU is null.
No total mapping-plus-scoring speedup is claimed.

## Alternative minimum-event patterns

At the {a['tied_cases']} AAM / SLAP-sweep tied scores, the catalogues contain
{a['verified_aam_patterns']} verified AAM patterns and {a['slap_patterns']} SLAP-sweep
patterns: **{a['shared_patterns']} shared, {a['aam_only_patterns']} AAM-only,
{a['slap_only_proven_patterns']} SLAP-only proved, and {a['unresolved_memberships']}
unresolved memberships**. Both methods have multiple verified alternatives in
{a['both_multiple_cases']} tied cases. Observed sets agree in
{a['equal_observed_catalogue_cases']} cases; equality is certified for
{a['complete_equal_catalogue_cases']} complete catalogues.

| Case | AAM / sweep score | Verified AAM patterns | Sweep patterns | Shared | Proved AAM-only / sweep-only |
| --- | ---: | ---: | ---: | ---: | ---: |
{patterns}

SLAP catalogues are complete for all 140 saved forward output families. AAM
catalogues are complete in **{alt['complete_catalogue_cases']['aam']}/140 cases**;
remaining counts are lower bounds. A complete bidirectional catalogue can bound
the forward subset at the same score: every class is then tested for forward
membership. Reverse witnesses are never accepted as forward detections. Other
cases use bounded enumeration of the forward AAM families. An incomplete scan
never proves exclusion. Comparisons of alternatives are inapplicable when the
two recorded minimum scores differ.

Against native SLAP, {n['tied_cases']} scores tie, with {n['shared_patterns']} shared
patterns, {n['aam_only_patterns']} AAM-only patterns and
{n['slap_only_proven_patterns']} SLAP-only patterns proved.

Patterns distinguish broken, formed, strengthened and weakened bonds, including
explicit H, modulo joint symmetries of the reactant event-response graph. Bond
presence uses WBO > 0.2 and order changes use magnitude > 0.5. These counts describe
candidate bond edits; they do not validate a physical reaction pathway.

## Effect of dropping the reverse direction

AAM retains its bidirectional best score on every case. SLAP sweep loses the
better scores contributed by reverse searches in cases 31 and 96. Consequently,
the earlier two AAM-lower scores become four in the forward-only comparison.

| Case | Name | AAM | Native SLAP | SLAP sweep |
| --- | --- | ---: | ---: | ---: |
{differences}

## Verification and artifacts

The audit scans {alt['full_aam_terminals_scanned']:,} forward full-terminal AAM
witnesses and independently rescored {alt['independently_rescored_witnesses']}
published method/pattern witnesses with the scalar bond-event routine. The
event-enumeration engine was also checked against 13,872 exhaustive assignments
in ten small real SLAP families in the preceding audit. Candidate origins and
all forward call journals are checked before accepting scores or times.
Direct forward-family enumeration independently confirms the subset-catalogue
results for cases 11, 64 and 101 (`directional_checks.json`).

[Open the offline bond-edit viewer](viewer.html). `comparison.pdf`, `.png` and
`.svg` show the scores, mapping CPU and alternative overlap. `per_case.csv` and
`per_case.json` contain all 140 records. Compressed case outputs, source hashes,
frozen analysis scripts, job records and the preceding finite checks accompany
the report. Full source run: `{run}`.

`summary.json:new_analysis_timing` records case-worker analysis separately
({s['new_analysis_timing']['total_new_analysis']['cpu']:.3f} CPU s; publication,
validation and process startup are excluded). It reuses
prior catalogues and exclusion proofs, so it is not a fresh alternative-extraction
speed benchmark. Mapping searches were not rerun, and the earlier bidirectional
report is preserved as a distinct protocol.
""")


def figure(rows,s,destination):
    os.environ.setdefault('MPLCONFIGDIR',str(destination/'.mpl-cache'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'pdf.fonttype':42,'svg.fonttype':'none',
        'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(1,3,figsize=(12,4.2),layout='constrained')
    for j,m in enumerate(METHODS[1:]):
        c=s['score_comparisons'][m];offset=0
        for key,color in [('aam_lower','#007F68'),('tie','#C2CDD1'),('slap_lower','#C36D30')]:
            count=c.get(key,0);axes[0].barh(j,count,left=offset,color=color)
            if count:axes[0].text(offset+count/2,j,str(count),ha='center',va='center',fontsize=8)
            offset+=count
    axes[0].set(yticks=[0,1],yticklabels=['Native SLAP','SLAP + sweep'],xlabel='Cases',title='a  AAM minimum-event scores')
    axes[0].set_ylim(-.5,1.95)
    axes[0].text(.02,.94,'Green: AAM lower · Gray: tie',transform=axes[0].transAxes,fontsize=8,va='top')
    values=[s['timing'][m]['cpu']['mean'] for m in METHODS]
    axes[1].barh(range(3),values,color=['#007F68','#8599A2','#3269A8'])
    for j,v in enumerate(values):axes[1].text(v+.08,j,f'{v:.3f}',va='center',fontsize=8)
    axes[1].set(yticks=range(3),yticklabels=['AAM + sweep','Native SLAP','SLAP + sweep'],
        xlim=(0,max(values)*1.18),xlabel='Mean mapping CPU seconds / case',title='b  Measured forward mapping CPU')
    axes[1].invert_yaxis()
    a=s['alternatives']['comparisons']['slap_sweep']
    counts=[a[k] for k in ['shared_patterns','aam_only_patterns','slap_only_proven_patterns','unresolved_memberships']]
    axes[2].barh(range(4),counts,color=['#007F68','#3269A8','#C36D30','#C2CDD1'])
    for j,v in enumerate(counts):axes[2].text(v+1,j,str(v),va='center',fontsize=8)
    axes[2].set(yticks=range(4),yticklabels=['Shared','AAM only','Sweep only','Unresolved'],
        xlim=(0,max(counts)*1.16),xlabel='Distinct case / bond-edit patterns',title=f"c  Alternatives on {a['tied_cases']} ties")
    axes[2].invert_yaxis()
    fig.get_layout_engine().set(rect=(0,.2,1,.8))
    fig.text(.02,.095,'Reactant → product only | AAM: stable engine, 1 seed, cap 1000 | Full-H WBO event scoring',fontsize=9)
    fig.text(.02,.035,f"AAM alternatives complete in {s['alternatives']['complete_catalogue_cases']['aam']}/140 cases; other counts are lower bounds. Mapping CPU excludes scoring.",fontsize=8)
    for ext in ('pdf','png','svg'):fig.savefig(destination/f'comparison.{ext}',dpi=200,facecolor='white')
    p=destination/'comparison.svg';p.write_text('\n'.join(line.rstrip() for line in p.read_text().splitlines())+'\n')
    plt.close(fig);shutil.rmtree(destination/'.mpl-cache',ignore_errors=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run',type=Path,required=True)
    parser.add_argument('--destination',type=Path,required=True)
    publish(parser.parse_args())
