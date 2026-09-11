"""Collect, audit, and plot all three completed holdout configurations."""
import argparse
from collections import Counter
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import shutil
import statistics
import subprocess

from holdout_slap_sweep import DIRECTIONS, PAIRED, label_key, read, records, sha


def publish(run,destination):
    from rxn_core import AAMProblem
    from rxn_core.domain import MolecularEndpoint
    from rxn_core.family_scoring import bond_events
    manifest = read(run/'manifest.json')
    baseline = read(PAIRED/'case_metrics.json')
    source = Path(manifest['source'])
    for name,digest in manifest['vendor_sha256'].items():
        assert sha(run/name) == digest
    for name,digest in manifest['drivers_sha256'].items():
        assert sha(run/name) == digest
    for name,digest in manifest['input_hashes'].items():
        assert sha(Path(name) if name.startswith('/') else source/'inputs'/name) == digest
    rows,witnesses,hashes = [],[],{}
    total_calls,total_candidates,certified_candidates = 0,0,0
    for index in range(140):
        status = read(run/f'status/{index}.json')
        assert 'finished' in status and all(p['exit'] == 0 for p in status['phases'].values()), index
        assert status['host'] == baseline[index]['environment']['host']
        sources = read(run/f'evaluations/{index}_sources.json')
        scored = read(run/f'slap_scoring/{index}.refined.json')
        native = read(run/f'slap_xyz/{index}.json')['candidates']
        assert len(native) == len(sources['sources']) == len(scored['slap'])
        raw = read(run/f'inputs/{index}/input.json')
        problem = AAMProblem(*(MolecularEndpoint(**raw[k]) for k in ('reactant','product')))
        n = len(raw['reactant']['elements'])
        method_cpu,mapping_cpu,search_wall,calls,native_count = 0.,0.,0.,0,0
        direction_rows = {}
        all_provenance = set()
        for name in DIRECTIONS:
            folder = run/f'outputs/{index}/{name}'
            header = read(folder/'input.json')
            journal = records(folder/'records.jsonl')
            assert read(folder/'complete.json')['complete']
            assert len(journal) == header['expected_calls'] == len(header['edges'])+1
            assert [r['ordinal'] for r in journal] == list(range(len(journal)))
            assert [r['cut'] for r in journal] == [None,*header['edges']]
            assert all(r['status'] == 'mapped' and r['candidate_count']>0 for r in journal)
            cpu = header['preparation_cpu']+sum(r['graph_cpu']+r['mapping_cpu']+r['export_cpu'] for r in journal)
            wall = header['preparation_wall']+sum(r['graph_wall']+r['mapping_wall']+r['export_wall'] for r in journal)
            method_cpu += cpu
            mapping_cpu += sum(r['mapping_cpu'] for r in journal)
            search_wall += wall
            calls += len(journal)
            native_count += sum(r['candidate_count'] for r in journal)
            for r in journal:
                all_provenance.update((name,r['ordinal'],j) for j in range(r['candidate_count']))
            direction_rows[name] = dict(calls=len(journal),workflow_cpu=cpu,workflow_wall=wall,
                uncut_cpu=sum(journal[0][k] for k in ('graph_cpu','mapping_cpu','export_cpu')))
            for file in ('input.json','complete.json','records.jsonl','native.bin'):
                hashes[str((folder/file).relative_to(run))] = sha(folder/file)
        observed = [(p['direction'],p['call_ordinal'],p['native_candidate'])
                    for entries in sources['sources'] for p in entries]
        assert len(observed) == len(set(observed)) == native_count and set(observed) == all_provenance
        uncut_r,uncut_both,scores = [],[],[]
        uncut_r_keys = set()
        for candidate,origins,prediction in zip(native,sources['sources'],scored['slap'],strict=True):
            mapping = dict(prediction['mapping'])
            assert sorted(mapping) == list(range(n)) and sorted(mapping.values()) == list(range(n))
            assert all(raw['reactant']['elements'][a] == raw['product']['elements'][b] for a,b in mapping.items())
            actual = bond_events(problem,mapping)
            assert actual == prediction['events'], (index,actual,prediction['events'])
            assert prediction['hydrogen_score_optimization']['optimal'], index
            scores.append(actual['total'])
            certified_candidates += 1
            if any(p['cut'] is None for p in origins):
                uncut_both.append(actual['total'])
            if any(p['cut'] is None and p['direction']=='R_to_P' for p in origins):
                uncut_r.append(actual['total'])
                uncut_r_keys.add(label_key(candidate))
            # Retain all candidates for independent per-case score inspection.
            witnesses.append(dict(index=index,candidate=prediction['candidate'],mapping=[mapping[i] for i in range(n)],
                events=actual,sources=origins))
        prior = read(PAIRED/f'slap_xyz/{index}.json')['candidates']
        assert uncut_r_keys == {label_key(c) for c in prior}, index
        assert min(uncut_r) == baseline[index]['slap']['best_events'], index
        best = min(scores)
        assert best == scored['slap_min_representative_events'] and best <= min(uncut_both)
        total_calls += calls
        total_candidates += native_count
        aam = baseline[index]['aam']
        winner = 'aam' if aam['best_events']<best else 'slap_sweep' if aam['best_events']>best else 'tie'
        rows.append(dict(index=index,name=baseline[index]['name'],host=status['host'],
            aam_events=aam['best_events'],native_slap_events=min(uncut_r),
            bidirectional_uncut_slap_events=min(uncut_both),slap_sweep_events=best,
            winner=winner,aam_cpu=aam['cpu_seconds'],native_slap_cpu=baseline[index]['slap']['cpu_seconds'],
            slap_sweep_workflow_cpu=method_cpu,slap_sweep_mapping_cpu=mapping_cpu,
            slap_sweep_workflow_wall=search_wall,
            aam_scoring_cpu=aam['scoring_cpu_seconds'],native_slap_scoring_cpu=baseline[index]['slap']['scoring_cpu_seconds'],
            slap_sweep_scoring_cpu=scored['scoring_cpu_seconds'],
            aggregation_cpu_including_read_decode=sources['aggregation_cpu_including_read_decode'],
            calls=calls,native_candidates=native_count,unique_label_partitions=len(native),directions=direction_rows))
        for name in (f'slap_xyz/{index}.json',f'slap_scoring/{index}.refined.json',f'evaluations/{index}_sources.json',f'status/{index}.json'):
            hashes[name] = sha(run/name)
    timing = {}
    for label,key,score_key in [('aam','aam_cpu','aam_scoring_cpu'),('native_slap','native_slap_cpu','native_slap_scoring_cpu'),
                              ('slap_sweep','slap_sweep_workflow_cpu','slap_sweep_scoring_cpu')]:
        times = [r[key] for r in rows]
        timing[label] = dict(total_cpu_seconds=sum(times),mean_cpu_seconds=statistics.mean(times),
            median_cpu_seconds=statistics.median(times),max_cpu_seconds=max(times),
            scoring_cpu_seconds=sum(r[score_key] for r in rows))
    summary = dict(cases=140,paired_complete_cases=140,full_mapping_cases=dict(aam=140,native_slap=140,slap_sweep=140),
        aam_vs_native=dict(Counter('aam' if r['aam_events']<r['native_slap_events'] else
                         'native_slap' if r['aam_events']>r['native_slap_events'] else 'tie' for r in rows)),
        aam_vs_sweep=dict(Counter(r['winner'] for r in rows)),
        sweep_vs_native=dict(Counter('sweep_lower' if r['slap_sweep_events']<r['native_slap_events'] else 'tie' for r in rows)),
        sweep_vs_bidirectional_uncut=dict(Counter('sweep_lower' if r['slap_sweep_events']<r['bidirectional_uncut_slap_events'] else 'tie' for r in rows)),
        timing=timing,aam_cpu_over_slap_sweep=timing['aam']['total_cpu_seconds']/timing['slap_sweep']['total_cpu_seconds'],
        sweep_cpu_over_native=timing['slap_sweep']['total_cpu_seconds']/timing['native_slap']['total_cpu_seconds'],
        total_sweep_calls=total_calls,total_native_candidates=total_candidates,
        unique_scored_label_partitions=certified_candidates,unresolved_scores=0,search_errors=0,unfinished_variants=0,
        matched_native_uncut_controls=140,matched_hosts=140,chemical_accuracy=None,
        graph_scope=manifest['graph_protocol'],cut_scope=manifest['cut_protocol'],timing_scope=manifest['timing'])
    for name,value in [('summary.json',summary),('case_metrics.json',rows)]:
        (run/name).write_text(json.dumps(value,indent=2)+'\n')
    destination.mkdir(parents=True,exist_ok=False)
    for name in ('manifest.json','submissions.json','summary.json','case_metrics.json'):
        shutil.copy2(run/name,destination/name)
    with (destination/'per_case.csv').open('w',newline='') as stream:
        keys = [k for k in rows[0] if k!='directions']
        writer = csv.DictWriter(stream,fieldnames=keys,lineterminator='\n')
        writer.writeheader()
        writer.writerows({k:r[k] for k in keys} for r in rows)
    with gzip.open(destination/'audited_witnesses.json.gz','wt') as stream:
        json.dump(witnesses,stream)
    for label,paths in [('native_call_records',[(i,d,run/f'outputs/{i}/{d}/records.jsonl') for i in range(140) for d in DIRECTIONS]),
                        ('scored_candidates',[(i,'scored',run/f'slap_scoring/{i}.refined.json') for i in range(140)])]:
        with gzip.open(destination/f'{label}.json.gz','wt') as stream:
            json.dump({f'{i}/{d}':records(p) if label=='native_call_records' else read(p) for i,d,p in paths},stream)
    (destination/'archive_hashes.json').write_text(json.dumps(hashes,indent=2)+'\n')
    jobs = ','.join(s['job'].split(';')[0] for s in read(run/'submissions.json'))
    (destination/'slurm_accounting.psv').write_text(subprocess.check_output(['sacct','-j',jobs,'--array','--noheader','--parsable2',
        '--format=JobID,State,ExitCode,ElapsedRaw,TotalCPU,AllocCPUS,MaxRSS,NodeList'],text=True))
    report(rows,summary,run,destination)
    figure(rows,summary,destination)
    print(json.dumps(dict(report=str(destination),summary=summary),indent=2),flush=True)


def report(rows,summary,run,destination):
    timing = summary['timing']
    comparison = summary['aam_vs_sweep']
    timings = '\n'.join(f"| {label} | {timing[key]['total_cpu_seconds']:.3f} | {timing[key]['mean_cpu_seconds']:.3f} | {timing[key]['median_cpu_seconds']:.3f} | {timing[key]['scoring_cpu_seconds']:.3f} |"
        for label,key in [('AAM: one seed, cap 1000','aam'),('Native XYZ SLAP','native_slap'),('SLAP + single-edge sweep','slap_sweep')])
    different = '\n'.join(f"| {r['index']} | {r['name']} | {r['aam_events']} | {r['native_slap_events']} | {r['bidirectional_uncut_slap_events']} | {r['slap_sweep_events']} |"
        for r in rows if len({r['aam_events'],r['native_slap_events'],r['bidirectional_uncut_slap_events'],r['slap_sweep_events']})>1)
    (destination/'README.md').write_text(f"""# 140 cases: one-seed cap-1000 AAM versus native SLAP and SLAP with cuts

All three configurations return full mappings on all 140 XYZ/WBO cases.
Against native XYZ SLAP, AAM has six lower minimum-change scores, 134 ties,
and no higher scores. Against **SLAP plus the single-edge sweep**, AAM has
**{comparison.get('aam',0)} lower scores, {comparison.get('tie',0)} ties,
and {comparison.get('slap_sweep',0)} higher scores**.

The sweep improves on native SLAP in {summary['sweep_vs_native'].get('sweep_lower',0)} cases.
Compared with the bidirectional uncut SLAP control, cuts improve
{summary['sweep_vs_bidirectional_uncut'].get('sweep_lower',0)} cases.

| Configuration | Total mapping CPU s | Mean CPU s/case | Median CPU s/case | Separate scoring CPU s |
| --- | ---: | ---: | ---: | ---: |
{timings}

AAM uses {summary['aam_cpu_over_slap_sweep']:.3f} times the recorded mapping CPU
of SLAP plus the sweep. The sweep uses {summary['sweep_cpu_over_native']:.3f}
times the native SLAP mapping CPU. Per-case timings are in `per_case.csv`.

## Protocol and timing scope

AAM is the unchanged stable engine `98b01b175eeed31f70d13e7cbf178b80bf07c9e0`,
with one deterministic seed order per cut, cap 1000, tolerance 1.0, both directions,
and the existing uncut/single-edge sweep. Its results and native SLAP timings are
from the fresh paired experiment in `../holdout_cap1000_seed1_20260910/`.

The added SLAP sweep keeps native XYZ graph perception (`geoms2lgp`, original
component files, bond scale 1.2) and `binary=True`, with heavy-atom symmetry
breaking. It runs the unmodified mapper on the uncut graph and every individual
source-edge deletion, including H bonds, in both directions. Initial WL labels
are rebuilt after each cut. Native XYZ edges are binary; this is a cut ablation
on the existing XYZ comparator, separate from the Golden SMILES binary/weighted
union. Every output is scored on the original full-H WBO matrices, with graph
threshold 0.2 and bond-order-change threshold 0.5. No reference mappings or
previous scores enter the search. These are best saved representative scores,
not reference accuracy or a proof of global optimality.

Each sweep case ran on the same host as the paired AAM/native-SLAP case. AAM
uses eight workers, while SLAP uses one; the table reports additive CPU, not
equal-core wall latency. AAM CPU excludes measured checkpoint persistence and
loading. Sweep CPU includes initial XYZ graph preparation, per-cut graph
construction, native mapping and label export; compressed archival I/O and
offline scoring are separate. Scoring workloads differ: AAM also canonicalizes
heavy mapping classes, while SLAP optimizes H assignments within its native
label groups. Do not interpret the mapping-only ratio as full application cost.

## Audit

All {summary['total_sweep_calls']:,} planned SLAP mapping calls completed without
errors or timeouts. All {summary['total_native_candidates']:,} native outputs
are accounted for through {summary['unique_scored_label_partitions']:,} unique
label partitions. Identical partitions are deduplicated across cuts/directions;
all discovering cut identities are retained. Every scored representative is a
full element-preserving bijection and was independently rescored with the scalar
bond-event implementation. Every H-label optimization certified its final score.
All 140 uncut R-to-P controls reproduce the exact native label partitions and
minimum scores of the preceding comparison.

`audited_witnesses.json.gz` contains every scored mapping and its discovering
cuts. `native_call_records.json.gz` records all committed cuts, timings and
native-frame hashes. Full native LAP archives remain at `{run}`;
`archive_hashes.json` identifies them. Frozen sources and input hashes are in
`manifest.json`. There are no unresolved calls or score optimizations.
Case 17 had one allocation stall before starting; it was canceled and submitted
again on the same host, with no mapping work repeated. Both allocation records
are retained in `submissions.json` and the Slurm accounting export.

## Cases with differing scores

Indices are zero-based, matching the earlier reports.

| Case | Name | AAM | Native SLAP | Uncut SLAP both directions | SLAP sweep |
| --- | --- | ---: | ---: | ---: | ---: |
{different}
""")


def figure(rows,summary,destination):
    os.environ.setdefault('MPLCONFIGDIR',str(destination/'.mpl-cache'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'pdf.fonttype':42,
        'axes.spines.top':False,'axes.spines.right':False})
    fig,axes = plt.subplots(1,3,figsize=(11,3.8),layout='constrained')
    points = Counter((r['slap_sweep_events'],r['aam_events']) for r in rows)
    limit = max(max(p) for p in points)+2
    axes[0].plot([0,limit],[0,limit],ls='--',color='#ABB6BA',lw=1)
    for (x,y),n in points.items():
        axes[0].scatter(x,y,s=20+16*np.sqrt(n),alpha=.75,color='#007F68' if y<=x else '#CC6B26')
    axes[0].set(xlabel='SLAP sweep: minimum changes',ylabel='AAM: minimum changes',
        xlim=(-.6,limit),ylim=(-.6,limit),title='a  Best saved scores')
    comparison = summary['aam_vs_sweep']
    values = [comparison.get('aam',0),comparison.get('tie',0),comparison.get('slap_sweep',0)]
    axes[1].barh([2,1,0],values,color=['#007F68','#ABB6BA','#CC6B26'])
    for y,n in zip([2,1,0],values):axes[1].text(n+2,y,str(n),va='center')
    axes[1].set(yticks=[2,1,0],yticklabels=['AAM lower','Tie','SLAP sweep lower'],xlim=(0,150),
        xlabel='Cases',title='b  Minimum-change comparison')
    labels = ['AAM\n1 seed, cap 1000','Native\nSLAP','SLAP\n+ sweep']
    times = [summary['timing'][k]['mean_cpu_seconds'] for k in ('aam','native_slap','slap_sweep')]
    axes[2].bar(range(3),times,color=['#007F68','#ABB6BA','#3269A8'])
    for x,n in enumerate(times):axes[2].text(x,n*1.15,f'{n:.3f}',ha='center',fontsize=9)
    axes[2].set(xticks=range(3),xticklabels=labels,yscale='log',ylim=(min(times)*.5,max(times)*2.2),
        ylabel='Mean mapping CPU seconds / case (log)',title='c  All 140 paired cases')
    for ax in axes:
        ax.grid(axis='y',alpha=.12)
        ax.set_axisbelow(True)
    fig.get_layout_engine().set(rect=(0,.15,1,.85))
    fig.text(.02,.07,'Native XYZ single-edge sweep, both directions | AAM: frozen stable engine, one seed, cap 1000',fontsize=9)
    fig.text(.02,.02,'Full-H WBO event scores; no reference accuracy labels | Mapping CPU excludes offline scoring and archival I/O',fontsize=8)
    for ext in ('png','pdf','svg'):fig.savefig(destination/f'comparison.{ext}',dpi=200,facecolor='white')
    svg = destination/'comparison.svg'
    svg.write_text('\n'.join(line.rstrip() for line in svg.read_text().splitlines())+'\n')
    plt.close(fig)
    shutil.rmtree(destination/'.mpl-cache',ignore_errors=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run',type=Path,required=True)
    parser.add_argument('--destination',type=Path,required=True)
    args = parser.parse_args()
    publish(args.run,args.destination)
