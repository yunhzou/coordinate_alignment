"""Audit and publish the completed cap/seed experiment without rerunning search."""
import argparse
from collections import Counter, defaultdict
import gzip
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess


def read(path):
    return json.loads(path.read_text())


def sha(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def publish(run, destination):
    from rxn_core import AAMProblem
    from rxn_core.domain import MolecularEndpoint
    from rxn_core.family_scoring import bond_events
    from rxn_core.artifacts import read_aam_checkpoint
    manifest = read(run/'manifest.json')
    rows, summary = read(run/'case_metrics.json'), read(run/'summary.json')
    assert len(rows) == summary['paired_complete_cases'] == 140
    assert summary['phase_exits'] == {'0':840}
    assert summary['full_mapping_cases'] == {'aam':140,'slap':140}
    assert summary['unresolved_slap_h_candidates'] == 0
    assert {k for k in manifest['original_config']
            if manifest['original_config'][k] != manifest['baseline_config'][k]} == {'seed_count','branch_limit'}
    assert manifest['original_config']['seed_count'] == 1 and manifest['original_config']['branch_limit'] == 1000
    for name, digest in manifest['frozen_sha256'].items():
        assert sha(run/name) == digest, name
    assert sha(run/'driver.py') == manifest['driver_sha256']
    package = Path(manifest['slap_python']).parent.parent/'lib/python3.10/site-packages/slapmapper'
    for name, digest in manifest['slap_source_sha256'].items():
        assert sha(package/name) == digest, name
    source = Path(manifest['source'])
    for name, digest in read(run/'input_hashes.json').items():
        path = Path(name) if name.startswith('/') else source/'inputs'/name
        assert sha(path) == digest, name
    witness_rows, hashes, context_checks = [], {}, []
    for row in rows:
        index = row['index']
        raw = read(run/f'inputs/{index}/input.json')
        problem = AAMProblem(*(MolecularEndpoint(**raw[k]) for k in ('reactant','product')))
        n = len(raw['reactant']['elements'])
        def audit(mapping, counts):
            assert sorted(mapping) == list(range(n)) and sorted(mapping.values()) == list(range(n))
            assert all(raw['reactant']['elements'][a] == raw['product']['elements'][b] for a,b in mapping.items())
            actual = bond_events(problem, mapping)
            assert actual == counts, (index,actual,counts)
            return actual
        for direction, metric in row['aam']['directions'].items():
            folder = run/f'results/holdout/{index}/{direction}/original'
            if index in (0,123):
                archive = read_aam_checkpoint(folder/'cuts/aam.pkl.gz')
                by_cut = defaultdict(set)
                for context in archive.graph.contexts:
                    assert context.branch_limit == 1000 and context.iso_tolerance == 1.0
                    by_cut[str(context.cuts)].add(context.seed_order)
                assert all(len(orders) == 1 for orders in by_cut.values())
                context_checks.append(dict(index=index,direction=direction,
                    contexts=len(archive.graph.contexts),cuts=len(by_cut),
                    branch_limit=1000,max_seed_orders_per_cut=max(map(len,by_cut.values()),default=0)))
            classes = read(folder/'full_sweep_classes.json')
            assert min((c['events'] for c in classes),default=None) == metric['best_events']
            if classes:
                best = min(classes, key=lambda c:c['events'])
                counts = dict(zip(('broken','formed','order_changed'),best['counts']))
                counts['total'] = best['events']
                audit(dict(enumerate(best['mapping'])),counts)
                witness_rows.append(dict(index=index, method='aam', direction=direction,
                    mapping=best['mapping'], events=counts, family_id=best['id']))
            for name in ('search.json','environment.json','full_sweep_evaluation.json','full_sweep_classes.json','cuts/aam.pkl.gz'):
                path = folder/name
                hashes[str(path.relative_to(run))] = sha(path)
        scored = read(run/f'slap_scoring/{index}.refined.json')
        assert min(c['events']['total'] for c in scored['slap']) == row['slap']['best_events']
        for candidate in scored['slap']:
            mapping = dict(candidate['mapping'])
            audit(mapping,candidate['events'])
            witness_rows.append(dict(index=index, method='slap', candidate=candidate['candidate'],
                mapping=[mapping[i] for i in range(n)], events=candidate['events']))
        for name in (f'slap_xyz/{index}.json',f'slap_scoring/{index}.refined.json',f'status/{index}.json'):
            hashes[name] = sha(run/name)
    destination.mkdir(parents=True, exist_ok=False)
    for name in ('manifest.json','submission.json','summary.json','case_metrics.json','per_case.csv','input_hashes.json'):
        shutil.copy2(run/name,destination/name)
    csv_path = destination/'per_case.csv'
    csv_path.write_text(csv_path.read_text())
    for name in ('slap_xyz','slap_scoring'):
        with gzip.open(destination/f'{name}.json.gz','wt') as stream:
            json.dump({str(i):read(run/name/(f'{i}.json' if name=='slap_xyz' else f'{i}.refined.json'))
                       for i in range(140)},stream)
    with gzip.open(destination/'audited_witnesses.json.gz','wt') as stream:
        json.dump(witness_rows,stream)
    (destination/'archive_hashes.json').write_text(json.dumps(hashes,indent=2)+'\n')
    (destination/'audit.json').write_text(json.dumps(dict(
        audited_witnesses=len(witness_rows), counts=dict(Counter(w['method'] for w in witness_rows)),
        full_archives=str(run), hashed_search_and_scoring_artifacts=len(hashes),
        frozen_files_verified=len(manifest['frozen_sha256']),
        slap_source_files_verified=len(manifest['slap_source_sha256']),
        archived_configuration_checks=context_checks,
        input_files_verified=len(read(run/'input_hashes.json')),
        scoring_check='Independent scalar bond_events validates each directional AAM minimum and every scored SLAP candidate.',
        failures=0),indent=2)+'\n')
    timing, result = summary['timing'], summary['event_comparison']
    case123 = rows[123]
    unequal = '\n'.join(f"| {r['index']} | {r['name']} | {r['aam']['best_events']} | {r['slap']['best_events']} |"
                        for r in rows if r['winner']!='tie')
    text = f"""# One seed, branch cap 1000: 140-case comparison with native XYZ SLAP

Fresh experiment on the same 140 XYZ/WBO inputs, using frozen stable AAM engine
`{manifest['original_commit']}`. AAM uses one deterministic seed order per cut,
branch cap 1000, tolerance 1.0, both directions, and the original uncut/single-edge
cut sweep. All other search settings match the stable ten-seed/cap-100 baseline.
SLAP is rerun through its unmodified native XYZ interface with original components,
binary adjacency, `break_sym='heavy'`, and default bond scale 1.2.
This comparator is separate from the Golden-set SLAP single-cut sweep ablation.

| Minimum full-H WBO change score | Cases |
| --- | ---: |
| AAM lower | {result.get('aam',0)} |
| Equal | {result.get('tie',0)} |
| SLAP lower | {result.get('slap',0)} |

Both methods return full mappings for all 140 cases. All 840 execution/scoring
phases succeed, and all SLAP H-label score optimizations finish with certified
bounds. These cases have no reference atom maps: this is an objective-score
comparison, not a chemical-accuracy estimate or proof of a global minimum.
Case 123 scores {case123['aam']['best_events']} for AAM and {case123['slap']['best_events']} for SLAP.

| Recorded time | AAM | Native XYZ SLAP |
| --- | ---: | ---: |
| Total search CPU seconds | {timing['aam']['total_cpu_seconds']:.3f} | {timing['slap']['total_cpu_seconds']:.3f} |
| Mean search CPU seconds/case | {timing['aam']['mean_cpu_seconds']:.3f} | {timing['slap']['mean_cpu_seconds']:.3f} |
| Median search CPU seconds/case | {timing['aam']['median_cpu_seconds']:.3f} | {timing['slap']['median_cpu_seconds']:.3f} |
| Separately measured scoring CPU seconds | {timing['aam']['scoring_cpu_seconds']:.3f} | {timing['slap']['scoring_cpu_seconds']:.3f} |

AAM uses **{summary['aam_cpu_over_slap']:.2f} times** the recorded mapping CPU of
native XYZ SLAP. Each paired case runs on the same host. Case 0 is the retained
eight-CPU local pilot; the other 139 pairs run in eight-CPU Slurm allocations.
AAM uses eight workers and SLAP one, so wall times do not establish equal-core
latency. AAM CPU includes parent and child search, process setup, IPC and profiling,
excluding measured checkpoint persistence/loading. SLAP CPU covers `map_3d`,
including native XYZ processing. Imports and scheduling are excluded.
Scoring is separate: AAM scoring also canonicalizes saved heavy mapping classes;
SLAP scoring resolves native H-label assignments. These are different output
processing workloads, so the main ratio uses the recorded search timers.

The prior cap-100/ten-seed comparison gave six AAM-lower cases, one SLAP-lower
case and 133 ties. Fresh SLAP scores differ from the archived comparison in
{len(summary['differs_from_previous_slap'])} cases.

| Case index (zero-based) | Name | AAM score | SLAP score |
| --- | --- | ---: | ---: |
{unequal}

`per_case.csv` contains all 140 comparisons and timings. `case_metrics.json`
retains directional completion, cap flags, phase status and hardware. The
compressed witnesses were independently rescored with the scalar bond-event
implementation; the published scores use the original vectorized evaluator.
`archive_hashes.json` identifies full local search archives at `{run}`.
The search driver is `bench/holdout_cap_seed_benchmark.py`; this report is built
by `bench/publish_holdout_cap_seed.py`. Existing manuscript baseline panels are
unchanged by this separate experiment.
"""
    (destination/'README.md').write_text(text)
    job = read(run/'submission.json')['job'].split(';')[0]
    command = ['sacct','-j',job,'--array','--noheader','--parsable2',
               '--format=JobID,State,ExitCode,ElapsedRaw,TotalCPU,AllocCPUS,MaxRSS,NodeList']
    (destination/'slurm_accounting.psv').write_text(subprocess.check_output(command,text=True))
    figure(rows,summary,destination)
    print(json.dumps(dict(report=str(destination),audited_witnesses=len(witness_rows),summary=summary),indent=2))


def figure(rows, summary, destination):
    os.environ.setdefault('MPLCONFIGDIR',str(destination/'.mpl-cache'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'pdf.fonttype':42,
                         'axes.spines.top':False,'axes.spines.right':False})
    fig, axes = plt.subplots(1,3,figsize=(11,3.7),layout='constrained')
    counts = Counter((r['slap']['best_events'],r['aam']['best_events']) for r in rows)
    limit = max(max(point) for point in counts)+2
    axes[0].plot([0,limit],[0,limit],ls='--',color='#ABB6BA',lw=1)
    for (x,y),n in counts.items():
        axes[0].scatter(x,y,s=20+15*np.sqrt(n),color='#007F68' if y<=x else '#CC6B26',alpha=.75)
    axes[0].set(xlabel='SLAP: minimum changes',ylabel='AAM: minimum changes',
                xlim=(-.6,limit),ylim=(-.6,limit),title='a  Best saved mapping scores')
    result = summary['event_comparison']
    values = [result.get('aam',0),result.get('tie',0),result.get('slap',0)]
    axes[1].barh([2,1,0],values,color=['#007F68','#ABB6BA','#CC6B26'])
    for y,n in zip([2,1,0],values):
        axes[1].text(n+2,y,str(n),va='center')
    axes[1].set(yticks=[2,1,0],yticklabels=['AAM lower','Tie','SLAP lower'],
                xlim=(0,150),xlabel='Cases',title='b  Minimum-change comparison')
    for row in rows:
        axes[2].scatter(row['slap']['cpu_seconds'],row['aam']['cpu_seconds'],
                        color='#3269A8',alpha=.6,s=17)
    times = [r[m]['cpu_seconds'] for r in rows for m in ('slap','aam')]
    bounds = [min(times)*.7,max(times)*1.4]
    axes[2].plot(bounds,bounds,ls='--',lw=1,color='#ABB6BA')
    axes[2].set(xscale='log',yscale='log',xlim=bounds,ylim=bounds,
                xlabel='SLAP CPU seconds / case',ylabel='AAM CPU seconds / case',
                title='c  Paired mapping CPU')
    for ax in axes:
        ax.grid(alpha=.12)
        ax.set_axisbelow(True)
    fig.get_layout_engine().set(rect=(0,.15,1,.85))
    fig.text(.02,.07,'140 XYZ/WBO cases | AAM: 1 seed order, cap 1000, tolerance 1.0, both directions and cut sweep',fontsize=9)
    fig.text(.02,.02,'Native XYZ SLAP | No reference accuracy labels | CPU excludes offline scoring; paired cases share a host',fontsize=8)
    for ext in ('png','pdf','svg'):
        fig.savefig(destination/f'comparison.{ext}',dpi=200,facecolor='white')
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
