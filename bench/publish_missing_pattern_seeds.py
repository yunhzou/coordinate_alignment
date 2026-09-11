"""Publish the three-case seed-budget diagnostic without changing full benchmarks."""
import argparse
import csv
import gzip
import json
from pathlib import Path
import shutil
import tarfile

from holdout_missing_pattern_seeds import CASES,SEEDS,AAM,EventPatterns,read,save,sha
from rxn_core import AAMProblem
from rxn_core.domain import MolecularEndpoint
from rxn_core.family_scoring import bond_events


def publish(args):
    targets=read(args.run/'targets.json');variants=read(args.run/'variants.json')
    rows=[];witness_checks=0
    base=read(AAM/'manifest.json')['original_config']
    for variant,spec in variants.items():
        if spec['fresh']:
            root=Path(spec['run']);manifest=read(root/'manifest.json')
            assert manifest['original_config']==dict(base,seed_count=spec['seeds'])
            assert read(root/'tasks.json')==[dict(dataset='holdout',index=i,direction='R_to_P') for i in CASES]
        for index in CASES:
            status=read(args.run/f'status/case_{index}.json');assert 'finished' in status
            assert status['phases'][f'{variant}/analyze']['exit']==0
            if spec['fresh']:assert status['phases'][f'{variant}/search']['exit']==0
            initial=read(args.run/f'analysis/{variant}/{index}.json')
            extended=args.run/f'analysis_extended/{variant}/{index}.json'
            row=read(extended) if extended.exists() else initial
            if extended.exists():
                assert all(row[k]==initial[k] for k in ['archive_sha256','search_sha256','search_cpu','target_pattern','minimum_events'])
                row['initial_check']=dict(status=initial['target_result']['status'],query_metrics=initial['query_metrics'],
                    analysis_cpu=initial['analysis_cpu'],analysis_wall=initial['analysis_wall'])
            coarse_path=args.run/f'coarse_checks/{variant}/{index}.json'
            if coarse_path.exists():
                coarse=read(coarse_path)
                assert coarse['archive_sha256']==row['archive_sha256'] and coarse['target_pattern']==row['target_pattern']
                row['coarse_certificate']=coarse
                if row['target_result']['status']=='unresolved' and coarse['complete_exclusion']:
                    assert not coarse['noninvariant_terminals'] and not coarse['matching_invariant_terminals']
                    assert coarse['full_terminals']==row['full_terminals_scanned']
                    row['initial_target_result']=row['target_result']
                    row['target_result']=dict(status='excluded_from_saved_families',method='binary_heavy_event_invariance_certificate',
                                             certificate=str(coarse_path))
            raw=read(AAM/f'inputs/{index}/input.json');canonical=EventPatterns(raw)
            problem=AAMProblem(*(MolecularEndpoint(**raw[k]) for k in ('reactant','product')))
            mapping=row['minimum_witness']['mapping'];actual=canonical.describe(mapping)
            assert actual['total']==row['minimum_events']
            assert bond_events(problem,dict(enumerate(mapping)))['total']==actual['total'];witness_checks+=1
            result=row['target_result']
            assert result['status'] in ('represented','excluded_from_saved_families','unresolved')
            if result['status']=='represented':
                mapping=result['witness']['mapping'];actual=canonical.describe(mapping)
                assert actual['id']==targets[str(index)]['pattern']['id'] and actual['total']==row['target_score']
                assert bond_events(problem,dict(enumerate(mapping)))['total']==actual['total'];witness_checks+=1
            elif result['status']=='excluded_from_saved_families':
                assert row['query_metrics']['exhausted'] or row.get('coarse_certificate',{}).get('complete_exclusion')
            row['host']=status['host'] if spec['fresh'] else read(Path(spec['run'])/f'results/holdout/{index}/R_to_P/original/environment.json').get('host')
            rows.append(row)
    fresh=[r for r in rows if r['fresh']]
    assert len(fresh)==12 and {r['index'] for r in fresh}==set(CASES)
    for index in CASES:
        case=[r for r in fresh if r['index']==index];assert len({r['host'] for r in case})==1
        assert next(r for r in case if r['seeds']==1)['target_result']['status']=='excluded_from_saved_families'
        assert next(r for r in case if r['seeds']==1)['minimum_events']==targets[str(index)]['baseline_minimum']
        ordered=sorted(case,key=lambda r:r['seeds'])
        assert all(b['minimum_events']<=a['minimum_events'] for a,b in zip(ordered,ordered[1:]))
    summary=dict(cases=list(CASES),direction='R_to_P',cap=1000,seeds=list(SEEDS),
        recovery={str(seed):dict(recovered=[r['index'] for r in fresh if r['seeds']==seed and r['target_result']['status']=='represented'],
            excluded=[r['index'] for r in fresh if r['seeds']==seed and r['target_result']['status']=='excluded_from_saved_families'],
            unresolved=[r['index'] for r in fresh if r['seeds']==seed and r['target_result']['status']=='unresolved'],
            search_cpu=sum(r['search_cpu'] for r in fresh if r['seeds']==seed)) for seed in SEEDS},
        fresh_mapping_runs=12,full_benchmark_rerun=False,independently_checked_witnesses=witness_checks,
        capped_fresh_runs=[dict(index=r['index'],seeds=r['seeds']) for r in fresh if r['capped']],
        fresh_search_cpu_total=sum(r['search_cpu'] for r in fresh),
        analysis_cpu_total=sum(r['analysis_cpu']+r.get('initial_check',{}).get('analysis_cpu',0)+
                               r.get('coarse_certificate',{}).get('cpu',0) for r in rows),
        analysis_cost_scope='Completed timed analysis stages and certificates. Cancelled general-query retry is recorded separately in Slurm accounting.',
        cancelled_general_query_cpu=None,
        cancelled_query_accounting_note='No validated inclusive descendant CPU total survived cancellation. Slurm counters are retained verbatim; completed-analysis CPU is not the whole campaign cost.',
        scope='Recovery of three specific missing bond-edit patterns at their previously recorded scores. '
              'Targeted diagnostic, not 140-case accuracy or total alternative coverage. '
              'Minimum seed count between tested budgets is not established.')
    dest=args.destination;dest.mkdir(parents=True,exist_ok=False)
    for name,value in [('summary.json',summary),('per_case.json',rows),('targets.json',targets)]:save(dest/name,value)
    with (dest/'per_case.csv').open('w',newline='') as stream:
        writer=csv.writer(stream,lineterminator='\n')
        writer.writerow(['case','seeds','cap','fresh','minimum_events','target_score','target_status','in_terminal_representatives',
                         'mapping_cpu','mapping_wall','analysis_cpu','capped','host'])
        for r in rows:writer.writerow([r['index'],r['seeds'],r['cap'],r['fresh'],r['minimum_events'],r['target_score'],
            r['target_result']['status'],r['target_in_terminal_representatives'],r['search_cpu'],r['search_wall'],r['analysis_cpu'],r['capped'],r['host']])
    for name in ['manifest.json','submission.json','variants.json','slurm_accounting.tsv']:shutil.copy2(args.run/name,dest/name)
    if (args.run/'extended_submission.json').exists():shutil.copy2(args.run/'extended_submission.json',dest/'extended_submission.json')
    for name in ['extended_cancellation.json','coarse_validation.json']:
        if (args.run/name).exists():shutil.copy2(args.run/name,dest/name)
    for p in args.run.glob('*relocation*.json'):shutil.copy2(p,dest/p.name)
    records={}
    for variant,spec in variants.items():
        for index in CASES:
            folder=Path(spec['run'])/f'results/holdout/{index}/R_to_P/original'
            records[f'{variant}/{index}']={name:read(folder/name) for name in ['search.json','environment.json']}
    with gzip.open(dest/'search_records.json.gz','wt') as stream:json.dump(records,stream)
    with gzip.open(dest/'worker_status.json.gz','wt') as stream:json.dump({p.stem:read(p) for p in (args.run/'status').glob('*.json')},stream)
    with tarfile.open(dest/'analysis_sources.tar.gz','w:gz') as archive:
        for p in args.run.glob('*.py'):archive.add(p,arcname=f'frozen/{p.name}')
        archive.add(__file__,arcname=f'published/{Path(__file__).name}')
    save(dest/'analysis_hashes.json',dict(published_driver=sha(Path(__file__)),
        source_results={str(p):sha(p) for folder in ['analysis','analysis_extended','coarse_checks'] for p in args.run.glob(f'{folder}/*/*.json')}))
    labels={'represented':'Recovered','excluded_from_saved_families':'Absent from saved families','unresolved':'Unresolved'}
    table='\n'.join(f"| {i} | {targets[str(i)]['name']} | {targets[str(i)]['pattern']['total']} | "+' | '.join(
        labels[next(r for r in fresh if r['index']==i and r['seeds']==seed)['target_result']['status']] for seed in SEEDS)+' |' for i in CASES)
    timing='\n'.join(f"| {i} | "+' | '.join(f"{next(r for r in fresh if r['index']==i and r['seeds']==seed)['search_cpu']:.3f}" for seed in SEEDS)+' |' for i in CASES)
    cached='\n'.join(f"- Case {r['index']}: {labels[r['target_result']['status']]}; capped={r['capped']}." for r in rows if not r['fresh'])
    (dest/'README.md').write_text(f"""# Can additional seeds recover the three missing AAM patterns?

This follow-up tests **only cases 11, 64 and 101**, the three known alternative
omissions in the forward one-seed/cap-1000 benchmark. It does not rerun the 140-case
benchmark. Each selected case uses 1, 3, 10 and 30 seed orders per cut, with cap
1000, tolerance 1.0, explicit H, the original uncut/single-edge sweep and
**reactant → product only**. The stable frozen engine remains unchanged.

## Recovery of the specific missing patterns

| Case | Name | Target event score | 1 seed | 3 seeds | 10 seeds | 30 seeds |
| --- | --- | ---: | --- | --- | --- | --- |
{table}

“Recovered” requires a validated mapping with the same joint symmetry-normalized
broken/formed/strengthened/weakened bond pattern. “Absent from saved families”
requires scanning all full terminal representatives and completing the correlated
family membership checks, or a complete necessary-condition exclusion certificate.
It concerns the saved output for that budget, not every
possible seed or every possible mapping. An unresolved query never proves absence.

The target event score is held fixed even if a new search finds a lower score.
`per_case.csv` records each run's best saved representative score separately.
The first successful tested budget does not identify the exact minimum seed count
between tested budgets. This is a selected-case recovery diagnostic, not a new
whole-dataset coverage estimate or a validation of physical reaction pathways.

## Mapping CPU seconds

| Case | 1 seed | 3 seeds | 10 seeds | 30 seeds |
| --- | ---: | ---: | ---: | ---: |
{timing}

All four fresh budgets for each case run sequentially on the same host with eight
workers. CPU includes the parent and all workers; **do not multiply it by eight**.
The frozen profiler excludes measured checkpoint persistence/loading. Mapping
CPU and observed wall time, which includes archive I/O, are recorded separately.
Analysis is also separate: exact target membership can cost more than mapping.
These are single measurements on three selected cases, not a full timing benchmark.

The {summary['fresh_mapping_runs']} fresh mapping runs consumed
{summary['fresh_search_cpu_total']:.3f} CPU s in total. Completed timed analysis
stages, including the three cached checks and the certificates, consumed
{summary['analysis_cpu_total']:.3f} CPU s. The interrupted general-query retry is
separately recorded in Slurm accounting and is excluded from this completed-stage
subtotal. Its inclusive CPU total is unavailable after cancellation; the recorded
Slurm counters are retained verbatim. The completed-stage subtotal is not the
whole campaign cost.
Searches with cap flags: {summary['capped_fresh_runs']}.

The initial 30-seed case-11 family check reached its 240-second analysis limit.
Its unresolved result is preserved; `extended_submission.json` records the longer
check on the same archive, with no mapping rerun. That redundant query was stopped
after an independent certificate excluded the target from all 5,637 full-terminal
families. Per-case records retain the initial unresolved status and its analysis cost.

The certificate checks every recorded family action against the product's binary
heavy-atom graph. These actions preserve each terminal's broken/formed heavy-bond
set, and none of the 5,637 sets matches the target under the benchmark's reactant
symmetry. This is sufficient to exclude the full target regardless of hydrogen
or bond-order choices. The certificate was checked on all 12 fresh outputs: it
agrees with the completed queries and does not exclude either known positive
case-101 output. `coarse_validation.json` records 366 literal invariance checks.

## Existing ten-seed, cap-100 outputs

The earlier saved ten-seed run was inspected without rerunning its search:

{cached}

That historical run changes both the seed count and the cap relative to the
one-seed/cap-1000 benchmark. Its results are contextual evidence; the controlled
seed comparison uses the fresh cap-1000 rows above.

## Verification and reproduction

The three targets are frozen before the diagnostic. The mapping subprocess uses
only endpoint inputs, seed configuration and the original cut policy. Target
patterns are read solely by the post-search analysis. Input hashes, configuration
checks, manifests and commands establish that only these three forward cases ran.
The one-seed control reproduces the three original omissions. Every published
minimum witness and every positive target witness was independently rescored
({witness_checks} scalar/canonical checks).

`targets.json` and `per_case.json` include the bond edits, actual mappings and
family-query provenance. Compressed search records, source hashes, frozen analysis
scripts and Slurm accounting accompany this report. Large search archives remain
at `{args.run}` and retain their hashes. Existing full benchmark reports are
unchanged. Run `bench/holdout_missing_pattern_seeds.py prepare` with a fresh
`--run`, then `submit`; its task list is restricted to these three cases. After
the case analyses finish, run `bench/certify_heavy_event_exclusion.py --run` on
that directory. This supplies the independent exclusion certificates used here.
Finally publish with `bench/publish_missing_pattern_seeds.py --run` and a fresh
`--destination`. The submission record supplies the frozen Python/source paths.
""")
    print(json.dumps(summary,indent=2),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run',type=Path,required=True)
    parser.add_argument('--destination',type=Path,required=True)
    publish(parser.parse_args())
