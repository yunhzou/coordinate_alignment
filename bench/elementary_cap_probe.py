"""Isolated cap ablation using the holdout's unchanged frozen engine and inputs."""
import argparse
from collections import Counter
import json
from pathlib import Path
import shlex
import shutil
import subprocess
import sys

from golden_competitors import save
from elementary_feasibility import digest


def prepare(args):
    args.run.mkdir(parents=True,exist_ok=False)
    (args.run/'inputs').mkdir();(args.run/'status').mkdir()
    shutil.copytree(args.source/'inputs'/str(args.index),args.run/'inputs'/str(args.index))
    manifest=json.loads((args.source/'inputs/manifest.json').read_text())
    before=dict(manifest['config'])
    manifest['config']['branch_limit']=args.cap
    manifest['records']=[r for r in manifest['records'] if r['index']==args.index]
    assert len(manifest['records'])==1
    save(args.run/'inputs/manifest.json',manifest)
    (args.run/'engine').symlink_to((args.source/'engine').resolve(),target_is_directory=True)
    save(args.run/'manifest.json',dict(source=str(args.source),index=args.index,
        config_before=before,config=manifest['config'],records=manifest['records'],
        frozen_engine=str((args.source/'engine').resolve()),
        source_manifest_sha256=digest(args.source/'manifest.json'),
        input_sha256=digest(args.run/'inputs'/str(args.index)/'input.json'),
        only_search_change='branch_limit',reference_available=False))
    engine=args.run/'engine'
    command=['env',f'PYTHONPATH={engine}/src:{engine}/bench','PYTHONHASHSEED=0',
        'RXN_CORE_NATIVE=1','OMP_NUM_THREADS=1','OPENBLAS_NUM_THREADS=1','MKL_NUM_THREADS=1',
        sys.executable,str(engine/'bench/elementary_feasibility.py'),'worker','--run',str(args.run),
        '--method','aam','--slot']
    options=['sbatch','--parsable','--partition=cpunodes','--exclude=bosque8',
        '--nodes=1','--cpus-per-task=16','--mem=32G','--time=00:10:00',
        f'--array={2*args.index},{2*args.index+1}',f'--job-name=elem_cap{args.cap}',
        f'--output={args.run}/status/slurm_%A_%a.out','--wrap',
        shlex.join(command)+' "$SLURM_ARRAY_TASK_ID"']
    job=subprocess.check_output(options,text=True).strip()
    save(args.run/'submission.json',dict(job=job,command=options))
    print(job)


def report(args):
    from rxn_core.artifacts import read_aam_checkpoint
    manifest=json.loads((args.run/'manifest.json').read_text());index=manifest['index']
    rows=[]
    for root,cap in [(Path(manifest['source']),manifest['config_before']['branch_limit']),
                     (args.run,manifest['config']['branch_limit'])]:
        for direction in ['R_to_P','P_to_R']:
            out=root/'directions'/str(index)/direction
            status=json.loads((root/'status'/f'aam_{index}_{direction}.json').read_text())
            row=dict(cap=cap,direction=direction,worker_exit=status['exit'])
            if (out/'feasibility.json').exists():
                row.update(json.loads((out/'feasibility.json').read_text()))
                timing=json.loads((out/'search.json').read_text())
                row.update(cpu_seconds=timing['compute_cpu_excluding_persistence_and_loading_seconds'],
                    elapsed_including_io_seconds=timing['elapsed_wall_including_io_seconds'])
                graph=read_aam_checkpoint(out/'cuts/aam.pkl.gz').graph
                row.update(contexts=len(graph.contexts),stops=dict(Counter(s.reason for s in graph.stops)),
                    cap_stages=dict(Counter(s.stage for s in graph.stops if s.reason=='capped')))
            rows.append(row)
    save(args.run/'comparison.json',dict(index=index,results=rows,accuracy=None))
    print(json.dumps(rows,indent=2))


def publish(args):
    report(args)
    out=Path('reports/elementary140_feasibility_20260908/case123_cap200')
    out.mkdir(exist_ok=True)
    for name in ('comparison.json','manifest.json','submission.json'):
        shutil.copy2(args.run/name,out/name)
    shutil.copy2(__file__,args.run/'probe_driver.py')
    for direction in ('R_to_P','P_to_R'):
        dest=out/direction;dest.mkdir(exist_ok=True)
        source=args.run/'directions/123'/direction
        for name in ('search.json','feasibility.json','environment.json','terminal_mappings.jsonl.gz'):
            shutil.copy2(source/name,dest/name)
    job=json.loads((args.run/'submission.json').read_text())['job']
    accounting=subprocess.check_output(['sacct','-j',job,'-P',
        '--format=JobID,State,ExitCode,ElapsedRaw,TotalCPU,AllocCPUS,MaxRSS,NodeList'],text=True)
    (out/'slurm_accounting.psv').write_text(accounting)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command',choices=('prepare','report','publish'))
    parser.add_argument('--source',type=Path)
    parser.add_argument('--run',type=Path,required=True)
    parser.add_argument('--index',type=int,default=123)
    parser.add_argument('--cap',type=int,default=200)
    args=parser.parse_args();globals()[args.command](args)
