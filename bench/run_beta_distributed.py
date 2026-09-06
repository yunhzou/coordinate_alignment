#!/usr/bin/env python3
"""Run the beta policy using disk-indexed, Slurm-sharded connected queries."""
import argparse
from collections.abc import Mapping
from dataclasses import dataclass
import csv
import gzip
import hashlib
from heapq import merge
import json
from pathlib import Path
import pickle
import sqlite3
import subprocess
import time
import os
import resource

from rxn_core.retrosynthesis.slurm_budget import SlurmBudget

from rxn_core.fragment_matching import FragmentDetectionConfig
from rxn_core.retrosynthesis.beta import FragmentQueryBank, recommend_big_blocks, proposal_rank
from rxn_core.smiles import smiles_to_weighted_graph


class CatalogGraphs(Mapping):
    """Materialize only selected source graphs, never all bank molecules."""
    def __init__(self, smiles):
        self.smiles = smiles
    def __iter__(self):
        return iter(self.smiles)
    def __len__(self):
        return len(self.smiles)
    def __getitem__(self, key):
        return smiles_to_weighted_graph(self.smiles[key], expand_hydrogens=True)


@dataclass(frozen=True)
class IndexedPlacement:
    """A proposal descriptor; its full AAM evidence stays in the source archive."""
    source_id: str
    mapping: tuple
    partition: tuple
    input_atom_count: int
    block_file: str
    block_index: int
    refined: bool = False
    @property
    def covered_atoms(self):
        return frozenset(p for _, p in self.mapping)
    @property
    def fragment_count(self):
        return len(self.partition)
    @property
    def key(self):
        return self.source_id, self.mapping, self.partition, self.refined


class DistributedBank(FragmentQueryBank):
    def __init__(self, root, excluded, budget=None, poll_seconds=2,
                 partitions='cpunodes_nia,cpunodes'):
        self.root = root
        self.manifest = json.loads((root/'query_full/manifest.json').read_text())
        target = smiles_to_weighted_graph(self.manifest['target_smiles'], expand_hydrogens=True)
        super().__init__([], target, config=FragmentDetectionConfig(**self.manifest['config']),
                         checkpoint=self.save_refinement)
        with gzip.open(self.manifest['catalog'], 'rt') as stream:
            rows = list(csv.DictReader(stream))
        self.sources = CatalogGraphs({r[self.manifest['id_column']]:r['SMILES'] for r in rows})
        self.excluded = excluded
        self.budget = budget
        self.jobs = []
        self.poll_seconds = poll_seconds
        self.partitions = partitions
        self.queried = set()
        self.completed_assemblies = []

    def save_refinement(self, event, evidence):
        directory = self.root/'refinements'
        directory.mkdir(exist_ok=True)
        path = directory/f"{event['source_id']}.pkl.gz"
        if path.exists():
            print(json.dumps(dict(event,checkpoint=str(path),reused=True)),flush=True)
            return
        with gzip.open(path,'wb',compresslevel=1) as stream:
            pickle.dump((event,evidence),stream,protocol=pickle.HIGHEST_PROTOCOL)
        print(json.dumps(dict(event,checkpoint=str(path))),flush=True)

    def detect_selected(self, source_id):
        path=self.root/'refinements'/f'{source_id}.pkl.gz'
        if path.exists():
            with gzip.open(path,'rb') as stream:
                _,evidence=pickle.load(stream)
            return evidence
        return super().detect_selected(source_id)

    def assembly_found(self, assembly, pattern, pattern_count, count):
        self.completed_assemblies.append(assembly)
        directory=self.root/'assembly_candidates'
        directory.mkdir(exist_ok=True)
        key=hashlib.sha256(repr(tuple(p.key for p in assembly.placements)).encode()).hexdigest()
        path=directory/f'{key}.pkl.gz'
        if not path.exists():
            with gzip.open(path,'wb',compresslevel=1) as stream:
                pickle.dump(assembly,stream,protocol=pickle.HIGHEST_PROTOCOL)
        print(json.dumps(dict(stage='assembled_cover',count=count,patterns=pattern_count,
            sources=[p.source_id for p in assembly.placements],
            fragments=sum(p.fragment_count for p in assembly.placements),
            checkpoint=str(path))),flush=True)

    def submit(self, script, directory, shards, label):
        if self.budget is None:
            raise ValueError('New Slurm work requires an explicit CPU budget')
        command=['sbatch','--parsable',*self.budget.arguments(shards),
            f'--partition={self.partitions}',
            f'--output={directory}/logs/{label}_%A_%a.out',
            f'--error={directory}/logs/{label}_%A_%a.err',script,str(directory)]
        if self.excluded:
            command.insert(2, f'--exclude={self.excluded}')
        environment=dict(os.environ)
        # Nested submission must not inherit the coordinator's memory mode.
        for name in ('SLURM_MEM_PER_CPU','SLURM_MEM_PER_GPU','SLURM_MEM_PER_NODE'):
            environment.pop(name,None)
        job=subprocess.check_output(command,text=True,env=environment).strip().split(';')[0]
        self.jobs.append(job)
        print(json.dumps(dict(stage='submitted',job=job,query=str(directory),kind=label)),flush=True)
        return job

    def wait_parts(self, directory, subdirectory, job):
        expected=self.manifest['shards']
        while len(list((directory/subdirectory).glob('*.json'))) != expected:
            # Pending jobs consume no CPU budget. The worker's own 580-second
            # timeout governs computation; queue time is not a failed search.
            states=subprocess.check_output(['squeue','-h','-j',job,'-o','%T'],text=True).splitlines()
            if not states:
                accounting=subprocess.check_output(['sacct','-X','-n','-P','-j',job,
                    '--format=State,ExitCode'],text=True).splitlines()
                if accounting and all(line.strip().split('|')[:2] in
                        (['COMPLETED','0:0'],['FAILED','75:0']) for line in accounting):
                    # Code 75 explicitly yields after saving completed sources.
                    # Restart only unfinished shards; do not rerun saved AAM.
                    return
                raise RuntimeError(f'Query jobs failed: {job}; checkpoints retained; {accounting}')
            print(json.dumps(dict(stage='waiting',query=str(directory),kind=subdirectory,
                finished=len(list((directory/subdirectory).glob('*.json'))),expected=expected)),flush=True)
            time.sleep(self.poll_seconds)

    def ensure_query(self, region):
        region=tuple(sorted(region))
        name='query_full' if len(region)==len(self.target) else 'query_'+hashlib.sha256(
            json.dumps(region).encode()).hexdigest()[:16]
        directory=self.root/name
        if not directory.exists():
            directory.mkdir()
            for part in ('sources','parts','logs','indexes'):
                (directory/part).mkdir()
            manifest=dict(self.manifest,region=region)
            (directory/'manifest.json').write_text(json.dumps(manifest)+'\n')
        missing=[str(i) for i in range(self.manifest['shards'])
                 if not (directory/'parts'/f'{i}.json').exists()
                 or not (directory/'indexes'/f'{i}.json').exists()]
        while missing:
            job=self.submit('hpc/query_beta_catalog.sbatch',directory,','.join(missing),'query')
            self.wait_parts(directory,'indexes',job)
            # Do not submit the next wave until allocations from this wave
            # have exited: array throttles are per array, not per run.
            while subprocess.check_output(['squeue','-h','-j',job,'-o','%T'],text=True).strip():
                time.sleep(self.poll_seconds)
            missing=[str(i) for i in range(self.manifest['shards'])
                     if not (directory/'indexes'/f'{i}.json').exists()]
        if name not in self.queried:
            parts=[json.loads(p.read_text()) for p in (directory/'parts').glob('*.json')]
            self.capped_searches+=sum(p['capped'] for p in parts)
            event=dict(stage='indexed_query',query=name,region=region,
                       rows=sum(p['rows'] for p in parts),blocks=sum(p['blocks'] for p in parts))
            self.events.append(event)
            print(json.dumps(event),flush=True)
            self.queried.add(name)
        return directory

    def ordered_query(self, region, placements):
        directory=self.ensure_query(region)
        selected=sorted({p.source_id for p in placements})
        def stream(path):
            with sqlite3.connect(f'file:{path}?mode=ro',uri=True) as db:
                preference=('CASE WHEN source IN ('+','.join('?' for _ in selected)+') THEN 0 ELSE 1 END,') if selected else ''
                sql='SELECT atoms,source,row,ordinal,mapping,partition FROM blocks ORDER BY coverage DESC,fragments,'+preference+'atoms,source,row,ordinal'
                for atoms,source,row,ordinal,mapping,partition in db.execute(sql,selected):
                    yield IndexedPlacement(source,tuple(map(tuple,json.loads(mapping))),
                        tuple(map(tuple,json.loads(partition))),atoms,
                        str(directory/'sources'/f'{row}.blocks.pkl.gz'),ordinal)
        streams=[stream(directory/'indexes'/f'{i}.sqlite') for i in range(self.manifest['shards'])]
        return merge(*streams,key=lambda block:proposal_rank(placements+(block,),len(self.target)))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run',type=Path,required=True)
    p.add_argument('--exclude',default='')
    p.add_argument('--cpu-budget',type=int,required=True)
    p.add_argument('--worker-cpus',type=int,default=1)
    p.add_argument('--catalog',type=Path)
    p.add_argument('--target-smiles')
    p.add_argument('--shards',type=int)
    p.add_argument('--poll-seconds',type=float,default=2)
    p.add_argument('--partitions',default='cpunodes_nia,cpunodes')
    p.add_argument('--recommendations',type=int,default=20)
    p.add_argument('--patterns',type=int,help='Construction patterns (default: min(4, recommendations))')
    p.add_argument('--result-stem',default='assemblies')
    a=p.parse_args()
    started=time.perf_counter()
    reserved = int(os.environ.get('SLURM_CPUS_ON_NODE', '1'))
    # CR_Core clusters can reserve all SMT threads for a one-CPU request.
    # Round workers to a whole core of the largest eligible topology, and use
    # those CPUs rather than accounting only for the requested one thread.
    topology=subprocess.check_output(['sinfo','-h','-N','-p',a.partitions,'-o','%Z'],text=True)
    quantum=max(int(value) for value in topology.split())
    budget = SlurmBudget(a.cpu_budget, a.worker_cpus, reserved, quantum)
    if not a.run.exists():
        if a.catalog is None or a.target_smiles is None or a.shards is None or a.shards < 1:
            p.error('A new run requires --catalog, --target-smiles and positive --shards')
        with gzip.open(a.catalog,'rt') as stream:
            bank_rows=sum(1 for _ in csv.DictReader(stream))
        query=a.run/'query_full'
        for directory in ('sources','parts','logs','indexes'):
            (query/directory).mkdir(parents=True,exist_ok=True)
        target=smiles_to_weighted_graph(a.target_smiles,expand_hydrogens=True)
        manifest=dict(workflow='big_blocks_beta/v1',catalog=str(a.catalog.resolve()),
            id_column='Bank ID',target_smiles=a.target_smiles,region=list(range(len(target.nodes))),
            shards=a.shards,bank_rows=bank_rows,
            config=dict(iso_tolerance=1.0,branch_limit=100,seed_mode='all'))
        (query/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    if a.poll_seconds <= 0:
        p.error('--poll-seconds must be positive')
    bank=DistributedBank(a.run,a.exclude,budget,a.poll_seconds,a.partitions)
    try:
        result=recommend_big_blocks(bank,recommendations=a.recommendations,pattern_limit=a.patterns)
    except BaseException:
        for job in bank.jobs:
            subprocess.run(['scancel',job],check=False)
        raise
    with gzip.open(a.run/f'{a.result_stem}.pkl.gz','wb',compresslevel=1) as stream:
        pickle.dump(result,stream,protocol=pickle.HIGHEST_PROTOCOL)
    def record(recommendation):
        return dict(uncovered_target_atoms=recommendation.uncovered_target_atoms,
            placements=[dict(source_id=p.source_id,mapping=p.mapping,
                retained_fragments=p.candidate.retained_fragments,
                target_fragments=[sorted(dict(p.mapping)[atom] for atom in fragment)
                                 for fragment in p.candidate.retained_fragments])
                for p in recommendation.placements])
    from rxn_core.retrosynthesis.beta_assembly import placement_pattern, pareto_assembly_ranks, assembly_key
    pareto_ranks=pareto_assembly_ranks(bank.completed_assemblies,bank.target)
    pattern_ids={}
    def patterned(recommendation):
        certificate=placement_pattern(recommendation.placements,bank.target)
        pattern_ids.setdefault(certificate,len(pattern_ids)+1)
        return dict(record(recommendation),construction_pattern=pattern_ids[certificate],
                    pareto_layer=pareto_ranks[assembly_key(recommendation)][0])
    report=dict(workflow='big_blocks_beta/v2',exhaustive=False,
        coordinator_seconds=time.perf_counter()-started,
        recommendation_count=len(result.recommendations),capped_searches=result.capped_searches,
        recommendations=[patterned(r) for r in result.recommendations],
        construction_patterns=[dict(pattern=i,certificate=k) for k,i in pattern_ids.items()],
        ranking_scope='Ranked among discovered beta assemblies; not globally certified over the bank',
        ranking_objective=['maximize explicit-H retention', 'minimize cuts + connections'],
        ranking_method='Pareto layers; within-layer order is display only; fewer fragments then species break equal-objective ties',
        nondominated_in_pool=sum(rank[0]==1 for rank in pareto_ranks.values()),
        nondominated_displayed=sum(pareto_ranks[assembly_key(r)][0]==1 for r in result.recommendations),
        best_partial=record(result.best_partial))
    report['resource_budget']=dict(cpu_limit=budget.cpu_limit,
        coordinator_cpus=budget.coordinator_cpus,worker_cpus=budget.worker_cpus,
        maximum_concurrent_workers=budget.concurrency, jobs=bank.jobs)
    report['resource_budget']['allocation_quantum']=quantum
    workers=[json.loads(path.read_text()) for path in a.run.glob('query_*/resources/*.json')]
    own=resource.getrusage(resource.RUSAGE_SELF)
    children=resource.getrusage(resource.RUSAGE_CHILDREN)
    report['resource_usage']=dict(
        coordinator_actual_cpu_seconds=own.ru_utime+own.ru_stime+children.ru_utime+children.ru_stime,
        worker_actual_cpu_seconds=sum(r['actual_cpu_seconds'] for r in workers),
        worker_process_allocated_cpu_seconds=sum(r['allocated_cpu_seconds'] for r in workers),
        worker_reports=len(workers),
        scope='Completed worker process lifetimes across attempts; excludes Slurm startup/cleanup. '
              'Hard-killed workers may lack reports; use Slurm accounting for those allocations.')
    (a.run/f'{a.result_stem}.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report),flush=True)


if __name__=='__main__':
    main()
