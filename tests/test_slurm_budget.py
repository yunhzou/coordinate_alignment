"""Resource admission and checkpointed execution stay independent of chemistry."""
import json
from pathlib import Path
import subprocess
import sys

import pytest

from rxn_core.retrosynthesis.slurm_budget import SlurmBudget


@pytest.mark.parametrize('limit', [2, 3, 8, 32, 64, 100])
def test_array_admission_includes_coordinator_and_fits_small_budgets(limit):
    for workers in range(1, limit):
        budget = SlurmBudget(limit, workers)
        assert budget.concurrency * workers + 1 <= limit
        assert (budget.concurrency + 1) * workers + 1 > limit
        assert f'--array=0-999%{budget.concurrency}' in budget.arguments('0-999')


def test_budget_cannot_silently_overallocate():
    with pytest.raises(ValueError):
        SlurmBudget(4, 4)


def test_submission_applies_budget_and_clears_parent_memory_mode(tmp_path, monkeypatch):
    import importlib.util
    spec=importlib.util.spec_from_file_location('budget_runner',
        Path(__file__).parents[1]/'bench/run_beta_distributed.py')
    runner=importlib.util.module_from_spec(spec)
    sys.modules[spec.name]=runner
    spec.loader.exec_module(runner)
    bank=object.__new__(runner.DistributedBank)
    bank.budget=SlurmBudget(8, 2)
    bank.excluded=''
    bank.jobs=[]
    bank.partitions='cpunodes'
    monkeypatch.setenv('SLURM_MEM_PER_NODE','4096')
    calls=[]
    def submit(command, *, text, env):
        calls.append(command)
        assert 'SLURM_MEM_PER_NODE' not in env
        return '12345\n'
    monkeypatch.setattr(runner.subprocess,'check_output',submit)
    bank.submit('worker.sbatch',tmp_path,'0-99','query')
    assert '--array=0-99%3' in calls[0]
    assert '--cpus-per-task=2' in calls[0]
    assert bank.jobs == ['12345']


def test_smt_allocations_are_budgeted_and_used():
    budget=SlurmBudget(3,worker_cpus=1,coordinator_cpus=1,allocation_quantum=2)
    assert budget.worker_cpus == 2
    assert budget.concurrency == 1
    assert '--cpus-per-task=2' in budget.arguments('0-63')
    with pytest.raises(ValueError):
        SlurmBudget(2,worker_cpus=1,coordinator_cpus=1,allocation_quantum=2)


def test_fused_worker_yields_and_reuses_evidence_on_resume(tmp_path):
    import gzip
    query = tmp_path/'query'
    for part in ('sources','indexes','parts'):
        (query/part).mkdir(parents=True)
    catalog = tmp_path/'bank.csv.gz'
    with gzip.open(catalog,'wt') as stream:
        stream.write('Bank ID,SMILES\nmethane,C\nwater,O\n')
    (query/'manifest.json').write_text(json.dumps(dict(
        target_smiles='CO',catalog=str(catalog),id_column='Bank ID',
        region=list(range(6)),shards=1,bank_rows=2,
        config={'iso_tolerance':1.0,'branch_limit':100})))
    command=[sys.executable,str(Path(__file__).parents[1]/'bench/execute_beta_query.py'),
        '--query',str(query),'--shard','0','--workers','1']
    paused=subprocess.run(command+['--wall-seconds','0'],capture_output=True,text=True)
    assert paused.returncode == 75, paused.stderr
    evidence=query/'sources/0.connected.pkl.gz'
    saved_time=evidence.stat().st_mtime_ns
    assert not (query/'parts/0.json').exists()
    completed=subprocess.run(command,capture_output=True,text=True)
    assert completed.returncode == 0, completed.stderr
    assert evidence.stat().st_mtime_ns == saved_time
    assert (query/'indexes/0.sqlite').exists()
    report=json.loads((query/'resources/0.local.json').read_text())
    assert [s['stage'] for s in report['stages']] == ['scan','index']
    assert report['actual_cpu_seconds'] > 0
    assert report['allocated_cpus'] == 1
