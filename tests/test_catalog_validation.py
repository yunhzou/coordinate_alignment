import gzip
import json
from pathlib import Path
import subprocess
import sys


def test_saved_bank_validation_checks_identity_and_reports_missing(tmp_path):
    bank = tmp_path / 'bank.csv.gz'
    with gzip.open(bank, 'wt') as stream:
        stream.write('SMILES,Inventory ID\nC,a\nO,b\n')
    parts = tmp_path / 'parts'
    parts.mkdir()
    record = dict(source_id='a', row_index=0, representation='C', status='matched',
                  complete=True, candidates=[{}], best_fragment_size=1,
                  capped_seed_count=0, branch_limit=100,
                  timing={'detection_seconds': 0.1})
    def save(records):
        with gzip.open(parts / 'part_0.jsonl.gz', 'wt') as stream:
            for item in records:
                stream.write(json.dumps(item) + '\n')
    report = tmp_path / 'validated.json'
    command = [sys.executable, str(Path(__file__).parents[1] / 'bench/validate_catalog_run.py'),
               '--catalog', str(bank), '--parts', str(parts), '--output', str(report),
               '--workers', '1']
    save([record])
    incomplete = subprocess.run(command, capture_output=True, text=True)
    assert incomplete.returncode == 1
    assert json.loads(report.read_text())['missing_rows'] == [1]
    second = {**record, 'source_id': 'b', 'row_index': 1, 'representation': 'O'}
    save([record, second])
    complete = subprocess.run(command, capture_output=True, text=True)
    assert complete.returncode == 0, complete.stderr
    assert json.loads(report.read_text())['complete']
    save([record, {**second, 'representation': 'N'}])
    mismatch = subprocess.run(command, capture_output=True, text=True)
    assert mismatch.returncode == 1
    assert 'source identity mismatch' in mismatch.stderr
    save([record, record])
    duplicate = subprocess.run(command, capture_output=True, text=True)
    assert duplicate.returncode == 1
    assert 'duplicate bank row' in duplicate.stderr
