import json
from dataclasses import asdict
from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'bench'))
from decode_saved_events import same_config, verified_resume, certificate_rows
from rxn_core import AAMSearchConfig


def test_json_roundtrip_is_the_same_search_configuration():
    config = AAMSearchConfig(seed_count=1, branch_limit=2000)
    recorded = json.loads(json.dumps(asdict(config)))
    assert same_config(recorded, config)
    assert not same_config(dict(recorded, seed_count=3), config)


def test_resume_rejects_different_archive_policy_or_window(tmp_path):
    policy = dict(threshold=.5, metal_threshold=.3)
    record = dict(archive_sha256='archive', event_policy=policy, max_events=5)
    (tmp_path / 'comparison.json').write_text(json.dumps(record))
    assert verified_resume(tmp_path, 'archive', policy, 5) == record
    for digest, other_policy, window in [('other',policy,5), ('archive',dict(policy,threshold=.2),5), ('archive',policy,None)]:
        with pytest.raises(ValueError):
            verified_resume(tmp_path, digest, other_policy, window)


def test_only_an_incomplete_final_certificate_line_is_ignored(tmp_path):
    path = tmp_path / 'families.jsonl'
    row = dict(complete=True, terminal=1, transitions=[0])
    path.write_text(json.dumps(row) + '\n{"complete":tr')
    assert list(certificate_rows(path)) == [row]
    path.write_text(json.dumps(row) + '\n{broken}\n')
    with pytest.raises(json.JSONDecodeError):
        list(certificate_rows(path))
