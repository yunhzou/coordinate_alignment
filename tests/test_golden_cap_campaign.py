import json
import signal
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'bench'))
import golden_cap_campaign as campaign


def test_watchdog_terminates_process_group_and_records_incomplete(tmp_path, monkeypatch):
    job = dict(index=1665, seeds=100, cap=2000, tolerance=1.0, workers=8, source='/unused')
    (tmp_path/'manifest.json').write_text(json.dumps(dict(jobs=[job], watchdog_seconds=0)))
    class Child:
        pid = 12345
        calls = 0

        def wait(self, timeout):
            self.calls += 1
            if self.calls == 1:
                raise subprocess.TimeoutExpired('diagnostic', timeout)
            return -signal.SIGTERM

    child = Child()
    kills = []
    monkeypatch.setattr(campaign.subprocess, 'Popen', lambda *a, **k: child)
    monkeypatch.setattr(campaign.os, 'killpg', lambda pid, sig: kills.append((pid, sig)))
    campaign.worker(tmp_path, 0)
    assert kills == [(12345, signal.SIGTERM)]
    status = json.loads((tmp_path/'slot0_status.json').read_text())
    assert status['exit_code'] == 'watchdog'
    assert not (tmp_path/'case1665_seeds100_cap2000/evaluation.json').exists()
