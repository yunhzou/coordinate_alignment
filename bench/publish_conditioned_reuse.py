"""Package saved reuse experiments without rerunning any matching."""
import argparse
import hashlib
import io
import json
from pathlib import Path
import shutil
import subprocess
import tarfile
import xml.etree.ElementTree as ET


def main(args):
    args.output.mkdir(parents=True, exist_ok=True)
    sources = {name: getattr(args, name) for name in ('pilot', 'refined', 'full', 'public')}
    for name in ('pilot', 'refined', 'full'):
        analysis = json.loads((sources[name] / 'analysis.json').read_text())
        assert not analysis['missing'] and not analysis['mismatches'], name
        shutil.copy2(sources[name] / 'analysis.json', args.output / f'{name}_analysis.json')
    pairs = []
    for slot, task in enumerate(json.loads((args.public / 'tasks.json').read_text())):
        modes = {}
        for mode in ('reference', 'reused_native'):
            folder = args.public / f'results/{slot}/{mode}'
            modes[mode] = dict(timing=json.loads((folder / 'timing.json').read_text()),
                               validation=json.loads((folder / 'validation.json').read_text()))
        assert modes['reference']['validation'] == modes['reused_native']['validation']
        pairs.append(dict(**task, methods=modes, exact=True))
    (args.output / 'public_analysis.json').write_text(json.dumps(pairs, indent=2) + '\n')
    tests = ET.parse(args.public / 'tests_complete.xml')
    assert all(int(s.get('failures', 0)) == int(s.get('errors', 0)) == 0
               for s in tests.iter('testsuite'))
    shutil.copy2(args.public / 'tests_complete.xml', args.output / 'tests.xml')
    root = Path(__file__).resolve().parents[1]
    hashes = {}
    with tarfile.open(args.output / 'provenance.tar.gz', 'w:gz') as archive:
        for name, run in sources.items():
            paths = set(run.glob('*.json')) | set(run.glob('*.out')) | set(run.glob('*.xml')) | set(run.glob('*.py'))
            paths |= set((run / 'inputs').glob('*.json')) | set((run / 'status').glob('*.out'))
            for pattern in ('results/*/*/summary.json', 'results/*/*/timing.json',
                            'results/*/*/validation.json'):
                paths |= set(run.glob(pattern))
            for path in sorted(paths):
                archive.add(path, arcname=f'{name}/{path.relative_to(run)}')
            for folder in ('baseline/src', 'engine/src', 'engine/native', 'engine/bench',
                           'engine/tests', 'engine/tools', 'engine/benchmarks', 'engine/docs'):
                for path in sorted((run / folder).rglob('*')):
                    if path.is_file() and '__pycache__' not in path.parts and '.pytest_cache' not in path.parts:
                        hashes[f'{name}/{path.relative_to(run)}'] = hashlib.sha256(path.read_bytes()).hexdigest()
        for driver in ('conditioned_reuse_pilot.py', 'public_reuse_pilot.py', 'publish_conditioned_reuse.py'):
            archive.add(root / 'bench' / driver, arcname=f'final_drivers/{driver}')
        jobs = '456603,456628,456729,456633,456733,456734,456634,456835,456836,456852,456862'
        accounting = subprocess.check_output(['sacct', '-j', jobs, '-P',
            '--format=JobID,State,ExitCode,Elapsed,TotalCPU,AllocCPUS,MaxRSS,NodeList'], text=True)
        blobs = {'source_sha256.json': json.dumps(hashes, indent=2), 'accounting.txt': accounting,
                 'locations.json': json.dumps({k: str(v) for k, v in sources.items()}, indent=2)}
        for name, content in blobs.items():
            data = content.encode()
            info = tarfile.TarInfo(name); info.size = len(data)
            archive.addfile(info, io.BytesIO(data))
    print('Packaged saved results; no matching performed.', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('pilot', 'refined', 'full', 'public', 'output'):
        parser.add_argument('--' + name, type=Path, required=True)
    main(parser.parse_args())
