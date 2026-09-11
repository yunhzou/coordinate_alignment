"""Build a verified 3D algorithm replay from a reusable archive-selection manifest.

Usage: python tools/build_search_trajectory.py manifest.json output_directory
The capture subprocess can load a frozen engine selected by engine_source.
Rendering uses current shared viewer assets and requires no search or replay.
"""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[1]


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('manifest',type=Path)
    p.add_argument('output',type=Path)
    p.add_argument('--capture-only',action='store_true',help=argparse.SUPPRESS)
    a=p.parse_args();a.manifest=a.manifest.resolve();a.output=a.output.resolve()
    manifest=json.loads(a.manifest.read_text())
    resolve=lambda value:str((a.manifest.parent/value).resolve())
    a.output.mkdir(parents=True,exist_ok=True)
    if a.capture_only:
        engine=resolve(manifest['engine_source']) if manifest.get('engine_source') else str(ROOT/'src')
        sys.path.insert(0,engine)
        import rxn_core
        # Load this pipeline module against the selected engine's matcher APIs.
        import importlib.util
        spec=importlib.util.spec_from_file_location('rxn_core.search_trajectory',ROOT/'src/rxn_core/search_trajectory.py')
        module=importlib.util.module_from_spec(spec);sys.modules[spec.name]=module;spec.loader.exec_module(module)
        selections=[{**s,'archive':resolve(s['archive'])} for s in manifest['selections']]
        document=module.build_trajectory(selections,title=manifest.get('title'),key_atoms=manifest.get('key_atoms'),
            watch_targets=manifest.get('watch_targets'),event_tolerance=manifest.get('event_tolerance',.5))
        document['engine_source']=str(Path(rxn_core.__file__).parent)
        document['requested_engine_commit']=manifest.get('engine_commit')
        (a.output/'trace.json').write_text(json.dumps(document,separators=(',',':'),default=lambda x:int(x))+'\n')
        (a.output/'validation.json').write_text(json.dumps(dict(status='passed',checks=document['checks']),indent=2)+'\n')
        print(f'Captured {len(document["checks"])} verified saved paths.',flush=True)
        return
    env={**os.environ,'PYTHONDONTWRITEBYTECODE':'1'}
    subprocess.run([sys.executable,str(Path(__file__).resolve()),str(a.manifest),str(a.output),'--capture-only'],
                   env=env,check=True)
    sys.path.insert(0,str(ROOT/'src'))
    from rxn_core.viewers import growth_trace_html
    (a.output/'algorithm_trajectory.html').write_text(growth_trace_html(json.loads((a.output/'trace.json').read_text())))
    print(a.output/'algorithm_trajectory.html')


if __name__=='__main__':main()
