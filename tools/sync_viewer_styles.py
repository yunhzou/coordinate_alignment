"""Refresh tracked HTML presentation from the two shared styles; retain payloads."""
import gzip
import hashlib
import json
from pathlib import Path
import re
import subprocess

from rxn_core.viewers import style_document

ROOT = Path(__file__).resolve().parents[1]


def sync():
    inventory = json.loads((ROOT / 'docs/viewer-style-migration.json').read_text())
    lookup = inventory['prior_style_layouts']
    changed = []
    names = subprocess.check_output(['git', 'ls-files', '*.html', '*.html.gz'], cwd=ROOT, text=True).splitlines()
    for name in names:
        path = ROOT / name
        if not path.exists() or name.startswith(('bench/', 'src/', 'manuscript/scripts/')):
            continue
        if name == 'tools/anchor_picker.html':
            source = path.read_text()
            if '../src/rxn_core/static/reaction_viewer.css' not in source:
                source = source.replace('</head>', '<link rel="stylesheet" href="../src/rxn_core/static/reaction_viewer.css"></head>', 1)
            source = source.replace('<body>', '<body data-viewer-layout="anchor">', 1)
            path.write_text(source)
            continue
        source = gzip.decompress(path.read_bytes()).decode() if name.endswith('.gz') else path.read_text()
        shared = re.search(r'<style data-viewer-style="([^"]+)"', source)
        body = re.search(r'<body data-viewer-layout="([^"]+)"', source)
        if body and body[1] == 'collection':
            # The renderer also refreshes the iframe template embedded in JavaScript.
            from render_mapping_comparison import render
            render(path.with_name('viewer_data.json'), path)
            continue
        layouts = [lookup[h] for css in re.findall(r'<style[^>]*>(.*?)</style>', source, re.S)
                   if (h := hashlib.sha256(css.encode()).hexdigest()) in lookup]
        if shared and body:
            layouts = [(shared[1], body[1])]
        if not layouts:
            continue
        family, layout = layouts[0]
        result = style_document(source, family, layout)
        if result == source:
            continue
        # Presentation updates must leave script bodies and molecular drawing data intact.
        assert re.findall(r'<script[^>]*>(.*?)</script>', source, re.S) == re.findall(r'<script[^>]*>(.*?)</script>', result, re.S)
        if name.endswith('.gz'):
            path.write_bytes(gzip.compress(result.encode(), mtime=0))
        else:
            path.write_text(result)
        changed.append(dict(path=name, family=family, layout=layout,
                            previous_sha256=hashlib.sha256(source.encode()).hexdigest(),
                            current_sha256=hashlib.sha256(result.encode()).hexdigest()))
    record = ROOT / 'docs/viewer-style-output-migration.json'
    previous = {row['path']: row for row in json.loads(record.read_text())} if record.exists() else {}
    for row in changed:
        if row['path'] in previous:
            row['previous_sha256'] = previous[row['path']]['previous_sha256']
        previous[row['path']] = row
    record.write_text(json.dumps(list(previous.values()), indent=2) + '\n')
    print(f'Restyled {len(changed)} tracked pages without changing embedded scripts or molecular data.')


if __name__ == '__main__':
    sync()
