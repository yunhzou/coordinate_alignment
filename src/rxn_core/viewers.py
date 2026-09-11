"""Shared molecular and catalog presentation; search code does not depend on it.

The molecular layout is the original white R/P/TS viewer from stable/legacy-aam.
Adapters supply recorded data, rather than introducing report-specific skins.
"""
import html
import json
from pathlib import Path
import re
from functools import wraps

import numpy as np

ASSETS = Path(__file__).with_name('static')


def stylesheet(family='reaction'):
    if family not in ('reaction', 'catalog'):
        raise ValueError('Supported viewer styles are reaction and catalog')
    return (ASSETS / f'{family}_viewer.css').read_text()


def style_document(document, family='reaction', layout='standard'):
    """Apply one registered style to an existing structural layout."""
    document = re.sub(r'<style[^>]*>.*?</style>', '', document, flags=re.S)
    style = f'<style data-viewer-style="{family}">{stylesheet(family)}</style>'
    if '</head>' in document:
        document = document.replace('</head>', style + '</head>', 1)
    else:
        document = re.sub(r'<!doctype[^>]*>', '', document, flags=re.I)
        document = '<!doctype html>' + style + document
    if re.search(r'<body(?:\s[^>]*)?>', document):
        document = re.sub(r' data-viewer-layout="[^"]*"', '', document, count=1)
        document = re.sub(r'<body(?=[\s>])', f'<body data-viewer-layout="{layout}"', document, count=1)
    else:
        document = re.sub(r'(?=<(?:header|main|h1|h2|div|section|article)\b)',
                          '<body data-viewer-layout="' + layout + '">', document, count=1)
    return document


def _json(value):
    return json.dumps(value, separators=(',', ':'), default=lambda x: x.tolist()).replace('</', '<\\/')


def viewer_style(family='reaction', layout='standard'):
    def decorate(function):
        @wraps(function)
        def render(*args, **kwargs):
            return style_document(function(*args, **kwargs), family, layout)
        return render
    return decorate


def align_product(reactant, product):
    """Rigid display fitting only; retain the recorded atom correspondence."""
    r, p = np.asarray(reactant, dtype=float), np.asarray(product, dtype=float)
    r0, p0 = r.mean(axis=0), p.mean(axis=0)
    u, _, vt = np.linalg.svd((p - p0).T @ (r - r0))
    correction = np.eye(3)
    correction[-1, -1] = np.linalg.det(u @ vt)
    return ((p - p0) @ (u @ correction @ vt) + r0).tolist()


def comparison_document(case, *, trajectory_frames=0):
    """Adapt complete saved R/P or endpoint/TS witnesses to the original viewer."""
    r, p = case['endpoints']
    mechanisms = []
    for index, record in enumerate(case['records'], 1):
        source = record['mapping']
        mapping = dict(source) if source and isinstance(source[0], (tuple, list)) else dict(enumerate(source))
        assert sorted(mapping) == sorted(mapping.values()) == list(range(len(r['elements'])))
        assert all(r['elements'][a] == p['elements'][b] for a, b in mapping.items())
        ordered = [p['coordinates'][mapping[a]] for a in range(len(mapping))]
        events = record.get('events', [])
        by_kind = {kind: [e['r'] for e in events if e['kind'] == kind]
                   for kind in ('broken', 'formed', 'strengthened', 'weakened', 'order_changed')}
        mechanisms.append(dict(
            id=index, label=record['name'], mapping_RP=mapping,
            product_xyz_in_R=ordered, product_xyz_in_R_aligned=align_product(r['coordinates'], ordered),
            broken_bonds_R=by_kind['broken'], formed_bonds_R=by_kind['formed'],
            formed_bonds_P=[[mapping[a], mapping[b]] for a, b in by_kind['formed']],
            order_events={k: by_kind[k] for k in ('strengthened', 'weakened', 'order_changed')},
            event_records=events, event_counts=record.get('counts'),
            pattern_id=record.get('pattern'), provenance=record.get('provenance', {}),
            core_atoms=sorted({a for e in events for a in e['r']}), igs=[], gt=None,
        ))
        if trajectory_frames:
            from .alignment.interpolation import internal_coordinate_interpolation
            n = len(mapping)
            reactant_bonds = {(a, b) for a in range(n) for b in range(a + 1, n)
                              if r['wbo'][a][b] > .2}
            product_bonds = {(a, b) for a in range(n) for b in range(a + 1, n)
                             if p['wbo'][mapping[a]][mapping[b]] > .2}
            mechanisms[-1]['endpoint_interpolation'] = internal_coordinate_interpolation(
                r['coordinates'], ordered, r['elements'],
                bonded_pairs=sorted(reactant_bonds | product_bonds),
                persistent_bonded_pairs=sorted(reactant_bonds & product_bonds),
                reactant_bonded_pairs=sorted(reactant_bonds),
                product_bonded_pairs=sorted(product_bonds), n_frames=trajectory_frames)
    return dict(step=str(case['name']), index=case['index'], n_atoms=len(r['elements']),
                default_mech_id=1, include_gt=False, note=case.get('note', ''),
                reactant=dict(elements=r['elements'], coords=r['coordinates'], wbo=r['wbo']),
                product=dict(elements=p['elements'], coords=p['coordinates'], wbo=p['wbo']),
                mechanisms=mechanisms)


def reaction_html(document):
    template = (ASSETS / 'reaction_viewer.html').read_text()
    return (template.replace('__TITLE__', html.escape(document['step']))
            .replace('__VIEWER_STYLE__', stylesheet())
            .replace('__LIBRARY__', (ASSETS / '3Dmol-min.js').read_text())
            .replace('__ZIP_LIBRARY__', (ASSETS / 'jszip.min.js').read_text())
            .replace('__DATA__', _json(document)))


def collection_html(documents, title='R/P and TS comparison'):
    """One copy of the shared renderer and libraries, with a saved-case selector."""
    if not documents:
        raise ValueError('No cases to display')
    template = ((ASSETS / 'reaction_viewer.html').read_text()
                .replace('__VIEWER_STYLE__', stylesheet())
                .replace('__LIBRARY__', (ASSETS / '3Dmol-min.js').read_text())
                .replace('__ZIP_LIBRARY__', (ASSETS / 'jszip.min.js').read_text()))
    page = '''<!doctype html><html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>__PAGE_TITLE__</title></head><body><header class="topbar"><b>__PAGE_TITLE__</b>
<label>Case <select id="case"></select></label></header><iframe id="reaction" title="Original R/P/TS viewer"></iframe>
<script>const CASES=__CASES__,TEMPLATE=__TEMPLATE__,select=document.getElementById('case'),frame=document.getElementById('reaction');
const escape=s=>String(s).replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
CASES.forEach((c,i)=>{const o=document.createElement('option');o.value=i;o.textContent=(c.index===undefined?'':c.index+' · ')+c.step;select.appendChild(o)});
function render(){const c=CASES[+select.value];frame.srcdoc=TEMPLATE.replaceAll('__TITLE__',escape(c.step)).replace('__DATA__',JSON.stringify(c).replaceAll('</','<\\/'));history.replaceState(null,'','#case='+encodeURIComponent(c.index===undefined?select.value:c.index));}
const query=new URLSearchParams(location.hash.slice(1)).get('case'),choice=CASES.findIndex(c=>String(c.index)===query);if(choice>=0)select.value=choice;
select.onchange=render;render();</script></body></html>'''
    page = page.replace('__PAGE_TITLE__', html.escape(title)).replace('__CASES__', _json(documents)).replace('__TEMPLATE__', _json(template))
    # Avoid treating the embedded template's encoded styles as outer styles.
    return page.replace('</head>', '<style data-viewer-style="reaction">' + stylesheet() + '</style></head>', 1).replace('<body>', '<body data-viewer-layout="collection">', 1)
