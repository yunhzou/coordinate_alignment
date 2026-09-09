"""Offline native-geometry viewer of saved mappings; no search or SMILES conversion."""
import argparse
import json
from pathlib import Path

import numpy as np
from rxn_core import AAMProblem
from rxn_core.domain import MolecularEndpoint
from rxn_core.family_scoring import bond_events


def build(base, output):
    cases = []
    for index in (135, 59, 64):
        row = json.loads((base / 'elementary140_output_comparison_20260908' / f'{index}.refined.json').read_text())
        raw = json.loads((Path(row['aam_source']) / 'inputs' / str(index) / 'input.json').read_text())
        endpoints = [raw[k] for k in ('reactant', 'product')]
        problem = AAMProblem(*(MolecularEndpoint(**e) for e in endpoints))
        records = []
        for name, record in [('Our AAM', row['aam_best_example']),
                             ('SLAP', min(row['slap'], key=lambda c: c['events']['total']))]:
            mapping = dict(record['mapping'])
            assert bond_events(problem, mapping) == record['events']
            assert set(mapping) == set(mapping.values()) == set(range(row['atoms']))
            events = []
            r, p = [np.asarray(e['wbo']) for e in endpoints]
            for a in range(len(r)):
                for b in range(a + 1, len(r)):
                    w, v = float(r[a,b]), float(p[mapping[a],mapping[b]])
                    kind = ('broken' if w > .2 and v <= .2 else
                            'formed' if w <= .2 and v > .2 else
                            'order_changed' if w > .2 and v > .2 and abs(w-v) > .5 else None)
                    if kind: events.append(dict(kind=kind, r=[a,b], p=[mapping[a],mapping[b]], wbo=[w,v]))
            assert len(events) == record['events']['total']
            records.append(dict(name=name, mapping=record['mapping'], events=events,
                                counts=record['events'], provenance=record))
        cases.append(dict(index=index, name=row['name'], endpoints=endpoints, records=records,
                          source=row['aam_source']))
    output.parent.mkdir(parents=True, exist_ok=True)
    library = (Path(__file__).resolve().parents[1] / 'src/rxn_core/static/3Dmol-min.js').read_text()
    output.write_text(HTML.replace('__LIBRARY__', library).replace('__DATA__', json.dumps(cases)))
    print(output.resolve())


HTML = r'''<!doctype html><html><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>AAM / SLAP — saved mapping comparison</title>
<style>
body{font:16px system-ui;margin:0;background:#eef2f6;color:#14283b}header{background:#15324b;color:white;padding:18px 24px}h1{font-size:24px;margin:0 0 8px}p{margin:8px 0}main{padding:16px}button,select{font:inherit;padding:7px}label{margin-right:16px}.controls{padding:12px;background:white;position:sticky;top:0;z-index:5}.row{display:grid;grid-template-columns:1fr 1fr;gap:12px}.panel{background:white;border:1px solid #c6d1df;margin:12px 0;padding:12px;min-width:0}.view{height:390px;position:relative}.heading{font-weight:700}.counts{padding:8px;background:#edf3fa}.note{background:#fff1d6;padding:12px}table{width:100%;border-collapse:collapse;font-size:14px}td,th{padding:6px;border-bottom:1px solid #ddd;text-align:left}tr[data-event]{cursor:pointer}tr[data-event]:hover{background:#e1edff}#inspect{padding:10px;background:#dff0ff}details{background:white;padding:12px;margin-top:12px}#pairs button{border:0;background:transparent;padding:0;color:#1459a1;cursor:pointer}@media(max-width:700px){.view{height:300px}main{padding:6px}.panel{padding:6px}body{font-size:14px}}
</style><header><h1>Our AAM ↔ SLAP: actual saved atom mappings</h1><p>Native XYZ coordinates + original WBO bonds. No rerun. Neither mapping is labeled ground truth.</p></header>
<div class="controls"><select id="case"></select> <label><input id="hydrogen" type="checkbox" checked>Show H</label><label><input id="labels" type="checkbox">Atom labels</label><label><input id="colors" type="checkbox" checked>Source identity colors</label><button id="reset">Reset views</button></div>
<main><div id="summary" class="note"></div><p>Within both methods, the same <b>r-number/color</b> identifies the same input reactant atom. Product atoms carry the assigned r-number. Colors identify atoms, not search fragments; use labels/clicks for exact identity. Element cores remain element-colored. All bonds drawn from WBO &gt; 0.2, with single sticks (not inferred integer bond orders).</p>
<div id="inspect">Click any atom or bond-event row to highlight its source atoms across all four views.</div><div id="panels"></div>
<details><summary>Full atom correspondence (zero-based input indices)</summary><table><thead><tr><th>Source</th><th>Our AAM target</th><th>SLAP target</th></tr></thead><tbody id="pairs"></tbody></table></details>
<details><summary>Provenance and scope</summary><p>Best saved AAM witness versus minimum shared-score SLAP candidate. SLAP hydrogen assignment refinement is included. Bond events count broken, formed and order-changed edges separately; order change means absolute WBO difference &gt; 0.5. Explicit H participates even when hidden. Heavy SLAP patterns shown here were excluded from the saved AAM families by the bounded overlap analysis. This does not prove unreachable under another search policy or establish chemical correctness. Coordinates are unchanged and are not interpolated transition states.</p><pre id="provenance" style="white-space:pre-wrap;overflow-wrap:anywhere"></pre></details></main>
<script>__LIBRARY__</script><script>
const DATA=__DATA__;let current,views=[],selected=[];
const $=id=>document.getElementById(id);
// Convert HSL ourselves: 3Dmol expects hex colors.
function color(i){let h=((i*137.508)%360)/60,c=.65*.84,x=c*(1-Math.abs(h%2-1)),m=.42-c/2,v=h<1?[c,x,0]:h<2?[x,c,0]:h<3?[0,c,x]:h<4?[0,x,c]:h<5?[x,0,c]:[c,0,x];return '#'+v.map(a=>Math.round((a+m)*255).toString(16).padStart(2,'0')).join('')}
function maps(rec){return Object.fromEntries(rec.mapping)}
function draw(){views.forEach(x=>x.v.clear());views=[];$('panels').innerHTML='';selected=[];$('inspect').textContent='Click any atom or bond-event row to highlight its source atoms across all four views.';
 current.records.forEach((rec,k)=>{
 const wrap=document.createElement('section');wrap.className='panel';wrap.innerHTML=`<h2>${rec.name} — ${rec.counts.total} events</h2><div class="counts">Broken ${rec.counts.broken} · Formed ${rec.counts.formed} · Order changed ${rec.counts.order_changed} · All atoms mapped</div><div class="row"><div><p class="heading">R — original reactants</p><div class="view" id="v${k}0"></div></div><div><p class="heading">P — product labeled by mapped R atoms</p><div class="view" id="v${k}1"></div></div></div><table><thead><tr><th>Event</th><th>R atoms → P atoms</th><th>WBO R → P</th></tr></thead><tbody>${rec.events.map((e,j)=>`<tr data-event="${j}"><td>${e.kind.replace('_',' ')}</td><td>r${e.r.join('–r')} → p${e.p.join('–p')}</td><td>${e.wbo.map(x=>x.toFixed(3)).join(' → ')}</td></tr>`).join('')}</tbody></table>`;
 $('panels').appendChild(wrap);wrap.querySelectorAll('[data-event]').forEach(el=>el.onclick=()=>highlight(rec.events[+el.dataset.event].r));
 const inv=Object.fromEntries(rec.mapping.map(([a,b])=>[b,a]));
 current.endpoints.forEach((ep,side)=>{const v=$3Dmol.createViewer($(`v${k}${side}`),{backgroundColor:'white'}),model=v.addModel();
 const atoms=ep.elements.map((elem,i)=>{const c=ep.coordinates[i];return {elem,x:c[0],y:c[1],z:c[2],serial:i,index:i,bonds:ep.wbo[i].flatMap((w,j)=>w>.2&&j!==i?[j]:[]),bondOrder:ep.wbo[i].flatMap((w,j)=>w>.2&&j!==i?[1]:[]),properties:{r:side?inv[i]:i,p:side?i:maps(rec)[i]}}});
 model.addAtoms(atoms);v.setClickable({},true,a=>highlight([a.properties.r]));views.push({v,atoms,k,side});});
 });style();views.forEach(x=>{x.v.zoomTo();x.v.render()});}
function style(){views.forEach(({v,atoms})=>{v.setStyle({},{stick:{radius:.095},sphere:{scale:.21}});v.removeAllLabels();
 atoms.forEach((a,i)=>{const r=a.properties.r,hidden=a.elem==='H'&&!$('hydrogen').checked;
 if(hidden){v.setStyle({index:i},{});return}
 if($('colors').checked)v.addStyle({index:i},{stick:{color:color(r),radius:.1}});
 if(selected.includes(r))v.addStyle({index:i},{sphere:{color:'#ffca00',scale:.38}});
 if($('labels').checked||selected.includes(r))v.addLabel(`${a.elem} r${r}${a.properties.p===undefined?'':' / p'+a.properties.p}`,{position:a,fontSize:12,backgroundColor:'#ffffff',fontColor:'#132e45',backgroundOpacity:.85});
 });v.render()});}
function highlight(rs){selected=rs;style();$('inspect').textContent=rs.map(r=>`r${r} (${current.endpoints[0].elements[r]}): AAM → p${maps(current.records[0])[r]}; SLAP → p${maps(current.records[1])[r]}`).join(' | ')}
function select(){current=DATA[+$('case').value];$('summary').textContent=`Case ${current.index} — ${current.name}. AAM ${current.records[0].counts.total} events; SLAP ${current.records[1].counts.total}. `+(current.index===135?'This is the confirmed lower-event SLAP pattern absent from our saved families.':'Here the absent SLAP pattern has MORE events than our best witness, not an equally good score.');
 $('pairs').innerHTML=current.endpoints[0].elements.map((e,r)=>`<tr><td><button data-r="${r}">${e} r${r}</button></td>${current.records.map(c=>`<td>p${maps(c)[r]}</td>`).join('')}</tr>`).join('');$('pairs').querySelectorAll('button').forEach(b=>b.onclick=()=>highlight([+b.dataset.r]));
 $('provenance').textContent=JSON.stringify({source:current.source,index:current.index,records:current.records.map(r=>({method:r.name,...r.provenance}))},null,2);draw();}
 $('case').innerHTML=DATA.map((r,i)=>`<option value="${i}">Case ${r.index} — ${r.name}</option>`).join('');$('case').onchange=select;['hydrogen','labels','colors'].forEach(id=>$(id).onchange=style);$('reset').onclick=()=>views.forEach(({v})=>{v.zoomTo();v.render()});select();
</script></html>'''


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--base', type=Path, default=Path('/project/yunhengzou/coordinate_alignment/aam_benchmarks'))
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    build(args.base, args.output)
