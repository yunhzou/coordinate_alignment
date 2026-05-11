"""Generate R↔P bond breaking/forming view for pr14, both mechanisms.

Re-runs the parallel cut_sweep, picks the 2 min-bond mechanisms,
renders a 3Dmol HTML showing R + P side-by-side with broken bonds
(red dashed) and formed bonds (green dashed) per mechanism."""
from __future__ import annotations
import sys, json, time
sys.path.insert(0, 'src')
import multiprocessing as mp
import random
from pathlib import Path
import numpy as np

from rxn_core import parse_xyz, classify_bonds
from rxn_core.pq import build_graph, find_islands_pq, expand_mapping

WORK = Path('appendix_perparation/xtb_frequency_calculations')
STEP = 'pr14.Pd_hydroamination_JOC2025_TS3_step2_alkene_inserion'
WBO_STRONG = 0.5
N_SEEDS_PER_CUT = 3


def load_step(d):
    xyz_path = next(p for p in d.glob('*.xyz') if 'xtbhess' not in p.name)
    el, xyz = parse_xyz(xyz_path)
    n = len(el); wbo = np.zeros((n, n))
    for ln in (d / 'wbo').read_text().splitlines():
        p = ln.split()
        if len(p) < 3: continue
        i, j = int(p[0])-1, int(p[1])-1
        wbo[i, j] = float(p[2]); wbo[j, i] = wbo[i, j]
    return el, np.asarray(xyz, float), wbo


_W = {}
def _winit(elR, wboR, elT, wboT):
    _W['elR'] = elR; _W['wboR'] = wboR
    _W['elT'] = elT; _W['wboT'] = wboT
    _W['g_P'] = build_graph(elT, wboT, bond_cut=0.2)
    _W['n'] = len(elR)


def _wrun(args):
    cut, order = args
    g_R = build_graph(_W['elR'], _W['wboR'], bond_cut=0.2)
    for (i, j) in cut:
        if g_R.has_edge(i, j): g_R.remove_edge(i, j)
    try:
        branches = find_islands_pq(g_R, _W['g_P'], list(order))
    except Exception:
        return []
    out = []
    for b in branches:
        mapping = expand_mapping(b.mapping, g_R, _W['g_P'])
        if len(mapping) < _W['n'] - 2: continue
        broken, formed, _, _ = classify_bonds(mapping, _W['wboR'], _W['wboT'])
        inv = {v: k for k, v in mapping.items()}
        br = tuple(sorted((min(a, b), max(a, b)) for (a, b, _, _) in broken))
        fm_R = tuple(sorted((min(inv.get(a, -1), inv.get(b, -1)),
                              max(inv.get(a, -1), inv.get(b, -1)))
                             for (a, b, _, _) in formed if a in inv and b in inv))
        # formed in P-frame
        fm_P = tuple(sorted((min(a, b), max(a, b)) for (a, b, _, _) in formed))
        out.append(((br, fm_R), tuple(sorted(mapping.items())), cut, fm_P))
    return out


def parallel_cut_sweep(elR, wboR, elT, wboT, n_workers=14):
    n = len(elR)
    strong = [(i, j) for i in range(n) for j in range(i+1, n)
              if wboR[i, j] >= WBO_STRONG]
    g_R = build_graph(elR, wboR, bond_cut=0.2)
    nodes = list(g_R.nodes())
    rng = random.Random(42)
    seed_orders = []
    for s in range(N_SEEDS_PER_CUT):
        order = list(nodes); rng.shuffle(order)
        seed_orders.append(tuple(order))
    work = []
    for s in seed_orders: work.append(((), s))
    for (i, j) in strong:
        for s in seed_orders: work.append((((i, j),), s))
    pool_chems = {}
    with mp.Pool(n_workers, initializer=_winit,
                  initargs=(elR, wboR, elT, wboT)) as pool:
        for results in pool.imap_unordered(_wrun, work, chunksize=4):
            for chem_sig, mapping_items, cut, fm_P in results:
                pool_chems.setdefault(chem_sig, {
                    'mapping': dict(mapping_items),
                    'cut': cut,
                    'fm_P': fm_P,
                })
    return pool_chems


HTML = """<!doctype html><html><head><meta charset="utf-8"><title>{title}</title>
<script src="https://3dmol.org/build/3Dmol-min.js"></script>
<style>
html,body{{margin:0;padding:0;font-family:-apple-system,sans-serif;background:#fafafa;color:#222}}
body{{padding:14px;box-sizing:border-box}}
h2{{margin:0 0 6px;font-size:20px}}
.note{{font-size:13px;color:#555;margin-bottom:12px}}
.mech-sel{{background:white;padding:10px;border:1px solid #ccc;border-radius:6px;margin-bottom:14px}}
.mech-sel button{{padding:7px 14px;margin-right:8px;border:1px solid #aaa;background:#f0f0f0;border-radius:4px;cursor:pointer;font-family:ui-monospace,monospace;font-size:12px}}
.mech-sel button.active{{background:#ffd700;border-color:#a90;font-weight:600}}
.row{{display:grid;grid-template-columns:1fr 1fr;gap:14px}}
.panel{{background:white;border:1px solid #ddd;border-radius:6px;padding:8px}}
.ph{{display:flex;justify-content:space-between;align-items:baseline;font-size:13px;margin-bottom:6px}}
.ph .lbl{{font-weight:600;font-size:14px}}
.vw{{position:relative;width:100%;height:520px}}
.vwbox{{position:absolute;inset:0}}
.meta{{font-family:ui-monospace,monospace;font-size:12px;color:#444;padding:6px 0 0;line-height:1.5;border-top:1px solid #eee;margin-top:6px}}
.bond{{display:inline-block;padding:1px 6px;border-radius:3px;font-family:ui-monospace,monospace;font-size:11px;margin-right:4px;margin-top:2px}}
.bond.br{{background:#fee;color:#a00;border:1px solid #faa}}
.bond.fm{{background:#efe;color:#080;border:1px solid #afa}}
</style></head><body>
<h2>{title}</h2>
<div class="note">R atoms numbered by R-frame index. Broken bonds (red dashed) drawn on R panel; formed bonds (green dashed) drawn on P panel using product atom indices. Two mechanisms found — both at br/fm=2/3.</div>
<div class="mech-sel" id="mech-sel"></div>
<div class="row">
  <div class="panel">
    <div class="ph"><span class="lbl">Reactant (R)</span><span style="color:#a00;font-family:ui-monospace,monospace;font-size:11px">broken bonds dashed red</span></div>
    <div class="vw"><div id="vw_R" class="vwbox"></div></div>
    <div class="meta" id="meta_R"></div>
  </div>
  <div class="panel">
    <div class="ph"><span class="lbl">Product (P)</span><span style="color:#080;font-family:ui-monospace,monospace;font-size:11px">formed bonds dashed green</span></div>
    <div class="vw"><div id="vw_P" class="vwbox"></div></div>
    <div class="meta" id="meta_P"></div>
  </div>
</div>
<script>
const DATA = {data_json};
let curMech = 0;
function xyzBody(els, xyz) {{
  let s = els.length + "\\nframe\\n";
  for (let i = 0; i < xyz.length; i++) s += els[i] + " " + xyz[i][0].toFixed(6) + " " + xyz[i][1].toFixed(6) + " " + xyz[i][2].toFixed(6) + "\\n";
  return s;
}}
function drawBonds(v, xyz, pairs, color) {{
  for (const [i, j] of pairs) {{
    if (i >= xyz.length || j >= xyz.length) continue;
    v.addCylinder({{start:{{x:xyz[i][0],y:xyz[i][1],z:xyz[i][2]}}, end:{{x:xyz[j][0],y:xyz[j][1],z:xyz[j][2]}}, color:color, radius:0.16, dashed:true}});
  }}
}}
function makeView(divId, els, xyz, bonds, color, label_atoms) {{
  document.getElementById(divId).innerHTML = "";
  const v = $3Dmol.createViewer(divId, {{backgroundColor:'white'}});
  v.addModel(xyzBody(els, xyz), 'xyz');
  v.setStyle({{}}, {{stick:{{radius:0.10}}, sphere:{{scale:0.20}}}});
  // label only atoms in the bonds list
  const lab = new Set();
  for (const [i,j] of bonds) {{ lab.add(i); lab.add(j); }}
  for (const i of lab) {{
    v.addLabel(String(i), {{
      position: {{x: xyz[i][0], y: xyz[i][1], z: xyz[i][2]}},
      backgroundColor: 'white', backgroundOpacity: 0.6,
      fontColor: '#333', fontSize: 11, borderThickness: 0
    }});
  }}
  drawBonds(v, xyz, bonds, color);
  v.zoomTo();
  v.render();
  return v;
}}
function fmtBonds(els, pairs, klass) {{
  return pairs.map(([i,j]) => `<span class="bond ${{klass}}">${{els[i]}}${{i}}-${{els[j]}}${{j}}</span>`).join('');
}}
function render() {{
  const m = DATA.mechs[curMech];
  document.querySelectorAll('.mech-sel button').forEach((b, idx) => b.classList.toggle('active', idx === curMech));
  makeView('vw_R', DATA.elements_R, DATA.xyzR, m.broken_R, 'red');
  makeView('vw_P', DATA.elements_P, DATA.xyzP, m.formed_P, 'green');
  document.getElementById('meta_R').innerHTML = "<b>broken bonds (R-frame):</b> " + fmtBonds(DATA.elements_R, m.broken_R, 'br');
  document.getElementById('meta_P').innerHTML = "<b>formed bonds (P-frame):</b> " + fmtBonds(DATA.elements_P, m.formed_P, 'fm');
}}
const ms = document.getElementById('mech-sel');
ms.innerHTML = "<span style='font-size:13px;margin-right:8px;color:#444'>Mechanism:</span>";
DATA.mechs.forEach((m, idx) => {{
  const b = document.createElement('button');
  b.textContent = `#${{idx+1}}  br/fm=${{m.broken_R.length}}/${{m.formed_P.length}}  (via cut ${{m.cut_label}})`;
  b.onclick = () => {{ curMech = idx; render(); }};
  ms.appendChild(b);
}});
window.addEventListener('load', render);
</script>
</body></html>
"""


def main():
    sd = WORK / STEP
    print(f'loading {STEP}...', flush=True)
    elR, xyzR, wboR = load_step(sd / 'R')
    elP, xyzP, wboP = load_step(sd / 'P')
    print(f'  {len(elR)} atoms', flush=True)

    print('running parallel cut_sweep (14 workers)...', flush=True)
    t0 = time.time()
    pool = parallel_cut_sweep(elR, wboR, elP, wboP, n_workers=14)
    print(f'  {len(pool)} chem classes in {time.time()-t0:.1f}s', flush=True)
    mn = min(len(k[0]) + len(k[1]) for k in pool)
    mn_mechs = [(k, v) for k, v in pool.items() if len(k[0]) + len(k[1]) == mn]
    print(f'  {len(mn_mechs)} mechanisms at min br+fm={mn}')

    mechs_for_html = []
    for (br, fm_R), info in mn_mechs:
        cut = info['cut']
        cut_label = 'none' if not cut else ','.join(f'{elR[a]}{a}-{elR[b]}{b}' for a, b in cut)
        mechs_for_html.append({
            'broken_R': [list(p) for p in br],
            'formed_R': [list(p) for p in fm_R],
            'formed_P': [list(p) for p in info['fm_P']],
            'cut_label': cut_label,
        })
        print(f'  br/fm={len(br)}/{len(fm_R)}  cut=[{cut_label}]')

    data = {
        'elements_R': elR,
        'xyzR': np.asarray(xyzR).tolist(),
        'elements_P': elP,
        'xyzP': np.asarray(xyzP).tolist(),
        'mechs': mechs_for_html,
    }
    out_dir = Path('out/bgcp_views') / STEP
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / 'RP_mechs.html'
    out_path.write_text(HTML.format(
        title=f'pr14 R↔P mechanism view ({len(mn_mechs)} mechanisms at br/fm=2/3)',
        data_json=json.dumps(data),
    ))
    print(f'\nwrote {out_path}', flush=True)
    print(f'view it: file://{out_path.resolve()}')


if __name__ == '__main__':
    main()
