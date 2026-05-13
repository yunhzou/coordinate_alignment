"""HTML template for the slider-driven alignment trace renderer."""

HTML = r"""<!doctype html>
<html><head><meta charset="utf-8"><title>Island growth: {title}</title>
<script src="https://3dmol.org/build/3Dmol-min.js"></script>
<style>
 html, body {{ margin: 0; padding: 0; }}
 body {{ font-family: -apple-system, sans-serif; background: #fafafa; padding: 10px; box-sizing: border-box; }}
 h2 {{ margin: 4px 0 8px; }}
 .ctl {{ background: white; padding: 8px 12px; border: 1px solid #ddd; border-radius: 6px; margin-bottom: 10px; display: flex; gap: 10px; align-items: center; flex-wrap: wrap; }}
 .row {{ display: flex; gap: 10px; margin-bottom: 10px; }}
 .pane {{ flex: 1; background: white; border: 1px solid #ddd; border-radius: 6px; padding: 8px; box-sizing: border-box; }}
 .viewer-wrap {{ position: relative; width: 100%; height: 480px; }}
 .viewer {{ position: absolute; top: 0; left: 0; right: 0; bottom: 0; }}
 input[type=range] {{ flex: 1; min-width: 200px; }}
 #event {{ font-family: ui-monospace, monospace; font-size: 12px; padding: 8px; background: #f4f4f4; border-radius: 4px; white-space: pre-wrap; max-height: 250px; overflow: auto; }}
 .legend {{ display: flex; gap: 6px; flex-wrap: wrap; }}
 .legend span {{ padding: 2px 8px; border-radius: 4px; font-size: 12px; }}
 .leg-seed {{ background: #ffd700; }}
 .leg-frag {{ background: #66c2a5; }}
 .leg-cut  {{ background: #f88; }}
 .leg-tried {{ background: #ffe066; }}
 .leg-cands-p {{ background: #d4a4ff; }}
 .leg-locked {{ background: #2ca02c; color: white; }}
 button {{ padding: 4px 10px; }}
 h3 {{ margin: 4px 0 6px; font-size: 14px; }}
 h4 {{ margin: 10px 0 4px; font-size: 13px; }}
 pre {{ font-size: 11px; margin: 0; max-height: 220px; overflow: auto; background: #fafafa; padding: 4px 6px; border-radius: 4px; }}
</style></head><body>

<h2>{title}</h2>
<div class="ctl">
  <button id="prev">◀</button>
  <button id="play">▶ Play</button>
  <button id="next">▶</button>
  <input type="range" id="slider" min="0" max="0" value="0">
  <span id="counter">0 / 0</span>
  <label style="font-size:13px"><input type="checkbox" id="showLabels" checked> labels</label>
  <span class="legend" style="margin-left: auto">
    <span class="leg-seed">seed</span>
    <span class="leg-frag">current fragment</span>
    <span style="background:#bbb; padding:2px 8px; border-radius:4px; font-size:12px">locked islands (each = unique color)</span>
  </span>
</div>

<div id="clickInfo" style="background:#fffacc; padding:8px 12px; border:1px solid #ccc; border-radius:6px; margin-bottom:10px; font-family:ui-monospace, monospace; font-size:13px;">
  Click two atoms in the same viewer (R or P) to see the WBO between them. <br>
  Selected R: <span id="selR_display">none</span>  |  Selected P: <span id="selP_display">none</span>
</div>

<div class="row">
  <div class="pane">
    <h3>{left_title}</h3>
    <div class="viewer-wrap"><div id="vR" class="viewer"></div></div>
  </div>
  <div class="pane">
    <h3>{right_title}</h3>
    <div class="viewer-wrap"><div id="vP" class="viewer"></div></div>
  </div>
</div>

<div class="row">
  <div class="pane">
    <h3>Event</h3>
    <div id="event"></div>
  </div>
  <div class="pane">
    <h3>Locked mapping ({left_short} → {right_short})</h3>
    <pre id="locked"></pre>
  </div>
</div>

<div class="row">
  <div class="pane" style="flex:1">
    <h3>All unique alignment mappings (ranked by ranker score S)</h3>
    {mappings_summary_html}
  </div>
</div>

<script>
const xyzR = {xyzR_json};
const xyzP = {xyzP_json};
const events = {events_json};
const wboR = {wboR_json};
const wboP = {wboP_json};
const elementsR = {elements_R_json};
const elementsP = {elements_P_json};

function parseXYZCoords(xyz) {{
  const lines = xyz.trim().split('\n');
  const n = parseInt(lines[0]);
  const out = [];
  for (let i = 0; i < n; i++) {{
    const parts = lines[2+i].trim().split(/\s+/);
    out.push([+parts[1], +parts[2], +parts[3]]);
  }}
  return out;
}}

const coordsR = parseXYZCoords(xyzR);
const coordsP = parseXYZCoords(xyzP);

const vR = $3Dmol.createViewer('vR', {{backgroundColor: 'white'}});
vR.addModel(xyzR, 'xyz');
vR.setStyle({{}}, {{stick: {{radius: 0.10}}, sphere: {{scale: 0.20}}}});

const vP = $3Dmol.createViewer('vP', {{backgroundColor: 'white'}});
vP.addModel(xyzP, 'xyz');
vP.setStyle({{}}, {{stick: {{radius: 0.10}}, sphere: {{scale: 0.20}}}});

// --- Atom labels (toggleable) ---
let labelHandlesR = [];
let labelHandlesP = [];
function addLabels() {{
  labelHandlesR.forEach(h => vR.removeLabel(h));
  labelHandlesP.forEach(h => vP.removeLabel(h));
  labelHandlesR = [];
  labelHandlesP = [];
  // Labels positioned exactly at atom centers, no background, so they
  // don't visually shift to neighboring atoms.
  for (let i = 0; i < coordsR.length; i++) {{
    const h = vR.addLabel(String(i), {{
      position: {{x: coordsR[i][0], y: coordsR[i][1], z: coordsR[i][2]}},
      fontSize: 9, fontColor: 'black',
      showBackground: false,
      borderThickness: 0,
      inFront: true,
    }});
    labelHandlesR.push(h);
  }}
  for (let i = 0; i < coordsP.length; i++) {{
    const h = vP.addLabel(String(i), {{
      position: {{x: coordsP[i][0], y: coordsP[i][1], z: coordsP[i][2]}},
      fontSize: 9, fontColor: 'black',
      showBackground: false,
      borderThickness: 0,
      inFront: true,
    }});
    labelHandlesP.push(h);
  }}
  vR.render(); vP.render();
}}
function removeLabels() {{
  labelHandlesR.forEach(h => vR.removeLabel(h));
  labelHandlesP.forEach(h => vP.removeLabel(h));
  labelHandlesR = []; labelHandlesP = [];
  vR.render(); vP.render();
}}
document.getElementById('showLabels').addEventListener('change', e => {{
  if (e.target.checked) addLabels(); else removeLabels();
}});

vR.zoomTo(); vR.render();
vP.zoomTo(); vP.render();
addLabels();

// --- Click selection for WBO display ---
let selR = [];
let selP = [];
function refreshSelectionDisplay() {{
  function fmt(arr, side) {{
    if (arr.length === 0) return 'none';
    if (arr.length === 1) return `${{side}}[${{arr[0]}}](${{(side==='R'?elementsR:elementsP)[arr[0]]}})`;
    if (arr.length === 2) {{
      const [i, j] = arr;
      const els = side === 'R' ? elementsR : elementsP;
      const wbo = side === 'R' ? wboR[i][j] : wboP[i][j];
      return `${{side}}[${{i}}](${{els[i]}}) ↔ ${{side}}[${{j}}](${{els[j]}})  WBO = ${{wbo.toFixed(3)}}`;
    }}
    return '';
  }}
  document.getElementById('selR_display').textContent = fmt(selR, 'R');
  document.getElementById('selP_display').textContent = fmt(selP, 'P');
}}

vR.setClickable({{}}, true, function(atom, viewer, event, container) {{
  if (event && event.preventDefault) event.preventDefault();
  if (event && event.stopPropagation) event.stopPropagation();
  const idx = atom.serial;
  if (selR.length >= 2) selR = [];
  if (!selR.includes(idx)) selR.push(idx);
  refreshSelectionDisplay();
  applyEvent(+slider.value);
}});
vP.setClickable({{}}, true, function(atom, viewer, event, container) {{
  if (event && event.preventDefault) event.preventDefault();
  if (event && event.stopPropagation) event.stopPropagation();
  const idx = atom.serial;
  if (selP.length >= 2) selP = [];
  if (!selP.includes(idx)) selP.push(idx);
  refreshSelectionDisplay();
  applyEvent(+slider.value);
}});

const slider = document.getElementById('slider');
slider.max = events.length - 1;
const counter = document.getElementById('counter');
const eventDiv = document.getElementById('event');
const lockedPre = document.getElementById('locked');

// Distinct color palette for islands (one color per island id; R-atom and
// P-atom in the same island share the same color so the two viewers stay
// visually synced).
const ISLAND_PALETTE = [
  '#1f77b4','#ff7f0e','#2ca02c','#d62728','#9467bd','#8c564b',
  '#e377c2','#7f7f7f','#bcbd22','#17becf','#aec7e8','#ffbb78',
  '#98df8a','#ff9896','#c5b0d5','#c49c94','#f7b6d2','#c7c7c7',
  '#dbdb8d','#9edae5','#3182bd','#e6550d','#31a354','#756bb1',
  '#636363','#969696','#fd8d3c','#74c476','#fdae6b','#a1d99b',
];
function islandColor(idx) {{ return ISLAND_PALETTE[idx % ISLAND_PALETTE.length]; }}

function applyEvent(idx) {{
  // mapping: r->p, atomIsland: r->island_idx (and inv for p->island_idx)
  let mapping = {{}};
  let atomIslandR = {{}};
  let atomIslandP = {{}};
  let currentSeed = null;
  let fragment = [];
  let pCands = [];
  let lastEvent = null;
  for (let i = 0; i <= idx; i++) {{
    const e = events[i];
    lastEvent = e;
    if (e.type === 'island_locked') {{
      for (const [r, p] of e.pairs) {{
        mapping[r] = p;
        atomIslandR[r] = e.island_idx;
        atomIslandP[p] = e.island_idx;
      }}
      // Re-color any previously-locked atoms that have been merged into
      // this island. `relabeled` lists (r_atom, old_island_id) pairs --
      // they keep their R-P mapping, just change island id (= color).
      if (e.relabeled) {{
        for (const [r, _oldId] of e.relabeled) {{
          atomIslandR[r] = e.island_idx;
          if (mapping[r] !== undefined) atomIslandP[mapping[r]] = e.island_idx;
        }}
      }}
      currentSeed = null; fragment = []; pCands = [];
    }} else if (e.type === 'seed_start') {{
      currentSeed = e.seed;
      fragment = e.fragment ? e.fragment.slice() : [e.seed];
      pCands = e.p_atoms || [];
    }} else if (e.type === 'commit') {{
      fragment = e.fragment.slice();
      pCands = e.p_atoms || pCands;
    }} else if (e.type === 'seed_end') {{
      if (e.result !== 'success') {{
        currentSeed = null; fragment = []; pCands = [];
      }}
    }}
  }}

  // ---- Reactant viewer ----
  vR.setStyle({{}}, {{stick: {{radius: 0.10}}, sphere: {{scale: 0.20}}}});
  // Group locked R atoms by island id and color each group
  const rByIsland = {{}};
  for (const r in atomIslandR) {{
    const k = atomIslandR[r];
    if (!rByIsland[k]) rByIsland[k] = [];
    rByIsland[k].push(+r);
  }}
  for (const k in rByIsland) {{
    const c = islandColor(+k);
    vR.setStyle({{serial: rByIsland[k]}}, {{stick: {{radius: 0.16, color: c}}, sphere: {{scale: 0.30, color: c}}}});
  }}
  if (fragment.length) {{
    const fragSer = fragment.slice();
    vR.setStyle({{serial: fragSer}}, {{stick: {{radius: 0.18, color: '#66c2a5'}}, sphere: {{scale: 0.32, color: '#66c2a5'}}}});
  }}
  if (currentSeed !== null) vR.setStyle({{serial: [currentSeed]}}, {{stick: {{radius: 0.22, color: '#ffd700'}}, sphere: {{scale: 0.44, color: '#ffd700'}}}});
  // Subtle highlight for user-selected atoms: no extra size, just magenta tint
  if (selR.length) {{
    vR.setStyle({{serial: selR.slice()}}, {{stick: {{radius: 0.12, color: '#e91e63'}}, sphere: {{scale: 0.24, color: '#e91e63'}}}});
  }}
  vR.render();

  // ---- Product viewer ----
  vP.setStyle({{}}, {{stick: {{radius: 0.10}}, sphere: {{scale: 0.20}}}});
  const pByIsland = {{}};
  for (const p in atomIslandP) {{
    const k = atomIslandP[p];
    if (!pByIsland[k]) pByIsland[k] = [];
    pByIsland[k].push(+p);
  }}
  for (const k in pByIsland) {{
    const c = islandColor(+k);
    vP.setStyle({{serial: pByIsland[k]}}, {{stick: {{radius: 0.16, color: c}}, sphere: {{scale: 0.30, color: c}}}});
  }}
  if (pCands.length) {{
    const pSer = pCands.slice();
    vP.setStyle({{serial: pSer}}, {{stick: {{radius: 0.18, color: '#66c2a5'}}, sphere: {{scale: 0.32, color: '#66c2a5'}}}});
  }}
  if (selP.length) {{
    vP.setStyle({{serial: selP.slice()}}, {{stick: {{radius: 0.12, color: '#e91e63'}}, sphere: {{scale: 0.24, color: '#e91e63'}}}});
  }}
  vP.render();

  // ---- Event log ---- (verbose diagnostic for alignment events)
  const e = lastEvent;
  let txt = `[${{idx + 1}}/${{events.length}}]  ${{e.type.toUpperCase()}}`;
  if (e.scenario) txt += `  (${{e.scenario}})`;
  txt += `\n`;

  function fmtCand(c) {{
    return '{{' + Object.keys(c).sort((a,b)=>+a-+b).map(k => `R[${{k}}]→P[${{c[k]}}]`).join(', ') + '}}';
  }}
  function fmtPatterns(patterns) {{
    if (!patterns || !patterns.length) return '';
    let out = '';
    for (let i = 0; i < patterns.length; i++) {{
      const p = patterns[i];
      out += `    [${{i}}] witness ${{fmtCand(p.witness || {{}})}}\n`;
      if (p.blocks && p.blocks.length) {{
        for (const b of p.blocks) {{
          out += `        block R[${{b.r_atoms.join(',')}}] ⇄ P[${{b.p_atoms.join(',')}}]`;
          out += `  assignments=${{b.assignments}}  extendable=${{b.extendable}}\n`;
        }}
      }}
    }}
    return out;
  }}
  function fmtHeapTop(h) {{
    if (!h || !h.length) return '  (heap empty)';
    return h.map(x => `    R[${{x.frag_atom}}]→R[${{x.ext_atom}}]  WBO=${{x.wbo}}  [${{x.ext_status}}]`).join('\n');
  }}
  function fmtPool(pool) {{
    if (!pool || !pool.length) return '  (no live edges in pool)';
    let total = pool.reduce((a, b) => a + b.edges.length, 0);
    let out = `  CANDIDATE POOL: ${{pool.length}} fragment atoms have ${{total}} live WBO≥0.2 edges\n`;
    for (const grp of pool) {{
      out += `    R[${{grp.frag_atom}}](${{grp.frag_element}}):  ${{grp.edges.map(ee => `→R[${{ee.ext_atom}}](${{ee.ext_element}})W${{ee.wbo}}[${{ee.ext_status}}]`).join('  ')}}\n`;
    }}
    return out;
  }}

  if (e.type === 'pass_start') {{
    txt += `  pass ${{e.pass}}, mapped so far = ${{e.mapped}}`;
  }} else if (e.type === 'seed_start') {{
    txt += `  seed = R[${{e.seed}}]\n`;
    txt += `  initial cands = ${{e.init_cands}}\n`;
    if (e.represented_assignments !== undefined) txt += `  represented concrete assignments = ${{e.represented_assignments}}\n`;
    txt += `  candidate P-atoms (set): [${{(e.p_atoms||[]).join(', ')}}]`;
    if (e.cand_patterns) txt += `\n  compressed candidate patterns:\n${{fmtPatterns(e.cand_patterns)}}`;
  }} else if (e.type === 'pop') {{
    const ed = e.edge;
    txt += `  POPPED edge: R[${{ed.frag_atom}}] → R[${{ed.ext_atom}}](${{ed.ext_element}})  WBO=${{ed.wbo}}\n`;
    txt += `  scenario: ${{e.scenario}}`;
    if (e.scenario === 'merge_island') txt += `   island_id=${{e.island_id_at_ext}}  →P[${{e.island_image}}]`;
    txt += `\n`;
    const ps = e.pre_state;
    txt += `  fragment (${{ps.fragment_size}} atoms): [${{ps.fragment.join(', ')}}]\n`;
    txt += `  compressed cands count: ${{ps.cands_count}}\n`;
    if (ps.represented_assignments !== undefined) txt += `  represented concrete assignments: ${{ps.represented_assignments}}\n`;
    txt += `  P-atoms claimed across cands: [${{ps.p_atoms_in_cands.join(', ')}}]\n`;
    txt += `  cand sample (top ${{ps.cands_sample.length}}):\n`;
    for (let i = 0; i < ps.cands_sample.length; i++) {{
      txt += `    [${{i}}] ${{fmtCand(ps.cands_sample[i])}}\n`;
    }}
    if (ps.cands_pattern_sample) txt += `  compressed patterns:\n${{fmtPatterns(ps.cands_pattern_sample)}}`;
    txt += `  HEAP next ${{e.heap_top_after_pop.length}}:\n${{fmtHeapTop(e.heap_top_after_pop)}}\n`;
    txt += fmtPool(e.pool_by_frag_atom);
  }} else if (e.type === 'pop_skip') {{
    txt += `  skipped edge R[${{e.edge.frag_atom}}]→R[${{e.edge.ext_atom}}] WBO=${{e.edge.wbo}}: ${{e.reason}}`;
  }} else if (e.type === 'commit') {{
    const ed = e.edge;
    txt += `  + R[${{e.added}}](${{e.element || '?'}}) added to fragment\n`;
    if (ed) txt += `  via edge R[${{ed.frag_atom}}]→R[${{ed.ext_atom}}] WBO=${{ed.wbo}}\n`;
    if (e.scenario === 'merge_island') {{
      txt += `  WHOLE-ISLAND MERGE: absorbed ${{e.island_size_absorbed}} atoms = [${{(e.island_atoms_absorbed||[]).join(', ')}}]\n`;
    }} else if (e.cand_n_value_set) {{
      txt += `  cand[R[${{e.added}}]] possible values: [${{e.cand_n_value_set.join(', ')}}]\n`;
    }}
    if (e.cands_after !== undefined) {{
      txt += `  compressed cands: ${{e.cands_before}} → ${{e.cands_after}}\n`;
      if (e.represented_assignments_after !== undefined) txt += `  represented concrete assignments after = ${{e.represented_assignments_after}}\n`;
    }} else if (e.cands !== undefined) {{
      txt += `  compressed cands = ${{e.cands}}\n`;
    }}
    txt += `  fragment now ${{e.fragment.length}} atoms\n`;
    if (e.cands_sample_after) {{
      txt += `  cand sample after:\n`;
      for (let i = 0; i < e.cands_sample_after.length; i++) {{
        txt += `    [${{i}}] ${{fmtCand(e.cands_sample_after[i])}}\n`;
      }}
    }}
    if (e.cands_pattern_after) txt += `  compressed patterns after:\n${{fmtPatterns(e.cands_pattern_after)}}`;
    if (e.bonds_to_fragment) {{
      txt += `  bonds to fragment: [${{e.bonds_to_fragment.map(([u,w]) => `R[${{u}}]:WBO=${{w}}`).join(', ')}}]\n`;
    }}
    if (e.distance_from_seed !== undefined) {{
      txt += `  distance from seed = ${{e.distance_from_seed}}\n`;
    }}
    if (e.step_info) {{
      const si = e.step_info;
      txt += `  --- step decision: shell=${{si.shell}}, top WBO=${{si.top_wbo}}, cands_before=${{si.cands_before}} ---\n`;
      if (si.tried) for (const t of si.tried) {{
        txt += `    R[${{t.atom}}](${{t.element}}) WBO=${{t.max_wbo_to_frag}} → ${{t.decision}}\n`;
      }}
    }}
    if (e.heap_remaining !== undefined) {{
      txt += `  HEAP remaining = ${{e.heap_remaining}}, next:\n${{fmtHeapTop(e.heap_top)}}\n`;
      txt += fmtPool(e.pool_by_frag_atom);
    }}
  }} else if (e.type === 'consumed') {{
    const ed = e.edge;
    txt += `  CONSUMED edge: R[${{ed.frag_atom}}] → R[${{ed.ext_atom}}](${{ed.ext_element}})  WBO=${{ed.wbo}}\n`;
    txt += `  scenario: ${{e.scenario}}, reason: ${{e.reason}}\n`;
    if (e.island_image !== undefined) {{
      txt += `  island_id=${{e.island_id}}, image=P[${{e.island_image}}]\n`;
    }}
    txt += `  fragment (${{e.fragment.length}} atoms) unchanged, compressed cands (${{e.cands_count}}) unchanged\n`;
    if (e.represented_assignments !== undefined) txt += `  represented concrete assignments = ${{e.represented_assignments}}\n`;
    if (e.cands_pattern_sample) txt += `  compressed patterns:\n${{fmtPatterns(e.cands_pattern_sample)}}`;
    if (e.why_per_cand) {{
      txt += `  WHY EACH CAND FAILED:\n`;
      for (const w of e.why_per_cand) {{
        txt += `    cand[${{w.cand_idx}}]`;
        if (w.variant_idx !== undefined) txt += `.variant[${{w.variant_idx}}]`;
        if (w.variant_multiplicity !== undefined) txt += ` x${{w.variant_multiplicity}}`;
        txt += `: `;
        if (w.cand_at_in_frag_neighbors) {{
          txt += `at neighbors: ${{Object.keys(w.cand_at_in_frag_neighbors).map(k => `R[${{k}}]→P[${{w.cand_at_in_frag_neighbors[k]}}]`).join(', ')}}\n`;
          const vc = w.candidate_v_count !== undefined ? w.candidate_v_count : w.common_v_set_size;
          txt += `             candidate P-target count = ${{vc}}\n`;
          for (const t of w.tried_v) {{
            txt += `               try v=P[${{t.v}}]: ${{t.rejected ? 'REJECT' : 'OK'}} — ${{t.reason}}\n`;
          }}
        }} else if (w.reasons) {{
          txt += `\n`;
          for (const r of w.reasons) txt += `      • ${{r}}\n`;
        }}
      }}
    }}
    txt += `  HEAP remaining = ${{e.heap_remaining}}, next:\n${{fmtHeapTop(e.heap_top)}}\n`;
    txt += fmtPool(e.pool_by_frag_atom);
  }} else if (e.type === 'seed_end') {{
    txt += `  result = ${{e.result}},  final compressed cands = ${{e.final_cands}}`;
    if (e.n_branches !== undefined) txt += `,  branches = ${{e.n_branches}}`;
    txt += `\n`;
    if (e.lock_reason) txt += `  lock reason = ${{e.lock_reason}}\n`;
    txt += `  fragment final size = ${{e.fragment.length}}\n`;
    txt += `  fragment atoms = [${{e.fragment.join(', ')}}]`;
    if (e.all_isos && e.all_isos.length > 1) {{
      txt += `\n  all ${{e.all_isos.length}} unique candidate mappings:`;
      for (let i = 0; i < e.all_isos.length; i++) {{
        txt += `\n    [#${{i+1}}] ${{fmtCand(e.all_isos[i])}}`;
      }}
    }} else if (e.iso) {{
      txt += `\n  iso = ${{fmtCand(e.iso)}}`;
    }}
    if (e.cand_patterns) txt += `\n  final compressed patterns:\n${{fmtPatterns(e.cand_patterns)}}`;
    if (e.heap_remaining !== undefined) txt += `\n  heap remaining at lock = ${{e.heap_remaining}}`;
  }} else if (e.type === 'island_locked') {{
    txt += `  island #${{e.island_idx}} locked\n`;
    txt += `  new pairs (${{e.pairs.length}}): ${{e.pairs.map(p => `R[${{p[0]}}]→P[${{p[1]}}]`).join(', ')}}\n`;
    txt += `  total mapped = ${{e.mapped_total}}`;
    if (e.merged_with && e.merged_with.length) {{
      txt += `\n  merged with islands: [${{e.merged_with.join(', ')}}]`;
    }}
  }} else if (e.type === 'done') {{
    txt += `  final mapped = ${{e.mapped}}`;
  }} else {{
    txt += `  ${{JSON.stringify(e)}}`;
  }}
  eventDiv.textContent = txt;
  counter.textContent = `${{idx + 1}} / ${{events.length}}`;

  let lockedTxt = '';
  for (const r of Object.keys(mapping).sort((a, b) => +a - +b)) lockedTxt += `R[${{r}}] → P[${{mapping[r]}}]\n`;
  lockedPre.textContent = lockedTxt || '(none yet)';
}}

slider.addEventListener('input', () => applyEvent(+slider.value));
document.getElementById('prev').onclick = () => {{ slider.value = Math.max(0, +slider.value - 1); applyEvent(+slider.value); }};
document.getElementById('next').onclick = () => {{ slider.value = Math.min(events.length - 1, +slider.value + 1); applyEvent(+slider.value); }};
let playing = false; let timer = null;
document.getElementById('play').onclick = (e) => {{
  if (playing) {{ clearInterval(timer); playing = false; e.target.textContent = '▶ Play'; }}
  else {{
    playing = true; e.target.textContent = '⏸ Pause';
    timer = setInterval(() => {{
      const cur = +slider.value;
      if (cur >= events.length - 1) {{ clearInterval(timer); playing = false; e.target.textContent = '▶ Play'; return; }}
      slider.value = cur + 1; applyEvent(cur + 1);
    }}, 250);
  }}
}};
window.addEventListener('keydown', e => {{
  if (e.key === 'ArrowRight') {{ slider.value = Math.min(events.length - 1, +slider.value + 1); applyEvent(+slider.value); }}
  if (e.key === 'ArrowLeft') {{ slider.value = Math.max(0, +slider.value - 1); applyEvent(+slider.value); }}
}});
applyEvent(0);
</script>
</body></html>
"""
