"""
Run grow_island with trace logging on a single step, then build an HTML
animation showing each seed attempt, fragment growth, cand count, and
cut/commit decisions.
"""

from __future__ import annotations
import json
import sys
from pathlib import Path

from rxn_core_frag import (
    run_xtb, build_graph, expand_mapping, classify_bonds,
    write_xyz_str,
)


def _p_atoms_from_cands(cands):
    """Return the union of P atoms across ALL candidate isos -- so the
    user can see every place in the product where the current R fragment
    could potentially be mapped. As the fragment grows and cands narrows,
    this set shrinks. When cands == 1, this equals the unique iso's
    image."""
    out = set()
    for c in cands:
        for v in c.values():
            out.add(int(v))
    return sorted(out)


def grow_island_traced(g_R, g_P, seed, mapping, inv,
                       wbo_tol=0.5, growth_min_wbo=0.6, top_degen=0.1,
                       max_lock_cands=100, min_lock_size=2):
    max_cands_hard = 2000
    max_iters = g_R.number_of_nodes()
    """Same logic as grow_island, but records every decision into `events`."""
    events = []
    if seed in mapping:
        return None, events
    seed_el = g_R.nodes[seed]['element']
    candidates = [{seed: v} for v in g_P.nodes()
                  if v not in inv and g_P.nodes[v]['element'] == seed_el]
    if not candidates:
        events.append({'type': 'seed_start', 'seed': seed, 'init_cands': 0,
                       'p_atoms': []})
        events.append({'type': 'seed_end', 'result': 'no_initial_cands'})
        return None, events
    fragment = {seed}
    distance = {seed: 0}
    events.append({
        'type': 'seed_start',
        'seed': seed,
        'init_cands': len(candidates),
        'fragment': sorted(fragment),
        'p_atoms': _p_atoms_from_cands(candidates),
    })
    for _ in range(max_iters):
        if len(candidates) == 1:
            events.append({
                'type': 'seed_end', 'result': 'success',
                'final_cands': len(candidates),
                'fragment': sorted(fragment),
                'iso': {int(k): int(v) for k, v in candidates[0].items()},
            })
            return candidates[0], events
        frontier = set()
        for u in fragment:
            for n in g_R.neighbors(u):
                if n not in fragment:
                    frontier.add(n)
        if not frontier:
            ok = (bool(candidates) and len(candidates) <= max_lock_cands
                  and len(fragment) >= min_lock_size)
            events.append({
                'type': 'seed_end',
                'result': 'success' if ok else 'no_frontier',
                'final_cands': len(candidates),
                'fragment': sorted(fragment),
                'iso': {int(k): int(v) for k, v in candidates[0].items()} if ok else None,
            })
            return (candidates[0] if ok else None), events

        # Compute frontier info: distance from seed + max WBO; track WHY
        # each frontier atom was filtered (so the trace can explain it).
        frontier_info = {}
        filter_reason = {}  # atom -> reason string
        for n in frontier:
            bonded = [u for u in g_R.neighbors(n) if u in fragment]
            max_w = max(g_R[u][n]['wbo'] for u in bonded)
            dist = 1 + min(distance[u] for u in bonded if u in distance)
            if max_w < growth_min_wbo:
                filter_reason[n] = (f'wbo<{growth_min_wbo}', max_w, dist)
                continue
            frontier_info[n] = (dist, max_w)
        if not frontier_info:
            ok = (bool(candidates) and len(candidates) <= max_lock_cands
                  and len(fragment) >= min_lock_size)
            events.append({
                'type': 'seed_end',
                'result': 'success' if ok else 'no_strong_frontier',
                'final_cands': len(candidates),
                'fragment': sorted(fragment),
                'iso': {int(k): int(v) for k, v in candidates[0].items()} if ok else None,
            })
            return (candidates[0] if ok else None), events
        # BFS by distance, then top-WBO within shell
        min_dist = min(d for (d, _) in frontier_info.values())
        for n, (d, w) in frontier_info.items():
            if d != min_dist:
                filter_reason[n] = (f'shell={d}>min_shell={min_dist}', w, d)
        same_shell = {n: w for n, (d, w) in frontier_info.items() if d == min_dist}
        top_w = max(same_shell.values())
        for n, w in same_shell.items():
            if w < top_w - top_degen:
                filter_reason[n] = (f'wbo {w:.2f} < top {top_w:.2f} - {top_degen}',
                                    w, min_dist)
        strong_frontier = [n for n, w in same_shell.items() if w >= top_w - top_degen]

        best_n = None
        best_cands = None
        tries = []
        for n in strong_frontier:
            n_el = g_R.nodes[n]['element']
            bonded = [u for u in g_R.neighbors(n) if u in fragment]
            r_wbos = [(u, g_R[u][n]['wbo']) for u in bonded]
            n_pinned = mapping.get(n, None)
            new_cands = []
            over = False
            for cand in candidates:
                used_p = set(cand.values())
                v_set = set(g_P.neighbors(cand[bonded[0]]))
                for u in bonded[1:]:
                    v_set &= set(g_P.neighbors(cand[u]))
                v_set -= used_p
                if n_pinned is not None:
                    v_set = v_set & {n_pinned}
                for v in v_set:
                    if g_P.nodes[v]['element'] != n_el:
                        continue
                    if all(abs(w - g_P[cand[u]][v]['wbo']) <= wbo_tol
                           for u, w in r_wbos):
                        nc = dict(cand); nc[n] = v
                        new_cands.append(nc)
                        if len(new_cands) > max_cands_hard:
                            over = True; break
                if over: break
            decision = 'CUT' if (over or not new_cands) else 'ok'
            # Connecting bond WBO(s) and distance from seed
            wbo_str = ', '.join(f'{u}:{w:.2f}' for u, w in r_wbos)
            max_w = max(w for _, w in r_wbos) if r_wbos else 0.0
            min_d = 1 + min(distance[u] for u, _ in r_wbos if u in distance)
            tries.append({
                'atom': int(n),
                'element': g_R.nodes[n]['element'],
                'new_cands': len(new_cands),
                'over': over,
                'decision': decision,
                'max_wbo_to_frag': round(max_w, 3),
                'wbo_bonds': wbo_str,
                'distance_from_seed': min_d,
            })
            if not new_cands or over:
                continue
            if best_cands is None or len(new_cands) < len(best_cands):
                best_n = n
                best_cands = new_cands
        # NO tries event -- user wants only commit/seed frames.
        # Capture the tried/filtered metadata to attach to the next commit.
        filtered = []
        for n, (reason, w, d) in filter_reason.items():
            filtered.append({
                'atom': int(n),
                'element': g_R.nodes[n]['element'],
                'max_wbo_to_frag': round(w, 3),
                'distance_from_seed': d,
                'filtered_reason': reason,
            })
        _step_info = {
            'cands_before': len(candidates),
            'shell': min_dist,
            'top_wbo': round(top_w, 3),
            'tried': tries,
            'filtered': filtered,
        }
        if best_n is None:
            ok = (bool(candidates) and len(candidates) <= max_lock_cands
                  and len(fragment) >= min_lock_size)
            events.append({
                'type': 'seed_end',
                'result': 'success' if ok else 'all_cut',
                'final_cands': len(candidates),
                'fragment': sorted(fragment),
                'iso': {int(k): int(v) for k, v in candidates[0].items()} if ok else None,
            })
            return (candidates[0] if ok else None), events
        fragment.add(best_n)
        bonded_in_frag = [u for u in g_R.neighbors(best_n) if u in fragment - {best_n}]
        distance[best_n] = 1 + min(distance[u] for u in bonded_in_frag)
        commit_bonds = [(u, round(g_R[u][best_n]['wbo'], 3)) for u in bonded_in_frag]
        candidates = best_cands
        events.append({
            'type': 'commit',
            'added': int(best_n),
            'element': g_R.nodes[best_n]['element'],
            'cands': len(candidates),
            'fragment': sorted(fragment),
            'p_atoms': _p_atoms_from_cands(candidates),
            'distance_from_seed': distance[best_n],
            'bonds_to_fragment': commit_bonds,
            # Step metadata: what was considered before this commit
            'step_info': _step_info,
        })
    ok = (bool(candidates) and len(candidates) <= max_lock_cands
          and len(fragment) >= min_lock_size)
    events.append({
        'type': 'seed_end',
        'result': 'success' if ok else 'max_iters',
        'final_cands': len(candidates),
        'fragment': sorted(fragment),
        'iso': {int(k): int(v) for k, v in candidates[0].items()} if ok else None,
    })
    return (candidates[0] if ok else None), events


def find_islands_traced(g_R, g_P, wbo_tol=0.5):
    """STRICTLY SEQUENTIAL trace: walk seed atoms in order, try each one
    exactly once, fully propagate or abort, move to next. No multi-pass
    retry, no parallel seed exploration."""
    mapping = {}
    inv = {}
    atom_island_R = {}
    atom_island_P = {}
    all_events = []
    n_islands_total = 0
    all_events.append({'type': 'pass_start', 'pass': 1, 'mapped': 0})
    if True:  # single pass
        for seed in sorted(g_R.nodes()):
            if seed in mapping:
                continue
            iso, events = grow_island_traced(
                g_R, g_P, seed, mapping, inv, wbo_tol=wbo_tol)
            all_events.extend(events)
            if iso is None:
                continue

            # Determine which existing islands this iso touches
            touched = set()
            for r in iso.keys():
                if r in atom_island_R:
                    touched.add(atom_island_R[r])
            # Pick a canonical id: smallest touched, or new id if none
            if touched:
                merged_id = min(touched)
            else:
                n_islands_total += 1
                merged_id = n_islands_total

            # Commit new atoms to the merged_id; relabel old atoms in touched
            committed_new = []
            relabeled = []
            for r, p in iso.items():
                if r not in mapping:
                    mapping[r] = p
                    inv[p] = r
                    atom_island_R[r] = merged_id
                    atom_island_P[p] = merged_id
                    committed_new.append((int(r), int(p)))
                else:
                    # Relabel if this atom was in a different (touched) island
                    if atom_island_R.get(r) != merged_id:
                        relabeled.append((int(r), int(atom_island_R[r])))
                        atom_island_R[r] = merged_id
                        atom_island_P[mapping[r]] = merged_id
            # Also relabel any other atoms whose island_id is in `touched`
            for r in list(atom_island_R.keys()):
                if atom_island_R[r] in touched and atom_island_R[r] != merged_id:
                    relabeled.append((int(r), int(atom_island_R[r])))
                    atom_island_R[r] = merged_id
                    atom_island_P[mapping[r]] = merged_id

            all_events.append({
                'type': 'island_locked',
                'island_idx': merged_id,
                'pairs': committed_new,
                'merged_with': sorted(touched - {merged_id}),
                'relabeled': relabeled,
                'mapped_total': len(mapping),
            })
        # end single sequential pass

    # Phase 2: expand_mapping -- pair up symmetric same-element neighbors
    # of mapped atoms (e.g. the three methyl Hs that strict fragment
    # growth refused to lock alone). The new atoms join the PARENT atom's
    # island (no new color), since they're chemically equivalent extensions
    # of an existing locked atom.
    from collections import defaultdict as _dd
    progressed = True
    pass_no_exp = 0
    while progressed:
        progressed = False
        pass_no_exp += 1
        for u in list(mapping.keys()):
            v = mapping[u]
            r_groups = _dd(list)
            for w in g_R.neighbors(u):
                if w in mapping:
                    continue
                r_groups[g_R.nodes[w]['element']].append(w)
            p_groups = _dd(list)
            for x in g_P.neighbors(v):
                if x in inv:
                    continue
                p_groups[g_P.nodes[x]['element']].append(x)
            for el, rs in r_groups.items():
                ps = p_groups.get(el, [])
                if len(ps) != len(rs):
                    continue
                # Paired atoms join u's island
                parent_island = atom_island_R.get(u)
                if parent_island is None:
                    n_islands_total += 1
                    parent_island = n_islands_total
                    atom_island_R[u] = parent_island
                    atom_island_P[v] = parent_island
                paired = []
                for w, x in zip(rs, ps):
                    mapping[w] = x
                    inv[x] = w
                    atom_island_R[w] = parent_island
                    atom_island_P[x] = parent_island
                    paired.append((int(w), int(x)))
                if paired:
                    all_events.append({
                        'type': 'island_locked',
                        'island_idx': parent_island,
                        'pairs': paired,
                        'merged_with': [],
                        'relabeled': [],
                        'mapped_total': len(mapping),
                        'expand_pass': pass_no_exp,
                        'parent_atom': int(u),
                    })
                    progressed = True
    # Phase 3: explicit ISLAND-TO-ISLAND MERGE pass.
    # For every pair of distinct islands that share a direct edge in g_R,
    # check whether the cross-island bond(s) are valid in g_P (with WBO
    # match). If so, merge the two islands into one (relabel atoms).
    while True:
        merged_this_round = False
        seen = set()
        for u, v in g_R.edges():
            iA = atom_island_R.get(u); iB = atom_island_R.get(v)
            if iA is None or iB is None or iA == iB:
                continue
            pair = (min(iA, iB), max(iA, iB))
            if pair in seen: continue
            seen.add(pair)
            ok = True
            cross_bonds = []
            for x, y in g_R.edges():
                ixA = atom_island_R.get(x); ixB = atom_island_R.get(y)
                if {ixA, ixB} != {iA, iB}: continue
                wR = g_R[x][y]['wbo']
                px, py = mapping[x], mapping[y]
                if not g_P.has_edge(px, py):
                    ok = False; break
                wP = g_P[px][py]['wbo']
                if abs(wR - wP) > wbo_tol:
                    ok = False; break
                cross_bonds.append((int(x), int(y), round(wR, 3), round(wP, 3)))
            if not ok:
                continue
            keep, drop = pair
            relabeled = []
            for r, idx in list(atom_island_R.items()):
                if idx == drop:
                    relabeled.append((int(r), int(drop)))
                    atom_island_R[r] = keep
                    atom_island_P[mapping[r]] = keep
            all_events.append({
                'type': 'island_locked',
                'island_idx': keep,
                'pairs': [],
                'merged_with': [drop],
                'relabeled': relabeled,
                'mapped_total': len(mapping),
                'island_island_merge': True,
                'cross_bonds': cross_bonds,
            })
            merged_this_round = True
            break
        if not merged_this_round:
            break

    all_events.append({'type': 'done', 'mapped': len(mapping)})
    return mapping, all_events


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
    <h3>Reactant</h3>
    <div class="viewer-wrap"><div id="vR" class="viewer"></div></div>
  </div>
  <div class="pane">
    <h3>Product</h3>
    <div class="viewer-wrap"><div id="vP" class="viewer"></div></div>
  </div>
</div>

<div class="row">
  <div class="pane">
    <h3>Event</h3>
    <div id="event"></div>
  </div>
  <div class="pane">
    <h3>Locked mapping (R → P)</h3>
    <pre id="locked"></pre>
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

  // ---- Event log ---- (full metadata)
  let txt = `[${{idx + 1}}/${{events.length}}] ${{lastEvent.type}}\n`;
  if (lastEvent.type === 'seed_start') {{
    txt += `  seed = R[${{lastEvent.seed}}]\n`;
    txt += `  initial cands = ${{lastEvent.init_cands}}\n`;
    txt += `  fragment = [${{lastEvent.fragment.join(', ')}}]\n`;
    txt += `  first iso P-atoms = [${{(lastEvent.p_atoms||[]).join(', ')}}]`;
  }} else if (lastEvent.type === 'commit') {{
    txt += `  + R[${{lastEvent.added}}](${{lastEvent.element || '?'}}),  cands=${{lastEvent.cands}},  fragment_size=${{lastEvent.fragment.length}}\n`;
    if (lastEvent.bonds_to_fragment) {{
      txt += `  bonds to fragment: [${{lastEvent.bonds_to_fragment.map(([u,w]) => `R[${{u}}]:WBO=${{w}}`).join(', ')}}]\n`;
    }}
    if (lastEvent.distance_from_seed !== undefined) {{
      txt += `  distance from seed = ${{lastEvent.distance_from_seed}}\n`;
    }}
    const si = lastEvent.step_info;
    if (si) {{
      txt += `  --- step decision: shell=${{si.shell}},  top WBO=${{si.top_wbo}},  cands_before=${{si.cands_before}} ---\n`;
      txt += `  TRIED (${{si.tried.length}}):\n`;
      for (const t of si.tried) {{
        txt += `    R[${{t.atom}}](${{t.element}})  WBO=${{t.max_wbo_to_frag}}  dist=${{t.distance_from_seed}}  bonds=[${{t.wbo_bonds}}]  new_cands=${{t.new_cands}} → ${{t.decision}}\n`;
      }}
      if (si.filtered && si.filtered.length) {{
        txt += `  FILTERED (${{si.filtered.length}}):\n`;
        for (const f of si.filtered) {{
          txt += `    R[${{f.atom}}](${{f.element}})  WBO=${{f.max_wbo_to_frag}}  dist=${{f.distance_from_seed}}  reason: ${{f.filtered_reason}}\n`;
        }}
      }}
    }}
  }} else if (lastEvent.type === 'seed_end') {{
    txt += `  result = ${{lastEvent.result}}\n  final cands = ${{lastEvent.final_cands}}\n  fragment final size = ${{lastEvent.fragment.length}}\n  fragment atoms = [${{lastEvent.fragment.join(', ')}}]`;
    if (lastEvent.iso) txt += `\n  iso = ${{JSON.stringify(lastEvent.iso)}}`;
  }} else if (lastEvent.type === 'island_locked') {{
    txt += `  island #${{lastEvent.island_idx}} locked\n`;
    txt += `  new pairs (${{lastEvent.pairs.length}}): ${{lastEvent.pairs.map(p => `R[${{p[0]}}]→P[${{p[1]}}]`).join(', ')}}\n`;
    txt += `  total mapped = ${{lastEvent.mapped_total}}`;
    if (lastEvent.merged_with && lastEvent.merged_with.length) {{
      txt += `\n  merged with islands: [${{lastEvent.merged_with.join(', ')}}]`;
      txt += `\n  relabeled atoms: ${{lastEvent.relabeled ? lastEvent.relabeled.map(([r,old]) => `R[${{r}}](was #${{old}})`).join(', ') : ''}}`;
    }}
    if (lastEvent.parent_atom !== undefined) {{
      txt += `\n  expand-merge with parent R[${{lastEvent.parent_atom}}] (pass ${{lastEvent.expand_pass}})`;
    }}
  }} else if (lastEvent.type === 'pass_start') {{
    txt += `  pass ${{lastEvent.pass}},  mapped so far = ${{lastEvent.mapped}}`;
  }} else if (lastEvent.type === 'done') {{
    txt += `  final mapped = ${{lastEvent.mapped}}`;
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


def main():
    step = sys.argv[1] if len(sys.argv) > 1 else 'pr13.Cyclobutane_JOC2023_TS-CD_step1'
    bench = Path('/Users/yunhengz/empty_for_claude/Benchmark') / step / 'plain' / 'stage0'
    work = Path('/Users/yunhengz/empty_for_claude/rxn_core/work_frag') / step
    out = Path('/Users/yunhengz/empty_for_claude/rxn_core/out') / f'animate_{step}.html'

    elR, xyzR_arr, wboR = run_xtb(bench / 'reactant.xyz', work / 'R')
    elP, xyzP_arr, wboP = run_xtb(bench / 'product.xyz', work / 'P')
    g_R = build_graph(elR, wboR)
    g_P = build_graph(elP, wboP)
    mapping, events = find_islands_traced(g_R, g_P, wbo_tol=0.5)

    print(f'{step}: {len(events)} events, {len(mapping)} atoms mapped')

    xyzR_str = write_xyz_str(elR, xyzR_arr, comment='reactant')
    xyzP_str = write_xyz_str(elP, xyzP_arr, comment='product')
    html = HTML.format(
        title=step,
        xyzR_json=json.dumps(xyzR_str),
        xyzP_json=json.dumps(xyzP_str),
        events_json=json.dumps(events),
        wboR_json=json.dumps(wboR.tolist()),
        wboP_json=json.dumps(wboP.tolist()),
        elements_R_json=json.dumps(elR),
        elements_P_json=json.dumps(elP),
    )
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(html)
    print(f'wrote {out}')


if __name__ == '__main__':
    main()
