"""Build interactive multi-mechanism views for all BGCP steps.

For each step:
  1. R<->P cut-sweep -> min-bondcount mechanisms.
  2. Under each mechanism: core-match GT + each IG, pick best R<->T core witness.
  3. Rank IGs per mech, mark top-2; union across mechs.
  4. Write out/bgcp_views/<step>/view.html with mechanism switcher.
  5. Dump out/bgcp_alignment_eval_v2.json for downstream CSV.

Parallelized via multiprocessing.
"""
from __future__ import annotations
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent / "src"))

import argparse
import json
import multiprocessing as mp
import os
import re
import time
import traceback
from pathlib import Path
import numpy as np

from rxn_core import (parse_xyz, classify_bonds, parse_g98_modes,
                      core_atoms_in_R_frame,
                      reaction_coord_delta, reindex_modes_to_R,
                      bond_overlap_per_mode, bond_reaction_vector,
                      rxn_overlap_per_mode,
                      build_graph, cut_sweep, select_min_mechanisms,
                      ts_core_pool)
from rxn_core.matcher import _nauty_orbits

PROJECT = _Path(__file__).resolve().parent
WORK = PROJECT / "appendix_perparation" / "xtb_frequency_calculations"
OUT_ROOT = PROJECT / "out" / "bgcp_views"
EVAL_JSON = PROJECT / "out" / "bgcp_alignment_eval_v2.json"
CUT_FLOOR = float(os.environ.get("BGCP_CUT_FLOOR", "0.2"))
N_SEEDS_PER_RUN = 3  # cut + seed are orthogonal diversity sources; keep both modest
VIEW_MAX_BRANCHES = int(os.environ.get("BGCP_VIEW_MAX_BRANCHES", "5000"))
CUTSWEEP_CHUNKSIZE = int(os.environ.get("BGCP_CUTSWEEP_CHUNKSIZE", "1"))
VIEW_ISO_TOL = float(os.environ.get("BGCP_ISO_TOL", "0.5"))
UNIT_TIMEOUT = float(os.environ.get("BGCP_UNIT_TIMEOUT", "10"))
BGCP_TIMING = os.environ.get("BGCP_TIMING", "0") == "1"
SYMMETRY_REPAIR = os.environ.get("BGCP_SYMMETRY_REPAIR", "1") != "0"
SYMMETRY_REPAIR_MIN_CHANGES = int(os.environ.get("BGCP_SYMMETRY_REPAIR_MIN_CHANGES", "5"))
SYMMETRY_REPAIR_MAX_EVALS = int(os.environ.get("BGCP_SYMMETRY_REPAIR_MAX_EVALS", "20000"))
TS_CORE_EDGE_FLOOR = float(os.environ.get("BGCP_TS_CORE_EDGE_FLOOR", "0.2"))
TS_CORE_MAX_CANDIDATES = int(os.environ.get("BGCP_TS_CORE_MAX_CANDIDATES", "20000"))
W_RXN, W_CORE, IMAG_PEN = 1.0, 0.2, 0.3


def load(d):
    xyz_path = next(p for p in d.glob("*.xyz") if "xtbhess" not in p.name)
    el, xyz = parse_xyz(xyz_path)
    n = len(el); wbo = np.zeros((n, n))
    for ln in (d / "wbo").read_text().splitlines():
        p = ln.split()
        if len(p) < 3: continue
        i, j = int(p[0])-1, int(p[1])-1
        wbo[i, j] = float(p[2]); wbo[j, i] = wbo[i, j]
    return el, np.asarray(xyz, float), wbo


def _bond_key(bonds, orbits=None):
    pairs = []
    for a, b in bonds:
        a = int(a); b = int(b)
        if orbits is not None:
            a = int(orbits[a]); b = int(orbits[b])
        if a > b:
            a, b = b, a
        pairs.append((a, b))
    return tuple(sorted(pairs))


def _mechanism_bond_key(mech, r_orbits):
    return (
        _bond_key(mech['broken_bonds_R'], r_orbits),
        _bond_key(mech['formed_bonds_R'], r_orbits),
    )


def _gt_score(mech):
    gt = mech.get('gt')
    return float(gt['S']) if gt and gt.get('S') is not None else float('-inf')


def dedupe_mechanisms_by_bond_changes(mechanisms, r_orbits):
    """Collapse final-view mechanisms with the same R-symmetry bond changes.

    The cut sweep may find multiple concrete alignments whose broken/formed
    R-index bonds differ only by swapping equivalent reactant atoms.  They are
    the same mechanism for the view, so keep the highest-GT-scoring
    representative and retain provenance for the slim JSON / button tooltip.
    """
    groups = {}
    for mech in mechanisms:
        key = _mechanism_bond_key(mech, r_orbits)
        groups.setdefault(key, []).append(mech)

    deduped = []
    for group in groups.values():
        rep = max(group, key=_gt_score)
        rep['dedup_count'] = sum(m.get('dedup_count', 1) for m in group)
        rep['dedup_source_ids'] = [int(m['id']) for m in group]
        rep['dedup_cuts'] = sorted({
            cut
            for m in group
            for cut in m.get('dedup_cuts', [m['cut']])
        })
        deduped.append(rep)

    for new_id, mech in enumerate(deduped, 1):
        suffix = re.sub(r"^#\d+:\s*", "", mech['label'])
        if mech['dedup_count'] > 1:
            suffix = f"{suffix} [dedup x{mech['dedup_count']}]"
        mech['id'] = new_id
        mech['label'] = f"#{new_id}: {suffix}"
    return deduped


def score_one(elR, xyzR, elT, xyzT, mapping_RT, freqs, modes_TS,
              broken_R, formed_R, core_R, delta_RP):
    modes_R = reindex_modes_to_R(modes_TS, mapping_RT, len(elR))
    mode_norms = np.linalg.norm(modes_TS.reshape(modes_TS.shape[0], -1), axis=1)
    sq = (modes_R**2).sum(axis=2)
    total = mode_norms ** 2
    core_e = sq[:, core_R].sum(axis=1) if core_R else np.zeros(modes_R.shape[0])
    kappa = np.where(total > 1e-12, core_e / total, 0.0)
    rho = rxn_overlap_per_mode(modes_R, delta_RP, core_R,
                                mode_norms=mode_norms)
    ts_in_R = np.asarray(xyzR, float).copy()
    for r, t in mapping_RT.items(): ts_in_R[r] = xyzT[t]
    V = bond_reaction_vector(ts_in_R, broken_R, formed_R)
    beta = bond_overlap_per_mode(modes_R, V, mode_norms=mode_norms)
    imag = list(np.where(freqs < 0)[0])
    if not imag: return None
    pk = max(imag, key=lambda k: beta[k])
    return {'S': float(beta[pk]*(1+W_RXN*rho[pk])*(1+W_CORE*kappa[pk])/max(len(imag),1)**IMAG_PEN),
            'beta': float(beta[pk]), 'rho': float(rho[pk]), 'kappa': float(kappa[pk]),
            'freq': float(freqs[pk]), 'k': int(pk), 'n_imag': len(imag),
            'core_map': {str(int(r)): int(mapping_RT[r])
                         for r in core_R if r in mapping_RT},
            'picked_disp': modes_R[pk].tolist(),
            'xyz_in_R': ts_in_R.tolist()}


def best_under_mech_using_pool(elR, xyzR, elT, xyzT, freqs, modes_TS,
                                 rt_pool, broken_R, formed_R, core_R, delta_RP):
    """Score every R<->T witness under one mech, with CORE-RESTRICTED DEDUP.

    Two witnesses that agree on (r -> mapping[r] for r in core_R) produce
    identical verifier scores under this mechanism: beta, rho, kappa, and
    the bond-reaction vector V only read TS coords / mode displacements at
    R-frame *core* indices. Spectator-atom permutations on the TS side
    don't affect those numerators; the denominator uses the full TS mode norm
    directly, so no spectator mapping is invented.

    We dedup the pool by per-mechanism core-restricted key and score one
    rep per equivalence class. Keep highest-S rep. This is option (a)
    from the design discussion — per-mechanism, not unioned across mechs,
    because different mechs have different core_R sets.
    """
    core_R_set = frozenset(core_R)
    seen_core = set()
    best = None
    for v in rt_pool.values():
        witness = v['mapping']
        # Per-mechanism core-restricted key
        core_key = frozenset((r, witness[r]) for r in core_R_set if r in witness)
        if core_key in seen_core:
            continue
        seen_core.add(core_key)
        mapping_RT = dict(witness)
        s = score_one(elR, xyzR, elT, xyzT, mapping_RT, freqs, modes_TS,
                      broken_R, formed_R, core_R, delta_RP)
        if s and (best is None or s['S'] > best['S']):
            best = s
    return best


HTML = r"""<!doctype html><html><head><meta charset="utf-8"><title>{title}</title>
<script src="https://3dmol.org/build/3Dmol-min.js"></script>
<style>
html,body{{margin:0;padding:0;font-family:-apple-system,sans-serif;background:#fafafa}}
body{{padding:14px;box-sizing:border-box}}
h2{{margin:0 0 4px;font-size:18px}}
.mech-sel{{background:white;padding:10px;border:1px solid #ccc;border-radius:6px;margin-bottom:14px}}
.mech-sel button{{padding:6px 12px;margin-right:6px;border:1px solid #aaa;background:#f0f0f0;border-radius:4px;cursor:pointer;font-family:ui-monospace,monospace;font-size:12px}}
.mech-sel button.active{{background:#ffd700;border-color:#a90;font-weight:600}}
.ref-row{{display:grid;grid-template-columns:repeat(3,1fr);gap:12px;margin-bottom:18px}}
.ig-grid{{display:grid;grid-template-columns:repeat(4,1fr);gap:12px}}
.panel{{background:white;border:1px solid #ddd;border-radius:6px;padding:6px 8px 8px}}
.panel.top2{{border:2px solid #d4af37}}
.panel.union{{box-shadow:0 0 0 2px #ff9}}
.ph{{display:flex;justify-content:space-between;align-items:baseline;font-size:12px;margin-bottom:4px}}
.ph .lbl{{font-weight:600;font-size:13px}}
.ph .rk{{font-family:ui-monospace,monospace;color:#024}}
.vw{{position:relative;width:100%;height:230px}}
.ref-row .vw{{height:300px}}
.vwbox{{position:absolute;inset:0}}
.meta{{font-family:ui-monospace,monospace;font-size:11px;color:#444;padding:3px 0 0;line-height:1.4}}
.meta b{{color:#024}}
</style></head><body>
<h2>{title}</h2>
<div class="mech-sel" id="mech-sel"></div>
<div class="ref-row">
  <div class="panel"><div class="ph"><span class="lbl">Reactant</span><span class="rk">static</span></div>
    <div class="vw"><div id="vw_R" class="vwbox"></div></div></div>
  <div class="panel"><div class="ph"><span class="lbl">Product</span><span class="rk" id="prod_label">static</span></div>
    <div class="vw"><div id="vw_P" class="vwbox"></div></div></div>
  <div class="panel"><div class="ph"><span class="lbl">Ground-truth TS</span><span class="rk" id="gt_S">S=?</span></div>
    <div class="vw"><div id="vw_GT" class="vwbox"></div></div>
    <div class="meta" id="gt_meta"></div></div>
</div>
<div class="ig-grid" id="grid"></div>
<script>
const DATA = {data_json};
let currentMechId = DATA.default_mech_id;
const elements = DATA.reactant.elements;
const xyzR_static = DATA.reactant.coords;
function findMech(id) {{ return DATA.mechanisms.find(m=>m.id===id); }}
function buildBody(els, xyz) {{ let s = els.length+"\nframe\n"; for (let i=0;i<xyz.length;i++) s += els[i]+"  "+xyz[i][0].toFixed(6)+"  "+xyz[i][1].toFixed(6)+"  "+xyz[i][2].toFixed(6)+"\n"; return s; }}
function buildBodyAt(els, xyz, disp, scale) {{ let s = els.length+"\nframe\n"; for (let i=0;i<xyz.length;i++) {{ const x=xyz[i][0]+scale*disp[i][0], y=xyz[i][1]+scale*disp[i][1], z=xyz[i][2]+scale*disp[i][2]; s += els[i]+"  "+x.toFixed(6)+"  "+y.toFixed(6)+"  "+z.toFixed(6)+"\n"; }} return s; }}
function xyzAt(xyz, disp, scale) {{ return xyz.map((p,i)=>[p[0]+scale*disp[i][0], p[1]+scale*disp[i][1], p[2]+scale*disp[i][2]]); }}
const animTimers = {{}};
function stopAnim(d) {{ if (animTimers[d]) {{ clearInterval(animTimers[d]); delete animTimers[d]; }} }}
function drawBonds(v, xyz, pairs, color) {{ for (const [i,j] of pairs) {{ if (i>=xyz.length||j>=xyz.length) continue; v.addCylinder({{start:{{x:xyz[i][0],y:xyz[i][1],z:xyz[i][2]}}, end:{{x:xyz[j][0],y:xyz[j][1],z:xyz[j][2]}}, color:color, radius:0.16, dashed:true}}); }} }}
function drawArrows(v, xyz, disp, core) {{ for (const i of core) {{ if (!disp||!disp[i]) continue; const d = disp[i]; const len = Math.hypot(d[0],d[1],d[2]); if (len<0.05) continue; v.addArrow({{start:{{x:xyz[i][0],y:xyz[i][1],z:xyz[i][2]}}, end:{{x:xyz[i][0]+d[0]*1.5,y:xyz[i][1]+d[1]*1.5,z:xyz[i][2]+d[2]*1.5}}, color:'#0066cc', radius:0.07}}); }} }}
function makeStatic(divId, els, xyz, broken, formed) {{ stopAnim(divId); document.getElementById(divId).innerHTML=""; const v = $3Dmol.createViewer(divId, {{backgroundColor:'white'}}); v.addModel(buildBody(els, xyz), 'xyz'); v.setStyle({{}}, {{stick:{{radius:0.10}}, sphere:{{scale:0.20}}}}); drawBonds(v, xyz, broken, 'red'); drawBonds(v, xyz, formed, 'green'); v.zoomTo(); v.render(); return v; }}
function makeAnimated(divId, els, xyz, disp, broken, formed, core) {{ stopAnim(divId); document.getElementById(divId).innerHTML=""; const v = $3Dmol.createViewer(divId, {{backgroundColor:'white'}}); v.addModel(buildBody(els, xyz), 'xyz'); v.setStyle({{}}, {{stick:{{radius:0.10}}, sphere:{{scale:0.20}}}}); drawBonds(v, xyz, broken, 'red'); drawBonds(v, xyz, formed, 'green'); drawArrows(v, xyz, disp, core); v.zoomTo(); v.render(); let t=0; const period=30, amp=0.6; animTimers[divId] = setInterval(()=>{{ t=(t+1)%period; const scale = amp*Math.sin(2*Math.PI*t/period); const cur = xyzAt(xyz, disp, scale); v.removeAllModels(); v.removeAllShapes(); v.addModel(buildBodyAt(els, xyz, disp, scale), 'xyz'); v.setStyle({{}}, {{stick:{{radius:0.10}}, sphere:{{scale:0.20}}}}); drawBonds(v, cur, broken, 'red'); drawBonds(v, cur, formed, 'green'); drawArrows(v, cur, disp, core); v.render(); }}, 60); return v; }}
function render() {{ const mech = findMech(currentMechId); document.querySelectorAll('.mech-sel button').forEach(b => {{ b.classList.toggle('active', parseInt(b.dataset.id)===currentMechId); }}); makeStatic('vw_R', elements, xyzR_static, mech.broken_bonds_R, []); makeStatic('vw_P', elements, mech.product_xyz_in_R, [], mech.formed_bonds_R); document.getElementById('prod_label').textContent = "static (mech #"+mech.id+")"; if (mech.gt && mech.gt.picked_disp) {{ makeAnimated('vw_GT', elements, mech.gt.xyz_in_R, mech.gt.picked_disp, mech.broken_bonds_R, mech.formed_bonds_R, mech.core_atoms); document.getElementById('gt_S').textContent = "S = "+mech.gt.S.toFixed(3); document.getElementById('gt_meta').innerHTML = "<b>&beta;</b>="+mech.gt.beta.toFixed(3)+" &nbsp; <b>&rho;</b>="+mech.gt.rho.toFixed(3)+" &nbsp; <b>&kappa;</b>="+mech.gt.kappa.toFixed(3)+" &nbsp; <b>n_imag</b>="+mech.gt.n_imag+" &nbsp; <b>freq</b>="+mech.gt.freq.toFixed(0)+"i cm&#x207B;&#xB9;"; }} const grid = document.getElementById('grid'); grid.innerHTML = ""; const igs = [...mech.igs].sort((a,b) => (b.S||0) - (a.S||0)); igs.forEach((ig, idx) => {{ const div = document.createElement('div'); let cls = 'panel'; if (ig.is_top2) cls += ' top2'; if (ig.is_union_top && !ig.is_top2) cls += ' union'; div.className = cls; const sStr = ig.S !== undefined ? "S = "+ig.S.toFixed(3) : "no score"; const tag = ig.is_top2 ? '<span style="background:#d4af37;color:white;padding:1px 6px;border-radius:3px;font-size:11px;margin-left:4px">TOP2</span>' : (ig.is_union_top ? '<span style="background:#ff9;color:#660;padding:1px 6px;border-radius:3px;font-size:11px;margin-left:4px">union</span>' : ''); div.innerHTML = '<div class="ph"><span class="lbl">'+ig.label+tag+'</span><span class="rk">'+sStr+'</span></div><div class="vw"><div id="vw_ig'+idx+'" class="vwbox"></div></div><div class="meta">'+(ig.beta!==undefined ? "<b>&beta;</b>="+ig.beta.toFixed(3)+" <b>&rho;</b>="+ig.rho.toFixed(3)+" <b>&kappa;</b>="+ig.kappa.toFixed(3)+" <b>n_imag</b>="+ig.n_imag+" <b>freq</b>="+ig.freq.toFixed(0)+"i" : "(no data)")+"</div>"; grid.appendChild(div); if (ig.picked_disp) makeAnimated("vw_ig"+idx, elements, ig.xyz_in_R, ig.picked_disp, mech.broken_bonds_R, mech.formed_bonds_R, mech.core_atoms); else if (ig.xyz_in_R) makeStatic("vw_ig"+idx, elements, ig.xyz_in_R, mech.broken_bonds_R, mech.formed_bonds_R); }}); }}
const ms = document.getElementById('mech-sel'); ms.innerHTML = "<span style='font-size:13px;margin-right:8px;color:#444'>Mechanism:</span>"; DATA.mechanisms.forEach(m => {{ const b = document.createElement('button'); b.dataset.id = m.id; b.textContent = m.label + "  GT S=" + (m.gt ? m.gt.S.toFixed(3) : '?'); if ((m.dedup_count||1) > 1) b.title = "Collapsed raw witnesses: "+m.dedup_count+"; source mechanisms: "+m.dedup_source_ids.join(", ")+"; cuts: "+m.dedup_cuts.join(", "); b.onclick = () => {{ currentMechId = m.id; render(); }}; ms.appendChild(b); }});
window.addEventListener('load', render);
</script>
</body></html>
"""


def process_step(step_name, inner_workers=0):
    """inner_workers: parallelism inside this step's cut_sweeps.
      0 / 1  → serial cut_sweep (safe when called inside outer mp.Pool worker)
      >= 2   → parallel cut_sweep on that many workers (use only when there
               is no outer Pool, i.e. single-step CLI mode)."""
    try:
        sd = WORK / step_name
        if not (sd / "R" / "wbo").exists() or not (sd / "P" / "wbo").exists():
            return {"step": step_name, "error": "missing R or P xtb cache"}
        if not (sd / "sp_groundtruth").exists() or not (sd / "hess_groundtruth" / "g98.out").exists():
            return {"step": step_name, "error": "missing GT"}
        elR, xyzR, wboR = load(sd / "R")
        elP, xyzP, wboP = load(sd / "P")
        elT_gt, xyzT_gt, wboT_gt = load(sd / "sp_groundtruth")
        freqs_gt, modes_gt = parse_g98_modes(sd / "hess_groundtruth" / "g98.out")

        def timed_cut_sweep(label, elT, wboT, core_R=None):
            t_sweep = time.time()
            pool = cut_sweep(elR, wboR, elT, wboT,
                             n_workers=inner_workers, core_R=core_R,
                             cut_floor=CUT_FLOOR,
                             graph_floor=0.2,
                             iso_tol=VIEW_ISO_TOL,
                             n_seeds=N_SEEDS_PER_RUN,
                             max_branches=VIEW_MAX_BRANCHES,
                             chunksize=CUTSWEEP_CHUNKSIZE,
                             unit_timeout=UNIT_TIMEOUT,
                             symmetry_repair=SYMMETRY_REPAIR,
                             symmetry_repair_min_changes=SYMMETRY_REPAIR_MIN_CHANGES,
                             symmetry_repair_max_evals=SYMMETRY_REPAIR_MAX_EVALS)
            if BGCP_TIMING:
                core_msg = f" core={len(core_R)}" if core_R else ""
                print(f"    {step_name} {label:>12s} cut_sweep: "
                      f"{len(pool):>4d} sigs{core_msg} "
                      f"in {time.time()-t_sweep:.1f}s",
                      flush=True)
            return pool

        def timed_ts_core_pool(label, elT, wboT, core_R,
                               broken_R=None, formed_R=None):
            t_match = time.time()
            pool = ts_core_pool(elR, wboR, elT, wboT, core_R,
                                broken_R=broken_R, formed_R=formed_R,
                                edge_floor=TS_CORE_EDGE_FLOOR,
                                iso_tol=VIEW_ISO_TOL,
                                max_candidates=TS_CORE_MAX_CANDIDATES)
            if BGCP_TIMING:
                if len(pool) >= TS_CORE_MAX_CANDIDATES:
                    print(f"    [warn] TS core pool hit cap={TS_CORE_MAX_CANDIDATES} "
                          f"core={list(core_R)}",
                          flush=True)
                print(f"    {step_name} {label:>12s} core_match: "
                      f"{len(pool):>4d} sigs core={len(core_R)} "
                      f"in {time.time()-t_match:.1f}s",
                      flush=True)
            return pool

        rp = timed_cut_sweep("R-P", elP, wboP)
        rp_min = select_min_mechanisms(rp)
        if not rp_min:
            return {"step": step_name, "error": "no min-bond mechanism"}

        mechanisms = []
        for mi, (_sig, info) in enumerate(rp_min.items(), 1):
            mapping_RP = info['mapping']
            inv_RP = {v: k for k, v in mapping_RP.items()}
            broken, formed, _, _ = classify_bonds(mapping_RP, wboR, wboP)
            broken_R = [(int(a), int(b)) for (a, b, _, _) in broken]
            formed_R = [(int(inv_RP[a]), int(inv_RP[b])) for (a, b, _, _) in formed
                        if a in inv_RP and b in inv_RP]
            core_R = list(core_atoms_in_R_frame(mapping_RP, broken, formed))
            delta_RP = reaction_coord_delta(np.asarray(xyzR), np.asarray(xyzP), mapping_RP)
            xyzP_in_R = np.asarray(xyzR, float).copy()
            for i_R, i_P in mapping_RP.items(): xyzP_in_R[i_R] = xyzP[i_P]
            cut = next(iter(info['cuts']), None)
            cut_name = f"{elR[cut[0]]}{cut[0]}-{elR[cut[1]]}{cut[1]}" if cut else "none"
            br_label = ",".join(f"{elR[a]}{a}-{elR[b]}{b}" for a, b in broken_R)
            mech = {
                'id': mi, 'cut': cut_name,
                'label': f"#{mi}: {br_label} (cut: {cut_name})",
                'dedup_count': info.get('dedup_count', 1),
                'dedup_cuts': [
                    f"{elR[a]}{a}-{elR[b]}{b}" for a, b in sorted(info['cuts'])
                ] or [cut_name],
                'broken_bonds_R': broken_R, 'formed_bonds_R': formed_R,
                'core_atoms': core_R,
                'product_xyz_in_R': xyzP_in_R.tolist(),
            }
            mech['_state'] = (broken_R, formed_R, core_R, delta_RP)
            mechanisms.append(mech)

        r_orbits = _nauty_orbits(build_graph(elR, wboR, bond_cut=0.2),
                                 wbo_tol=0.2)
        mechanisms = dedupe_mechanisms_by_bond_changes(mechanisms, r_orbits)
        # R<->T core matching is mechanism-local.  A target-side symmetry
        # representative that is arbitrary for one core can become unique once
        # a particular mechanism's broken/formed bonds define the scoring core.
        # Spectator alternatives are not enumerated.
        for mech in mechanisms:
            br_R, fm_R, core_R, dRP = mech['_state']
            gt_rt_pool = timed_ts_core_pool(
                f"GT:m{mech['id']}", elT_gt, wboT_gt, core_R,
                broken_R=br_R, formed_R=fm_R)
            mech['gt'] = best_under_mech_using_pool(
                elR, xyzR, elT_gt, xyzT_gt, freqs_gt, modes_gt, gt_rt_pool,
                br_R, fm_R, core_R, dRP)

        # IGs: enumerate mechanism-local core alternatives, then score under
        # each mechanism. Spectators are not greedily materialized; unmapped
        # viewer coordinates remain at the R-frame placeholder.
        iter_dirs = sorted([d for d in sd.iterdir()
                            if d.is_dir() and re.match(r"hess_iter(\d+)$", d.name)],
                           key=lambda d: int(re.match(r"hess_iter(\d+)$", d.name).group(1)))
        for mech in mechanisms: mech['igs'] = []
        for hess_dir in iter_dirs:
            label = hess_dir.name.replace("hess_", "")
            sp_dir = sd / f"sp_{label}"
            try:
                elI, xyzI, wboI = load(sp_dir)
                freqs_i, modes_i = parse_g98_modes(hess_dir / "g98.out")
            except Exception:
                for mech in mechanisms: mech['igs'].append({'label': label})
                continue
            for mech in mechanisms:
                br_R, fm_R, core_R, dRP = mech['_state']
                ig_rt_pool = timed_ts_core_pool(
                    f"{label}:m{mech['id']}", elI, wboI, core_R,
                    broken_R=br_R, formed_R=fm_R)
                s = best_under_mech_using_pool(elR, xyzR, elI, xyzI,
                                                 freqs_i, modes_i, ig_rt_pool,
                                                 br_R, fm_R, core_R, dRP)
                entry = {'label': label}
                if s: entry.update(s)
                mech['igs'].append(entry)

        union_top = set()
        for mech in mechanisms:
            ranked = sorted([(i, ig) for i, ig in enumerate(mech['igs']) if ig.get('S') is not None],
                            key=lambda x: -x[1]['S'])
            top2 = {i for i, _ in ranked[:2]}
            for i, ig in enumerate(mech['igs']):
                ig['is_top2'] = (i in top2)
                if i in top2: union_top.add(ig['label'])
        for mech in mechanisms:
            for ig in mech['igs']:
                ig['is_union_top'] = ig['label'] in union_top
            del mech['_state']

        default_id = max(mechanisms, key=lambda m: m['gt']['S'] if m['gt'] else 0)['id']
        data = {'step': step_name, 'n_atoms': len(elR),
                'reactant': {'elements': elR, 'coords': np.asarray(xyzR).tolist()},
                'mechanisms': mechanisms, 'default_mech_id': default_id}

        run_dir = OUT_ROOT / step_name
        run_dir.mkdir(parents=True, exist_ok=True)
        (run_dir / "view.html").write_text(HTML.format(
            title=f"BGCP &mdash; {step_name}  ({len(mechanisms)} mechanisms)",
            data_json=json.dumps(data),
        ))
        # stripped down for the eval JSON (no big xyz/disp arrays)
        slim = {'step': step_name, 'n_atoms': len(elR), 'n_mechs': len(mechanisms),
                'mechanisms': []}
        for mech in mechanisms:
            slim['mechanisms'].append({
                'id': mech['id'], 'cut': mech['cut'],
                'dedup_count': mech.get('dedup_count', 1),
                'dedup_source_ids': mech.get('dedup_source_ids', [mech['id']]),
                'dedup_cuts': mech.get('dedup_cuts', [mech['cut']]),
                'broken_R': mech['broken_bonds_R'], 'formed_R': mech['formed_bonds_R'],
                'core_R': mech['core_atoms'],
                'gt': {k: mech['gt'][k] for k in ['S', 'beta', 'rho', 'kappa', 'freq', 'n_imag', 'core_map']} if mech['gt'] else None,
                'igs': [{k: ig.get(k) for k in ['label', 'S', 'beta', 'rho', 'kappa', 'freq', 'n_imag', 'core_map', 'is_top2']}
                        for ig in mech['igs']],
            })
        # Per-step slim record (so parallel Slurm array tasks don't race on
        # the global EVAL_JSON). The post-run merge step reads these and
        # builds the global EVAL_JSON. See nrt_verification_workflow.sh.
        (run_dir / "_eval_v2_slim.json").write_text(json.dumps(slim))
        return {'step': step_name, 'slim': slim,
                'top1_label': max(mechanisms, key=lambda m: m['gt']['S'] if m['gt'] else 0)['igs'][0]['label'] if mechanisms[0]['igs'] else "?"}
    except Exception as e:
        return {"step": step_name, "error": f"{type(e).__name__}: {e}",
                "trace": traceback.format_exc()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=max(1, os.cpu_count() - 1),
                    help="Outer parallelism: how many steps run concurrently. "
                         "Ignored when --inner-workers > 1.")
    ap.add_argument("--inner-workers", type=int, default=0,
                    help="Inner parallelism: how many workers each step's "
                         "cut_sweep uses. Default 0 = serial inside step. "
                         "Setting > 1 disables --workers (outer) to avoid "
                         "nested daemonic multiprocessing.")
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--steps", nargs="+", default=None)
    args = ap.parse_args()

    all_steps = sorted(d.name for d in WORK.iterdir() if d.is_dir())
    if args.steps: steps = [s for s in all_steps if s in set(args.steps)]
    elif args.limit: steps = all_steps[:args.limit]
    else: steps = all_steps

    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    eval_records = []
    n_ok = n_err = 0

    def _record(i, rec):
        nonlocal n_ok, n_err
        if rec.get("error"):
            print(f"  [{i:>3d}/{len(steps)}] {rec['step']:60s}  ERROR: {rec['error'][:80]}", flush=True)
            n_err += 1
            eval_records.append({'step': rec['step'], 'error': rec['error']})
        else:
            slim = rec['slim']
            gt_best = max((m['gt']['S'] for m in slim['mechanisms'] if m['gt']), default=0)
            print(f"  [{i:>3d}/{len(steps)}] {rec['step']:60s}  mechs={slim['n_mechs']}  best_GT_S={gt_best:.3f}", flush=True)
            eval_records.append(slim)
            n_ok += 1

    if args.inner_workers and args.inner_workers > 1:
        # Inner-parallel mode: steps run serially in main; each step's
        # cut_sweep uses inner_workers cores. Best for a single step or a
        # few large steps where the cut_sweep itself dominates cost.
        print(f"Processing {len(steps)} steps serially; each step uses "
              f"{args.inner_workers} inner workers "
              f"(cut_sweep chunksize={CUTSWEEP_CHUNKSIZE}, "
              f"iso_tol={VIEW_ISO_TOL}, "
              f"unit_timeout={UNIT_TIMEOUT}s)")
        for i, step in enumerate(steps, 1):
            rec = process_step(step, inner_workers=args.inner_workers)
            _record(i, rec)
    else:
        # Outer-parallel mode: args.workers steps run concurrently; each
        # step's cut_sweep is serial (no nested daemonic Pool). Best for
        # the full 155-step benchmark.
        print(f"Processing {len(steps)} steps with {args.workers} outer workers")
        with mp.Pool(args.workers) as pool:
            for i, rec in enumerate(pool.imap_unordered(process_step, steps), 1):
                _record(i, rec)

    print(f"\n{n_ok} ok, {n_err} errors in {time.time()-t0:.0f}s")

    EVAL_JSON.write_text(json.dumps(eval_records))
    print(f"wrote {EVAL_JSON}  ({n_ok} step records)")


if __name__ == "__main__":
    main()
