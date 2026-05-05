"""
Compute Hessian + mode ranking for every (step, TS) in
Benchmark_Guesses_Collective_Package, cache to disk, and emit:

  out_initial_posing_mode/
    index.html              — step picker linking to per-step pages
    <step_id>.html          — mode viewer for one step (all 21 TSes)
    <step_id>.json          — raw cached data per step (for re-use)

xtb hess + sp results are cached in work_modes/ via the global
run_xtb / run_xtb_hess caching layer (xyz-content match). Per-step
JSON caches the alignment + ranking so a future change of alignment
or ranking algorithm only needs to re-run the cheap parts.

JSON schema per step:
  {
    "step_id":      str,
    "charge":       int,
    "uhf":          int,
    "n_atoms_R":    int,
    "core_R":       [int, ...],         # core atoms in R-frame
    "broken":       [[i, j, wR, wP], ...],
    "formed":       [[ip, jp, wR, wP], ...],
    "ts_list": [
      {
        "label":            "groundtruth" | "iter1" | ...,
        "ts_file":          str,
        "elements":         [str, ...],
        "coords":           [[x, y, z], ...],          # n_atoms x 3
        "core_atoms_ts":    [int, ...],                # core in TS-frame
        "map_R_to_TS":      {r_idx: ts_idx, ...},
        "n_imag":           int,
        "freqs":            [float, ...],              # all freqs
        "modes_ranked":     [
          {
            "mode_idx":      int,
            "freq":          float,
            "core_fraction": float,
            "displacement":  [[dx, dy, dz], ...],      # n_atoms x 3
          },
          ...
        ],
      },
      ...
    ],
  }
"""
from __future__ import annotations
import argparse
import json
import re
import sys
import time
import traceback
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from rxn_core_frag import (
    run_xtb, run_xtb_hess, build_graph, find_islands, expand_mapping,
    classify_bonds, _generate_seed_orders,
)
from build_bgcp_viewer import (
    BGCP_ROOT, LOOKUP, list_initial_guesses, iter_num, read_xyzs,
)


ROOT = Path("/Users/yunhengz/empty_for_claude/rxn_core")
OUT_DIR = ROOT / "out_initial_posing_mode"
WORK = ROOT / "work_modes"
OUT_DIR.mkdir(parents=True, exist_ok=True)
WORK.mkdir(parents=True, exist_ok=True)


def best_mapping(g_a, g_b, wbo_a, wbo_b, n_seeds=10):
    orders = _generate_seed_orders(g_a, n_seeds)
    best = None
    for order in orders:
        m, _ = find_islands(g_a, g_b, seed_order=order)
        m = expand_mapping(m, g_a, g_b)
        br, fm, _, _ = classify_bonds(m, wbo_a, wbo_b)
        score = (len(br) + len(fm), -len(m))
        if best is None or score < best[0]:
            best = (score, m)
    return best[1]


def core_atoms_R_frame(map_R_to_P, broken, formed):
    inv = {v: k for k, v in map_R_to_P.items()}
    core = set()
    for (i, j, _, _) in broken:
        core.add(i); core.add(j)
    for (ip, jp, _, _) in formed:
        if ip in inv: core.add(inv[ip])
        if jp in inv: core.add(inv[jp])
    return sorted(core)


def core_in_TS_frame(core_R, map_R_to_TS):
    return sorted(map_R_to_TS[r] for r in core_R if r in map_R_to_TS)


def rank_modes(modes, freqs, core_indices_ts):
    """Imag first (by core_fraction desc), then real (by core_fraction desc)."""
    if modes.shape[0] == 0:
        return []
    sq = (modes ** 2).sum(axis=2)
    total = sq.sum(axis=1)
    if not core_indices_ts:
        # No core atoms — fall back to ranking by freq magnitude
        return [{
            "mode_idx": int(i),
            "freq": float(freqs[i]),
            "core_fraction": 0.0,
            "displacement": modes[i].astype(float).tolist(),
        } for i in range(len(freqs))]
    core = np.array(core_indices_ts, dtype=int)
    core_e = sq[:, core].sum(axis=1)
    fraction = np.where(total > 1e-9, core_e / total, 0.0)
    out = []
    for i in range(len(freqs)):
        out.append({
            "mode_idx": int(i),
            "freq": float(freqs[i]),
            "core_fraction": float(fraction[i]),
            "displacement": modes[i].astype(float).tolist(),
        })
    out.sort(key=lambda d: (0 if d["freq"] < 0 else 1, -d["core_fraction"]))
    return out


def process_ts(label, ts_path, wd, chg, uhf, base):
    ts_local = wd / f"{label}__{ts_path.stem}.xyz"[:200]
    ts_local.write_text(ts_path.read_text())
    elTS, xyzTS, freqs, modes = run_xtb_hess(
        ts_local, wd / f"hess_{label}", charge=chg, uhf=uhf)
    elTS2, xyzTS2, wboTS = run_xtb(
        ts_local, wd / f"sp_{label}", charge=chg, uhf=uhf)
    g_TS = build_graph(elTS2, wboTS)
    map_R_to_TS = best_mapping(base["g_R"], g_TS, base["wboR"], wboTS)
    core_ts = core_in_TS_frame(base["core_R"], map_R_to_TS)
    ranked = rank_modes(modes, freqs, core_ts)
    # Limit: keep all imag, top-24 real
    imag = [m for m in ranked if m["freq"] < 0]
    real_top = [m for m in ranked if m["freq"] >= 0][:24]
    return {
        "label": label,
        "ts_file": ts_path.name,
        "elements": elTS,
        "coords": xyzTS.tolist(),
        "core_atoms_ts": core_ts,
        "map_R_to_TS": {int(r): int(t) for r, t in map_R_to_TS.items()},
        "n_imag": int((freqs < 0).sum()),
        "freqs": freqs.tolist(),
        "modes_ranked": imag + real_top,
    }


def analyze_step(step_dir, force=False):
    name = step_dir.name
    cache_path = OUT_DIR / f"{name}.json"
    if cache_path.exists() and not force:
        return json.loads(cache_path.read_text()), True   # cached

    chg, uhf = LOOKUP.get(name, (0, 0))
    sanitized = re.sub(r"[^A-Za-z0-9._-]", "_", name)
    wd = WORK / sanitized
    wd.mkdir(parents=True, exist_ok=True)

    rxyz = read_xyzs(step_dir / "reactants")
    pxyz = read_xyzs(step_dir / "products")
    if rxyz is None or pxyz is None:
        raise RuntimeError("missing reactant or product")
    (wd / "reactant.xyz").write_text(rxyz)
    (wd / "product.xyz").write_text(pxyz)

    elR, xyzR, wboR = run_xtb(wd / "reactant.xyz", wd / "R", charge=chg, uhf=uhf)
    elP, xyzP, wboP = run_xtb(wd / "product.xyz", wd / "P", charge=chg, uhf=uhf)
    g_R = build_graph(elR, wboR); g_P = build_graph(elP, wboP)
    map_R_to_P = best_mapping(g_R, g_P, wboR, wboP)
    br, fm, _, _ = classify_bonds(map_R_to_P, wboR, wboP)
    core_R = core_atoms_R_frame(map_R_to_P, br, fm)
    base = dict(g_R=g_R, wboR=wboR, core_R=core_R)

    ts_list = []
    gt = sorted((step_dir / "groundtruth").glob("*.xyz"))
    if gt:
        ts_list.append(process_ts("groundtruth", gt[0], wd, chg, uhf, base))
    for g in sorted(list_initial_guesses(step_dir), key=iter_num):
        ts_list.append(process_ts(f"iter{iter_num(g)}", g, wd, chg, uhf, base))

    data = {
        "step_id":   name,
        "charge":    chg,
        "uhf":       uhf,
        "n_atoms_R": len(elR),
        "core_R":    core_R,
        "broken":    [[int(i), int(j), float(wR), float(wP) if wP is not None else None]
                       for (i, j, wR, wP) in br],
        "formed":    [[int(i), int(j), float(wR) if wR is not None else None, float(wP)]
                       for (i, j, wR, wP) in fm],
        "ts_list":   ts_list,
    }
    cache_path.write_text(json.dumps(data))
    return data, False


PER_STEP_HTML = r"""<!doctype html>
<html><head><meta charset="utf-8">
<title>Mode viewer: __STEP_ID__</title>
<script src="https://3dmol.org/build/3Dmol-min.js"></script>
<style>
 body { font-family: -apple-system, sans-serif; margin: 12px; background: #fafafa; }
 .ctl { background: white; padding: 8px 12px; border: 1px solid #ddd; border-radius: 6px; margin-bottom: 12px; display: flex; gap: 12px; align-items: center; flex-wrap: wrap; }
 select { padding: 4px 6px; font-size: 13px; min-width: 220px; }
 .row { display: flex; gap: 12px; }
 .pane { background: white; border: 1px solid #ddd; border-radius: 6px; padding: 8px; }
 .viewer { width: 600px; height: 520px; position: relative; }
 .modetab { border-collapse: collapse; font-size: 12px; }
 .modetab td, .modetab th { border: 1px solid #ccc; padding: 3px 8px; text-align: right; cursor: pointer; }
 .modetab th { background: #eee; cursor: default; }
 .modetab tr.sel td { background: #d4ecff; font-weight: bold; }
 .meta { color: #444; font-size: 13px; }
 a.back { color: #2368a2; text-decoration: none; }
</style></head><body>

<p><a class="back" href="index.html">← all steps</a></p>
<h2>__STEP_ID__</h2>
<p class="meta">N atoms in R: __N_ATOMS__ | core: __CORE__ | charge/uhf: __CHG_UHF__</p>

<div class="ctl">
  <label><b>TS:</b><select id="tsSel"></select></label>
  <label style="font-size:13px"><input type="checkbox" id="highlightCore" checked> highlight core atoms (orange)</label>
  <span class="meta" id="tsMeta"></span>
</div>

<div class="row">
  <div class="pane">
    <div id="viewer" class="viewer"></div>
  </div>
  <div class="pane">
    <h3 style="margin:6px 0">Modes (imag first, then real; each tier sorted by core fraction)</h3>
    <table class="modetab" id="modeTab">
      <thead><tr><th>rank</th><th>mode#</th><th>ν (cm⁻¹)</th><th>core frac</th></tr></thead>
      <tbody></tbody>
    </table>
  </div>
</div>

<script>
const STEP = __STEP_DATA__;
const tsSel = document.getElementById('tsSel');
STEP.ts_list.forEach((ts, i) => {
  const opt = document.createElement('option');
  opt.value = i;
  opt.textContent = `${ts.label}  (n_imag=${ts.n_imag})`;
  tsSel.appendChild(opt);
});

let viewer = $3Dmol.createViewer('viewer', {backgroundColor: 'white'});

function buildVibFramesXYZ(elements, coords, displacement, n_frames=40, peakAmp=0.7) {
  const n = elements.length;
  let maxLen = 0;
  for (let a = 0; a < n; a++) {
    const d = displacement[a];
    const m = Math.sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);
    if (m > maxLen) maxLen = m;
  }
  if (maxLen < 1e-9) maxLen = 1;
  const scale = peakAmp / maxLen;
  const frames = [];
  for (let f = 0; f < n_frames; f++) {
    const phase = Math.sin(2 * Math.PI * f / n_frames);
    let lines = [String(n), `frame ${f}`];
    for (let a = 0; a < n; a++) {
      const x = coords[a][0] + scale * phase * displacement[a][0];
      const y = coords[a][1] + scale * phase * displacement[a][1];
      const z = coords[a][2] + scale * phase * displacement[a][2];
      lines.push(`${elements[a]}  ${x.toFixed(4)}  ${y.toFixed(4)}  ${z.toFixed(4)}`);
    }
    frames.push(lines.join('\n'));
  }
  return frames.join('\n');
}

function highlightModeRow(modeIdx) {
  const rows = document.querySelectorAll('#modeTab tbody tr');
  rows.forEach(r => r.classList.toggle('sel', +r.dataset.mode === modeIdx));
}

function loadMode(tsIdx, modeRank) {
  const ts = STEP.ts_list[tsIdx];
  if (!ts.modes_ranked || modeRank >= ts.modes_ranked.length) return;
  const mode = ts.modes_ranked[modeRank];
  const framesXyz = buildVibFramesXYZ(ts.elements, ts.coords, mode.displacement);
  if (viewer.isAnimated && viewer.isAnimated()) viewer.pauseAnimate();
  viewer.removeAllModels();
  viewer.removeAllLabels();
  viewer.addModelsAsFrames(framesXyz, 'xyz');
  viewer.setStyle({}, {stick: {radius: 0.10}, sphere: {scale: 0.20}});
  const showCore = document.getElementById('highlightCore').checked;
  if (showCore && ts.core_atoms_ts && ts.core_atoms_ts.length) {
    viewer.setStyle({serial: ts.core_atoms_ts}, {stick: {radius: 0.15, color: 'orange'}, sphere: {scale: 0.32, color: 'orange'}});
  }
  viewer.zoomTo();
  viewer.render();
  viewer.animate({loop: 'forward', interval: 200});
  document.getElementById('tsMeta').textContent =
    `${ts.label}: mode#${mode.mode_idx}  ν=${mode.freq.toFixed(2)} cm⁻¹  core=${mode.core_fraction.toFixed(3)}  (${ts.elements.length} atoms)`;
  highlightModeRow(mode.mode_idx);
}

function rebuildModeTable(tsIdx) {
  const ts = STEP.ts_list[tsIdx];
  const tbody = document.querySelector('#modeTab tbody');
  tbody.innerHTML = '';
  (ts.modes_ranked || []).forEach((m, rank) => {
    const tr = document.createElement('tr');
    tr.dataset.mode = m.mode_idx;
    tr.innerHTML = `<td>${rank+1}</td><td>${m.mode_idx}</td>` +
                   `<td style="color:${m.freq < 0 ? 'red' : 'black'}">${m.freq.toFixed(2)}</td>` +
                   `<td>${m.core_fraction.toFixed(3)}</td>`;
    tr.onclick = () => loadMode(tsIdx, rank);
    tbody.appendChild(tr);
  });
}

tsSel.addEventListener('change', e => {
  const idx = +e.target.value;
  rebuildModeTable(idx);
  loadMode(idx, 0);
});
document.getElementById('highlightCore').addEventListener('change', () => {
  const tsIdx = +tsSel.value;
  const rows = document.querySelectorAll('#modeTab tbody tr');
  let rank = 0;
  rows.forEach((r, i) => { if (r.classList.contains('sel')) rank = i; });
  loadMode(tsIdx, rank);
});
rebuildModeTable(0);
loadMode(0, 0);
</script>
</body></html>
"""


def write_per_step_html(data, out_path):
    html = (PER_STEP_HTML
            .replace("__STEP_ID__", data["step_id"])
            .replace("__N_ATOMS__", str(data["n_atoms_R"]))
            .replace("__CORE__", str(data["core_R"]))
            .replace("__CHG_UHF__", f"{data['charge']}/{data['uhf']}")
            .replace("__STEP_DATA__", json.dumps(data)))
    out_path.write_text(html)


INDEX_HTML = r"""<!doctype html>
<html><head><meta charset="utf-8">
<title>Initial-posing mode viewer — index</title>
<style>
 body { font-family: -apple-system, sans-serif; margin: 20px; background: #fafafa; }
 input { padding: 6px; font-size: 14px; width: 360px; margin-bottom: 12px; }
 table { border-collapse: collapse; font-size: 13px; }
 th, td { border: 1px solid #ccc; padding: 4px 10px; }
 th { background: #eee; }
 a { color: #2368a2; text-decoration: none; }
 a:hover { text-decoration: underline; }
</style></head><body>

<h2>Initial-posing mode viewer</h2>
<p>Per-step pages show all 21 TSes (groundtruth + 20 initial guesses) with vibration modes ranked by core-atom contribution (imag first, then real).</p>
<input id="filter" placeholder="filter (substring)">
<table id="t">
  <thead><tr><th>step</th><th>N atoms (R)</th><th>core size</th><th>charge/uhf</th><th>n TSes</th></tr></thead>
  <tbody>__ROWS__</tbody>
</table>
<script>
document.getElementById('filter').addEventListener('input', e => {
  const f = e.target.value.toLowerCase();
  document.querySelectorAll('#t tbody tr').forEach(r => {
    r.style.display = r.dataset.name.toLowerCase().includes(f) ? '' : 'none';
  });
});
</script>
</body></html>
"""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--start", type=int, default=0)
    ap.add_argument("--steps", nargs="+", default=None)
    ap.add_argument("--force", action="store_true",
                    help="re-run even if cache exists")
    args = ap.parse_args()

    step_dirs = sorted(d for d in BGCP_ROOT.iterdir()
                       if d.is_dir() and not d.name.startswith("."))
    if args.steps:
        wanted = set(args.steps)
        step_dirs = [d for d in step_dirs if d.name in wanted]
    elif args.limit is not None:
        step_dirs = step_dirs[args.start:args.start + args.limit]
    else:
        step_dirs = step_dirs[args.start:]

    print(f"[posing] {len(step_dirs)} steps")
    rows = []
    for k, sd in enumerate(step_dirs, 1):
        t = time.time()
        try:
            data, was_cached = analyze_step(sd, force=args.force)
            html_path = OUT_DIR / f"{sd.name}.html"
            write_per_step_html(data, html_path)
            tag = "CACHED" if was_cached else "OK"
            best = data["ts_list"][0]["modes_ranked"][0] if data["ts_list"] and data["ts_list"][0]["modes_ranked"] else None
            best_f = f"{best['freq']:.0f}" if best else "—"
            best_cf = f"{best['core_fraction']:.2f}" if best else "—"
            print(f"[{k:>3}/{len(step_dirs)}]  {time.time()-t:5.1f}s  {tag:7s}  "
                  f"{sd.name:<60s}  N={data['n_atoms_R']}  core={len(data['core_R'])}  "
                  f"GT_top: ν={best_f} cf={best_cf}")
            rows.append({
                "name": sd.name,
                "n_atoms": data["n_atoms_R"],
                "core_size": len(data["core_R"]),
                "chg_uhf": f"{data['charge']}/{data['uhf']}",
                "n_ts": len(data["ts_list"]),
            })
        except Exception as e:
            print(f"[{k:>3}/{len(step_dirs)}]  FAIL {sd.name}: {e}")
            traceback.print_exc()

    # Index
    rows.sort(key=lambda r: r["name"])
    body = "".join(
        f'<tr data-name="{r["name"]}"><td><a href="{r["name"]}.html">{r["name"]}</a></td>'
        f'<td>{r["n_atoms"]}</td><td>{r["core_size"]}</td>'
        f'<td>{r["chg_uhf"]}</td><td>{r["n_ts"]}</td></tr>'
        for r in rows
    )
    (OUT_DIR / "index.html").write_text(INDEX_HTML.replace("__ROWS__", body))
    print(f"[posing] wrote {OUT_DIR}/index.html  ({len(rows)} steps)")


if __name__ == "__main__":
    main()
