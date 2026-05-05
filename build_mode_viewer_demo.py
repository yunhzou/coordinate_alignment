"""
Single-step demo: vibration-mode viewer with modes ranked by core-atom
contribution. For one BGCP step, process GT + first 3 initial guesses,
compute Hessian, identify core atoms via R/P alignment, project mode
displacements onto core atoms, sort modes by core fraction, animate.

Usage:
  python build_mode_viewer_demo.py [step_id]
  default step: Cyclobutane_JOC2023_updated_TS-CD_step1
"""
from __future__ import annotations
import json
import re
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from rxn_core_frag import (
    run_xtb, run_xtb_hess, build_graph, find_islands, expand_mapping,
    classify_bonds, _generate_seed_orders, write_xyz_str,
)
from build_bgcp_viewer import (
    BGCP_ROOT, OUT, LOOKUP, list_initial_guesses, iter_num, read_xyzs,
)


WORK = Path("/Users/yunhengz/empty_for_claude/rxn_core/work_modedemo")
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
    """Translate core atom indices from R-frame to TS-frame."""
    return sorted(map_R_to_TS[r] for r in core_R if r in map_R_to_TS)


def rank_modes_by_core(modes, freqs, core_indices_ts):
    """Return list of dicts sorted with imaginary frequencies first
    (by core_fraction desc), then real frequencies (by core_fraction
    desc). Imaginary modes are the chemistry-relevant reaction modes
    of a TS, so they should be displayed first regardless of whether
    a real mode happens to have higher core fraction."""
    if modes.shape[0] == 0:
        return []
    sq = (modes ** 2).sum(axis=2)
    total = sq.sum(axis=1)
    if not core_indices_ts:
        return []
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
    # Two-tier sort: imag (freq < 0) first, then real. Within each
    # tier, sort by core fraction descending.
    out.sort(key=lambda d: (0 if d["freq"] < 0 else 1, -d["core_fraction"]))
    return out


def process_ts(label, ts_path, wd, chg, uhf, base):
    """Run xtb hess on TS, align to R, compute mode rankings."""
    ts_local = wd / f"{label}_{ts_path.stem}.xyz"[:120]
    ts_local.write_text(ts_path.read_text())
    elTS, xyzTS, freqs, modes = run_xtb_hess(
        ts_local, wd / f"hess_{label}", charge=chg, uhf=uhf)
    # Need a SP wbo for graph mapping (hess dir also writes wbo but
    # be defensive)
    elTS2, xyzTS2, wboTS = run_xtb(
        ts_local, wd / f"sp_{label}", charge=chg, uhf=uhf)
    g_TS = build_graph(elTS2, wboTS)
    map_R_to_TS = best_mapping(base["g_R"], g_TS, base["wboR"], wboTS)
    core_ts = core_in_TS_frame(base["core_R"], map_R_to_TS)
    ranked = rank_modes_by_core(modes, freqs, core_ts)
    return {
        "label": label,
        "ts_file": ts_path.name,
        "elements": elTS,
        "coords": xyzTS.tolist(),
        "core_atoms_ts": core_ts,
        "n_imag": int((freqs < 0).sum()),
        # Keep all imag modes plus top 24 real modes (so imag are never
        # truncated even when the molecule has many of them).
        "modes": ([m for m in ranked if m["freq"] < 0]
                  + [m for m in ranked if m["freq"] >= 0][:24]),
    }


def analyze_step(step_dir, n_guesses=3):
    name = step_dir.name
    chg, uhf = LOOKUP.get(name, (0, 0))
    sanitized = re.sub(r"[^A-Za-z0-9._-]", "_", name)
    wd = WORK / sanitized
    wd.mkdir(parents=True, exist_ok=True)

    rxyz = read_xyzs(step_dir / "reactants")
    pxyz = read_xyzs(step_dir / "products")
    (wd / "reactant.xyz").write_text(rxyz)
    (wd / "product.xyz").write_text(pxyz)

    elR, xyzR, wboR = run_xtb(wd / "reactant.xyz", wd / "R", charge=chg, uhf=uhf)
    elP, xyzP, wboP = run_xtb(wd / "product.xyz", wd / "P", charge=chg, uhf=uhf)
    g_R = build_graph(elR, wboR); g_P = build_graph(elP, wboP)
    map_R_to_P = best_mapping(g_R, g_P, wboR, wboP)
    br, fm, _, _ = classify_bonds(map_R_to_P, wboR, wboP)
    core_R = core_atoms_R_frame(map_R_to_P, br, fm)

    base = dict(g_R=g_R, wboR=wboR, core_R=core_R)
    ts_results = []

    # Ground truth
    gt = sorted((step_dir / "groundtruth").glob("*.xyz"))
    if gt:
        ts_results.append(process_ts("groundtruth", gt[0], wd, chg, uhf, base))

    # First N initial guesses (lowest iters)
    guesses = sorted(list_initial_guesses(step_dir), key=iter_num)[:n_guesses]
    for g in guesses:
        ts_results.append(process_ts(f"iter{iter_num(g)}", g, wd, chg, uhf, base))

    return {
        "step_id": name,
        "n_atoms_R": len(elR),
        "core_R": core_R,
        "br_fm": (len(br), len(fm)),
        "ts_list": ts_results,
        "charge": chg, "uhf": uhf,
    }


HTML = r"""<!doctype html>
<html><head><meta charset="utf-8">
<title>Mode viewer: {step_id}</title>
<script src="https://3dmol.org/build/3Dmol-min.js"></script>
<style>
 body {{ font-family: -apple-system, sans-serif; margin: 12px; background: #fafafa; }}
 .ctl {{ background: white; padding: 8px 12px; border: 1px solid #ddd; border-radius: 6px; margin-bottom: 12px; display: flex; gap: 12px; align-items: center; flex-wrap: wrap; }}
 select {{ padding: 4px 6px; font-size: 13px; min-width: 220px; }}
 .row {{ display: flex; gap: 12px; }}
 .pane {{ background: white; border: 1px solid #ddd; border-radius: 6px; padding: 8px; }}
 .viewer {{ width: 600px; height: 520px; position: relative; }}
 .modetab {{ border-collapse: collapse; font-size: 12px; }}
 .modetab td, .modetab th {{ border: 1px solid #ccc; padding: 3px 8px; text-align: right; cursor: pointer; }}
 .modetab th {{ background: #eee; cursor: default; }}
 .modetab tr.sel td {{ background: #d4ecff; font-weight: bold; }}
 .meta {{ color: #444; font-size: 13px; }}
</style></head><body>

<h2>Mode viewer — {step_id}</h2>
<p class="meta">N atoms in R: {n_atoms_R} | core: {core_summary} | br/fm: {br_fm_str} | charge/uhf: {chg_uhf}</p>

<div class="ctl">
  <label><b>TS:</b><select id="tsSel"></select></label>
  <label style="font-size:13px"><input type="checkbox" id="highlightCore" checked> highlight core atoms (orange)</label>
  <span class="meta" id="tsMeta"></span>
  <span class="meta" style="margin-left:auto">Click any row in the mode list to select; default = highest core fraction.</span>
</div>

<div class="row">
  <div class="pane">
    <div id="viewer" class="viewer"></div>
  </div>
  <div class="pane">
    <h3 style="margin:6px 0">Modes (sorted by core-atom contribution)</h3>
    <table class="modetab" id="modeTab">
      <thead><tr><th>rank</th><th>mode#</th><th>ν (cm⁻¹)</th><th>core frac</th></tr></thead>
      <tbody></tbody>
    </table>
  </div>
</div>

<script>
const STEP = {step_data};

const tsSel = document.getElementById('tsSel');
STEP.ts_list.forEach((ts, i) => {{
  const opt = document.createElement('option');
  opt.value = i;
  opt.textContent = `${{ts.label}}  (n_imag=${{ts.n_imag}})`;
  tsSel.appendChild(opt);
}});

let viewer = $3Dmol.createViewer('viewer', {{backgroundColor: 'white'}});

function buildVibFramesXYZ(elements, coords, displacement, n_frames=40, peakAmp=0.7) {{
  // Normalize the displacement vector so the largest atom displacement
  // has unit length. Then multiply by peakAmp Å so all modes are
  // visualized at the same physical amplitude regardless of xtb's raw
  // normalization. n_frames=40 + slow interval below gives a smooth
  // visible oscillation; raw mode magnitudes can vary 100x and would
  // otherwise cause some modes to look frantic and others static.
  const n = elements.length;
  let maxLen = 0;
  for (let a = 0; a < n; a++) {{
    const dx = displacement[a][0], dy = displacement[a][1], dz = displacement[a][2];
    const m = Math.sqrt(dx*dx + dy*dy + dz*dz);
    if (m > maxLen) maxLen = m;
  }}
  if (maxLen < 1e-9) maxLen = 1;
  const scale = peakAmp / maxLen;
  const frames = [];
  for (let f = 0; f < n_frames; f++) {{
    const phase = Math.sin(2 * Math.PI * f / n_frames);
    let lines = [String(n), `frame ${{f}}`];
    for (let a = 0; a < n; a++) {{
      const x = coords[a][0] + scale * phase * displacement[a][0];
      const y = coords[a][1] + scale * phase * displacement[a][1];
      const z = coords[a][2] + scale * phase * displacement[a][2];
      lines.push(`${{elements[a]}}  ${{x.toFixed(4)}}  ${{y.toFixed(4)}}  ${{z.toFixed(4)}}`);
    }}
    frames.push(lines.join('\n'));
  }}
  return frames.join('\n');
}}

let currentTSIdx = 0;
let currentModeIdx = 0;

function highlightModeRow(modeIdx) {{
  const rows = document.querySelectorAll('#modeTab tbody tr');
  rows.forEach(r => r.classList.toggle('sel', +r.dataset.mode === modeIdx));
}}

function loadMode(tsIdx, modeRank) {{
  const ts = STEP.ts_list[tsIdx];
  if (modeRank >= ts.modes.length) return;
  const mode = ts.modes[modeRank];
  currentTSIdx = tsIdx;
  currentModeIdx = mode.mode_idx;

  const framesXyz = buildVibFramesXYZ(ts.elements, ts.coords, mode.displacement);
  // Pause any running animation before reloading the model. Without
  // this, the previous animation timer keeps firing render() while
  // the new frames are being installed, producing a blinking flash.
  if (viewer.isAnimated && viewer.isAnimated()) {{
    viewer.pauseAnimate();
  }}
  viewer.removeAllModels();
  viewer.removeAllLabels();
  viewer.addModelsAsFrames(framesXyz, 'xyz');
  viewer.setStyle({{}}, {{stick: {{radius: 0.10}}, sphere: {{scale: 0.20}}}});
  // Highlight core atoms (toggleable). When highlight is off, core
  // atoms render with the default coloring like everything else.
  const showCore = document.getElementById('highlightCore').checked;
  if (showCore && ts.core_atoms_ts && ts.core_atoms_ts.length) {{
    viewer.setStyle({{serial: ts.core_atoms_ts}}, {{stick: {{radius: 0.15, color: 'orange'}}, sphere: {{scale: 0.32, color: 'orange'}}}});
  }}
  viewer.zoomTo();
  // Render once to draw the new structure cleanly, then start the
  // animation. Order matters: animate() schedules its first frame on
  // a timer, but the model has to already be drawn or we get a
  // blink between the cleared scene and the first animate frame.
  viewer.render();
  // Slow animation: 40 frames * 200ms = 8s per oscillation cycle.
  // The displacement is also normalized so each mode reads at the same
  // visible amplitude, independent of xtb's raw mode magnitudes.
  viewer.animate({{loop: 'forward', interval: 200}});

  document.getElementById('tsMeta').textContent =
    `${{ts.label}}: mode#${{mode.mode_idx}}  ν=${{mode.freq.toFixed(2)}} cm⁻¹  core=${{mode.core_fraction.toFixed(3)}}  (${{ts.elements.length}} atoms)`;
  highlightModeRow(mode.mode_idx);
}}

function rebuildModeTable(tsIdx) {{
  const ts = STEP.ts_list[tsIdx];
  const tbody = document.querySelector('#modeTab tbody');
  tbody.innerHTML = '';
  ts.modes.forEach((m, rank) => {{
    const tr = document.createElement('tr');
    tr.dataset.mode = m.mode_idx;
    tr.dataset.rank = rank;
    tr.innerHTML = `<td>${{rank+1}}</td><td>${{m.mode_idx}}</td>` +
                   `<td style="color:${{m.freq < 0 ? 'red' : 'black'}}">${{m.freq.toFixed(2)}}</td>` +
                   `<td>${{m.core_fraction.toFixed(3)}}</td>`;
    tr.onclick = () => loadMode(tsIdx, rank);
    tbody.appendChild(tr);
  }});
}}

tsSel.addEventListener('change', e => {{
  const idx = +e.target.value;
  rebuildModeTable(idx);
  loadMode(idx, 0);
}});
document.getElementById('highlightCore').addEventListener('change', () => {{
  // Reload current mode to apply new highlight state
  const tsIdx = +tsSel.value;
  const rows = document.querySelectorAll('#modeTab tbody tr');
  let rank = 0;
  rows.forEach((r, i) => {{ if (r.classList.contains('sel')) rank = i; }});
  loadMode(tsIdx, rank);
}});

rebuildModeTable(0);
loadMode(0, 0);
</script>
</body></html>
"""


def main():
    step_id = sys.argv[1] if len(sys.argv) > 1 else "Cyclobutane_JOC2023_updated_TS-CD_step1"
    n_guesses = int(sys.argv[2]) if len(sys.argv) > 2 else 20

    step_dir = BGCP_ROOT / step_id
    if not step_dir.is_dir():
        raise SystemExit(f"step dir not found: {step_dir}")

    print(f"[demo] processing {step_id} (GT + first {n_guesses} guesses)")
    data = analyze_step(step_dir, n_guesses=n_guesses)

    html = HTML.format(
        step_id=step_id,
        n_atoms_R=data["n_atoms_R"],
        core_summary=str(data["core_R"]),
        br_fm_str=f"{data['br_fm'][0]}/{data['br_fm'][1]}",
        chg_uhf=f"{data['charge']}/{data['uhf']}",
        step_data=json.dumps(data),
    )
    out = OUT / f"mode_demo_{step_id}.html"
    out.write_text(html)
    print(f"[demo] wrote {out}")
    print(f"[demo] Top modes per TS:")
    for ts in data["ts_list"]:
        print(f"  {ts['label']:<14s}  n_imag={ts['n_imag']}  best_mode_freq="
              f"{ts['modes'][0]['freq']:>8.1f}  core_frac={ts['modes'][0]['core_fraction']:.3f}")


if __name__ == "__main__":
    main()
