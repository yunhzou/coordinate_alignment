"""
Guess-quality table: compare each initial-guess TS's default reaction
mode to its groundtruth's default mode by 3N-Cartesian cosine
similarity. This is an EVALUATION metric (uses GT) downstream of
the bond_overlap ranking (which doesn't use GT).

For every BGCP step:
  GT mode  = GT TS's imag mode with highest bond_overlap
  IG mode  = each IG TS's imag mode with highest bond_overlap (its own)
  alignment = |d_GT · d_IG| / (||d_GT|| · ||d_IG||)

Both modes are already reindexed into R-atom-index frame in the
per-step viewer HTMLs, so atom i is the same chemical atom in both.

Outputs:
  out/mode_viewer/guess_quality.html — combined table with sortable
                                        rows, per-step expansion
  out/mode_analysis/guess_quality.csv — long-format rows
"""
from __future__ import annotations
import csv
import json
import re
import time
from pathlib import Path

import numpy as np


SRC_DIR = Path(__file__).parent / "out" / "mode_viewer"
OUT_HTML = SRC_DIR / "guess_quality.html"
OUT_CSV = Path(__file__).parent / "out" / "mode_analysis" / "guess_quality.csv"


def cos_sim(a, b):
    """|a·b| / (||a||·||b||) for 3D arrays, sign-blind. ∈ [0,1]."""
    a = np.asarray(a).reshape(-1)
    b = np.asarray(b).reshape(-1)
    na = float(np.linalg.norm(a)); nb = float(np.linalg.norm(b))
    if na < 1e-9 or nb < 1e-9:
        return 0.0
    return abs(float(a @ b)) / (na * nb)


def load_step_payload(html_path):
    text = html_path.read_text()
    m = re.search(r"const DATA = (\{.*?\});\n", text, re.DOTALL)
    if not m:
        return None
    return json.loads(m.group(1))


def step_rows(payload):
    """For one step: list of dicts, one per IG TS, with comparison
    against the GT TS. Skips steps with no GT or no imag GT mode."""
    ts_by_label = {ts['label']: ts for ts in payload['ts_list']}
    gt = ts_by_label.get('groundtruth')
    if gt is None or not gt.get('modes'):
        return []
    gt_default_idx = gt.get('default_mode_idx', 0)
    gt_mode = gt['modes'][gt_default_idx]
    if gt_mode['freq'] >= 0:
        return []  # GT default isn't an imag mode (TS is bad)
    gt_disp = np.asarray(gt_mode['disp'])
    rows = []
    for ts in payload['ts_list']:
        if ts['label'] == 'groundtruth' or not ts.get('modes'):
            continue
        ig_default_idx = ts.get('default_mode_idx', 0)
        ig_mode = ts['modes'][ig_default_idx]
        ig_disp = np.asarray(ig_mode['disp'])
        align = cos_sim(ig_disp, gt_disp)
        rows.append({
            'step': payload['step'],
            'ig_label': ts['label'],
            'ig_freq': ig_mode['freq'],
            'ig_is_imag': ig_mode['freq'] < 0,
            'ig_bond_ovlp': ig_mode.get('bond_overlap', 0.0),
            'ig_rxn_ovlp':  ig_mode.get('rxn_overlap', 0.0),
            'ig_n_imag': ts['n_imag'],
            'gt_freq': gt_mode['freq'],
            'gt_bond_ovlp': gt_mode.get('bond_overlap', 0.0),
            'gt_n_imag': gt['n_imag'],
            'gt_alignment': align,
            'n_atoms': payload['n_atoms'],
            'n_core': len(payload['core_atoms']),
        })
    return rows


HTML = r"""<!doctype html>
<html><head><meta charset="utf-8">
<title>BGCP — initial-guess quality (vs GT reaction-mode direction)</title>
<style>
 body { font-family: -apple-system, sans-serif; margin: 16px; max-width: 1300px; }
 .ctl { display: flex; gap: 12px; align-items: center; flex-wrap: wrap; margin-bottom: 12px; }
 input { padding: 4px 6px; font-size: 13px; }
 select { padding: 4px 6px; font-size: 13px; }
 table { border-collapse: collapse; font-size: 13px; width: 100%; }
 th, td { border: 1px solid #ccc; padding: 4px 8px; text-align: left; }
 th { background: #eee; cursor: pointer; }
 th.asc::after  { content: ' ▲'; }
 th.desc::after { content: ' ▼'; }
 tr.imag-no { background: #fff0e0; }
 .bar { display: inline-block; height: 10px; background: #4a90e2; vertical-align: middle; }
 .gt-mark { background: #fffacc; font-weight: 600; }
 a { color: #06c; text-decoration: none; }
 a:hover { text-decoration: underline; }
 details { margin-bottom: 8px; }
 details summary { padding: 4px 8px; background: #f4f4f4; border: 1px solid #ddd; border-radius: 4px; cursor: pointer; font-weight: 600; }
</style></head><body>

<h2>Guess quality: cosine alignment of each IG's default reaction mode vs the GT's</h2>
<p>For every BGCP step, the GT TS's default mode (highest-<code>bond_overlap</code>
imaginary mode) is the reference. Each initial guess TS contributes its own
default mode; <code>gt_alignment</code> is |cos similarity| in 3N-Cartesian
displacement space. <code>1.0</code> = same direction as GT, <code>0.0</code>
= orthogonal. Sign-blind (modes are sign-arbitrary).</p>

<div class="ctl">
  <input id="filter" placeholder="filter by step or label" oninput="filterRows()">
  <label>Min gt_alignment <input type="number" id="minA" min="0" max="1" step="0.05" value="0" oninput="filterRows()" style="width:60px"></label>
  <label>Min imag count for IG <input type="number" id="minI" min="0" max="20" step="1" value="0" oninput="filterRows()" style="width:60px"></label>
  <span class="ctl"><a href="index.html">← all steps</a> · <a href="flat_view.html">flat view</a></span>
</div>

<h3>Per-IG long-format table (sortable: click headers)</h3>
<table id="long"><thead><tr>
<th data-k="step">step</th>
<th data-k="ig_label">IG</th>
<th data-k="ig_freq" data-num="1">IG freq</th>
<th data-k="ig_n_imag" data-num="1">IG n_imag</th>
<th data-k="ig_bond_ovlp" data-num="1">IG bond_ovlp</th>
<th data-k="gt_freq" data-num="1">GT freq</th>
<th data-k="gt_bond_ovlp" data-num="1">GT bond_ovlp</th>
<th data-k="gt_alignment" data-num="1">gt_alignment</th>
</tr></thead><tbody id="longBody">__LONG_ROWS__</tbody></table>

<script>
const ROWS = __JSON_ROWS__;

function fmt(n) { return (typeof n === 'number') ? n.toFixed(3) : n; }
function bar(v) {
  if (typeof v !== 'number') return '';
  const w = Math.max(0, Math.min(1, v)) * 80;
  return `<span class="bar" style="width:${w}px"></span> ${v.toFixed(3)}`;
}

function rowHtml(r) {
  const cls = r.ig_is_imag ? '' : ' class="imag-no"';
  const stepLink = `<a href="${r.step.replace(/[^A-Za-z0-9._-]/g, '_')}.html" target="_blank">${r.step}</a>`;
  return `<tr${cls} data-step="${r.step}" data-label="${r.ig_label}" data-align="${r.gt_alignment}" data-imag="${r.ig_n_imag}">
    <td>${stepLink}</td>
    <td>${r.ig_label}</td>
    <td>${r.ig_freq.toFixed(2)}</td>
    <td>${r.ig_n_imag}</td>
    <td>${r.ig_bond_ovlp.toFixed(3)}</td>
    <td>${r.gt_freq.toFixed(2)}</td>
    <td>${r.gt_bond_ovlp.toFixed(3)}</td>
    <td>${bar(r.gt_alignment)}</td></tr>`;
}

function render(rows) {
  document.getElementById('longBody').innerHTML = rows.map(rowHtml).join('');
}

let curRows = ROWS.slice();
let sortKey = 'gt_alignment';
let sortDesc = true;

function applySort() {
  curRows.sort((a, b) => {
    const av = a[sortKey], bv = b[sortKey];
    if (typeof av === 'string') return sortDesc ? bv.localeCompare(av) : av.localeCompare(bv);
    return sortDesc ? (bv - av) : (av - bv);
  });
}
applySort();

function filterRows() {
  const f = document.getElementById('filter').value.toLowerCase();
  const minA = parseFloat(document.getElementById('minA').value || '0');
  const minI = parseInt(document.getElementById('minI').value || '0');
  const filt = ROWS.filter(r => {
    if (f && !((r.step + ' ' + r.ig_label).toLowerCase().includes(f))) return false;
    if (r.gt_alignment < minA) return false;
    if (r.ig_n_imag < minI) return false;
    return true;
  });
  curRows = filt.slice(); applySort(); render(curRows);
}

document.querySelectorAll('th').forEach(th => {
  th.addEventListener('click', () => {
    const k = th.dataset.k;
    if (sortKey === k) sortDesc = !sortDesc; else { sortKey = k; sortDesc = !!th.dataset.num; }
    document.querySelectorAll('th').forEach(x => { x.classList.remove('asc'); x.classList.remove('desc'); });
    th.classList.add(sortDesc ? 'desc' : 'asc');
    applySort(); render(curRows);
  });
});

render(curRows);
</script>

</body></html>
"""


def main():
    files = sorted(SRC_DIR.glob("*.html"))
    files = [f for f in files
             if f.name not in ('index.html', 'flat_view.html', 'guess_quality.html')]
    print(f"Computing guess-quality from {len(files)} per-step HTMLs...")
    all_rows = []
    t0 = time.time()
    for hp in files:
        payload = load_step_payload(hp)
        if payload is None:
            continue
        rows = step_rows(payload)
        all_rows.extend(rows)

    # CSV
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    if all_rows:
        with OUT_CSV.open('w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=list(all_rows[0].keys()))
            w.writeheader(); w.writerows(all_rows)

    # HTML — embed pre-rendered rows + JSON for client-side sort/filter
    json_rows = json.dumps(all_rows)
    long_html = "".join(
        f'<tr><td>{r["step"]}</td><td>{r["ig_label"]}</td>'
        f'<td>{r["ig_freq"]:.2f}</td><td>{r["ig_n_imag"]}</td>'
        f'<td>{r["ig_bond_ovlp"]:.3f}</td>'
        f'<td>{r["gt_freq"]:.2f}</td><td>{r["gt_bond_ovlp"]:.3f}</td>'
        f'<td>{r["gt_alignment"]:.3f}</td></tr>'
        for r in sorted(all_rows, key=lambda x: -x['gt_alignment']))
    OUT_HTML.write_text(
        HTML.replace('__JSON_ROWS__', json_rows).replace('__LONG_ROWS__', long_html))

    # Quick stats
    n_steps = len({r['step'] for r in all_rows})
    aligns = [r['gt_alignment'] for r in all_rows]
    print(f"\nDone in {time.time()-t0:.1f}s. {len(all_rows)} (step, IG) rows across {n_steps} steps.")
    if aligns:
        a = np.asarray(aligns)
        print(f"gt_alignment: mean={a.mean():.3f}  median={np.median(a):.3f}  "
              f"max={a.max():.3f}  min={a.min():.3f}")
        for thr in [0.9, 0.7, 0.5, 0.3, 0.1]:
            print(f"  ≥{thr}: {(a >= thr).sum():>4} / {len(a)} IGs ({100*(a >= thr).mean():.1f}%)")
    print(f"\nCSV:  {OUT_CSV}")
    print(f"HTML: {OUT_HTML}")


if __name__ == "__main__":
    main()
