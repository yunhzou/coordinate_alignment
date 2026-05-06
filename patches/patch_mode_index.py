"""Rebuild out/mode_viewer/index.html from the patched per-step HTMLs.
Sorts steps by best imag bond_overlap across their TS structures."""
from __future__ import annotations
import sys as _sys
from pathlib import Path as _Path
_HERE = _Path(__file__).resolve().parent
_sys.path.insert(0, str(_HERE.parent / "src"))
_sys.path.insert(0, str(_HERE))
PROJECT_ROOT = _HERE.parent  # _RXN_CORE_PATH_SETUP

import json
import re
from pathlib import Path

OUT_DIR = PROJECT_ROOT / "out" / "mode_viewer"


INDEX_HTML = r"""<!doctype html>
<html><head><meta charset="utf-8"><title>BGCP vibration mode viewer</title>
<style>
 body { font-family: -apple-system, sans-serif; margin: 20px; max-width: 1100px; }
 input { padding: 4px 6px; font-size: 13px; width: 320px; }
 table { border-collapse: collapse; font-size: 13px; width: 100%; margin-top: 12px; }
 td, th { border: 1px solid #ccc; padding: 4px 8px; text-align: left; }
 a { color: #06c; text-decoration: none; }
 a:hover { text-decoration: underline; }
</style></head><body>
<h2>BGCP — vibration modes ranked by bond_overlap (concerted bond breaking/forming alignment)</h2>
<p>Per-step page shows TS structure, broken (red) / formed (green) bonds,
and animated vibration of the selected mode. Default selection is the
imaginary mode whose displacement projects most strongly onto the
bond-stretching/-compressing direction at TS coords (= the actual
reaction coordinate, weighted by mode amplitude).</p>

<input id="filter" placeholder="filter (substring)" oninput="filt()">

<table id="t"><thead><tr>
<th>step</th><th>N</th><th>core</th><th>n_TS</th><th>n_TS_with_imag</th>
<th>best imag bond_overlap</th><th>best label</th>
</tr></thead><tbody>
__ROWS__
</tbody></table>

<script>
function filt() {
  const f = document.getElementById('filter').value.toLowerCase();
  document.querySelectorAll('#t tbody tr').forEach(tr => {
    tr.style.display = tr.dataset.name.toLowerCase().includes(f) ? '' : 'none';
  });
}
</script>
</body></html>
"""


def main():
    rows = []
    for hp in sorted(OUT_DIR.glob("*.html")):
        if hp.name in ('index.html', 'flat_view.html'): continue
        text = hp.read_text()
        m = re.search(r"const DATA = (\{.*?\});\n", text, re.DOTALL)
        if not m: continue
        data = json.loads(m.group(1))
        best = 0.0; best_label = '—'; n_with_imag = 0
        for ts in data['ts_list']:
            for mm in ts['modes']:
                if mm['freq'] < 0 and mm.get('bond_overlap', 0) > best:
                    best = mm['bond_overlap']; best_label = ts['label']
            if ts['n_imag'] > 0: n_with_imag += 1
        rows.append({
            'name': data['step'], 'file': hp.name,
            'n_atoms': data['n_atoms'], 'n_core': len(data['core_atoms']),
            'n_ts': len(data['ts_list']), 'n_with_imag': n_with_imag,
            'best': best, 'best_label': best_label,
        })
    rows.sort(key=lambda r: -r['best'])
    body = ""
    for r in rows:
        body += (f"<tr data-name='{r['name']}'>"
                 f"<td><a href='{r['file']}' target='_blank'>{r['name']}</a></td>"
                 f"<td>{r['n_atoms']}</td><td>{r['n_core']}</td>"
                 f"<td>{r['n_ts']}</td><td>{r['n_with_imag']}</td>"
                 f"<td>{r['best']:.3f}</td><td>{r['best_label']}</td></tr>")
    (OUT_DIR / "index.html").write_text(INDEX_HTML.replace('__ROWS__', body))
    print(f"Index regenerated: {OUT_DIR / 'index.html'}  ({len(rows)} steps)")


if __name__ == "__main__":
    main()
