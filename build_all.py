"""
Run analysis on every benchmark step that has plain/stage0/{reactant,product}.xyz,
write per-step HTML, then build a single index.html that links to all of them.
"""

from __future__ import annotations
import json
import re
import sys
import time
import traceback
from pathlib import Path

from rxn_core_wbo import analyze
from viz_3dmol import render_html


BENCH = Path("/Users/yunhengz/empty_for_claude/Benchmark")
ROOT = Path("/Users/yunhengz/empty_for_claude/rxn_core")
OUT = ROOT / "out"
WORK = ROOT / "work"
OUT.mkdir(parents=True, exist_ok=True)


def parse_charge_uhf(xyz_path):
    """Best-effort parse of 'charge=N' / 'multiplicity=M' from xyz title line."""
    title = Path(xyz_path).read_text().splitlines()[1] if Path(xyz_path).exists() else ""
    chg = 0
    uhf = 0
    m = re.search(r"charge\s*=\s*(-?\d+)", title)
    if m:
        chg = int(m.group(1))
    m = re.search(r"multiplicity\s*=\s*(\d+)", title)
    if m:
        uhf = max(0, int(m.group(1)) - 1)
    return chg, uhf


def list_steps():
    steps = []
    for d in sorted(BENCH.iterdir()):
        if not d.is_dir():
            continue
        r = d / "plain" / "stage0" / "reactant.xyz"
        p = d / "plain" / "stage0" / "product.xyz"
        if r.exists() and p.exists():
            steps.append((d.name, r, p))
    return steps


def run_one(name, r, p):
    chg_r, uhf_r = parse_charge_uhf(r)
    chg_p, uhf_p = parse_charge_uhf(p)
    # use reactant's; mismatched charge/multiplicity is its own kind of failure
    chg, uhf = chg_r, uhf_r
    res = analyze(r, p, WORK / name, charge=chg, uhf=uhf)
    render_html(res, title=name, out_path=OUT / f"{name}.html")
    return {
        "step": name, "status": "ok",
        "natoms": len(res["elements_R"]),
        "anchors": len(res["anchors"]),
        "spectator": res["n_spectator"],
        "mapped": len(res["mapping"]),
        "broken": len(res["broken"]),
        "formed": len(res["formed"]),
        "core_R": len(res["core_R"]),
        "core_P": len(res["core_P"]),
        "charge": chg,
        "uhf": uhf,
    }


def write_index(records):
    OK = [r for r in records if r["status"] == "ok"]
    FAIL = [r for r in records if r["status"] != "ok"]
    rows = []
    for r in OK:
        coverage = r["spectator"] / max(r["natoms"], 1)
        bad_flag = ""
        if r["broken"] != r["formed"]:
            bad_flag = "asym"
        if r["spectator"] < 0.3 * r["natoms"]:
            bad_flag = "low-coverage"
        if r["broken"] + r["formed"] == 0:
            bad_flag = "no-change"
        cls = ""
        if bad_flag:
            cls = " warn"
        rows.append(
            f'<tr class="row{cls}" data-natoms="{r["natoms"]}" data-broken="{r["broken"]}" '
            f'data-formed="{r["formed"]}" data-coverage="{coverage:.2f}">'
            f'<td><a href="{r["step"]}.html">{r["step"]}</a></td>'
            f'<td>{r["natoms"]}</td>'
            f'<td>{r["anchors"]}</td>'
            f'<td>{r["spectator"]} ({coverage*100:.0f}%)</td>'
            f'<td>{r["broken"]}</td>'
            f'<td>{r["formed"]}</td>'
            f'<td>{r["core_R"]}/{r["core_P"]}</td>'
            f'<td>{r["charge"]}/{r["uhf"]}</td>'
            f'<td>{bad_flag}</td>'
            f'</tr>'
        )
    fail_rows = []
    for r in FAIL:
        msg = r.get("error", "?")[:200]
        fail_rows.append(
            f'<tr><td>{r["step"]}</td><td colspan="8" class="err">{msg}</td></tr>'
        )

    html = f"""<!doctype html>
<html><head><meta charset="utf-8"><title>Reaction-core analyses</title>
<style>
 body {{ font-family: -apple-system, sans-serif; margin: 16px; background: #fafafa; }}
 h2 {{ margin: 6px 0; }}
 table {{ border-collapse: collapse; width: 100%; background: white; font-size: 12px; }}
 th, td {{ border: 1px solid #ddd; padding: 4px 8px; text-align: left; }}
 th {{ background: #eee; cursor: pointer; user-select: none; position: sticky; top: 0; }}
 tr.warn {{ background: #fff7e0; }}
 td.err {{ color: #a00; font-family: monospace; font-size: 11px; }}
 .stats {{ background: white; padding: 8px 12px; border: 1px solid #ddd; margin-bottom: 12px; border-radius: 6px; }}
 input {{ padding: 4px; font-size: 13px; width: 320px; }}
</style></head><body>
<h2>Reaction-core analyses</h2>
<div class="stats">
  Total: {len(records)} | OK: {len(OK)} | Failed: {len(FAIL)} |
  Total broken bonds: {sum(r["broken"] for r in OK)} |
  Total formed bonds: {sum(r["formed"] for r in OK)}
</div>
<input id="filter" placeholder="filter by step name...">
<table id="t">
<thead><tr>
 <th onclick="sortBy(0,'s')">Step</th>
 <th onclick="sortBy(1,'n')">N atoms</th>
 <th onclick="sortBy(2,'n')">Anchors</th>
 <th onclick="sortBy(3,'n')">Spectator (cov)</th>
 <th onclick="sortBy(4,'n')">Broken</th>
 <th onclick="sortBy(5,'n')">Formed</th>
 <th onclick="sortBy(6,'s')">Core R/P</th>
 <th onclick="sortBy(7,'s')">chg/uhf</th>
 <th onclick="sortBy(8,'s')">Flag</th>
</tr></thead>
<tbody>
{''.join(rows)}
</tbody>
</table>

<h3 style="margin-top: 24px">Failed steps ({len(FAIL)})</h3>
<table>
<thead><tr><th>Step</th><th colspan="8">Error</th></tr></thead>
<tbody>{''.join(fail_rows)}</tbody>
</table>

<script>
let asc = true;
function sortBy(col, type) {{
  const tb = document.querySelector('#t tbody');
  const rows = [...tb.querySelectorAll('tr')];
  rows.sort((a, b) => {{
    let A = a.cells[col].textContent.trim();
    let B = b.cells[col].textContent.trim();
    if (type === 'n') {{ A = parseFloat(A) || 0; B = parseFloat(B) || 0; }}
    return asc ? (A > B ? 1 : -1) : (A < B ? 1 : -1);
  }});
  asc = !asc;
  rows.forEach(r => tb.appendChild(r));
}}
document.getElementById('filter').addEventListener('input', e => {{
  const q = e.target.value.toLowerCase();
  document.querySelectorAll('#t tbody tr').forEach(r => {{
    r.style.display = r.cells[0].textContent.toLowerCase().includes(q) ? '' : 'none';
  }});
}});
</script>
</body></html>
"""
    (OUT / "index.html").write_text(html)


def main():
    steps = list_steps()
    print(f"[build_all] {len(steps)} steps to process")
    records = []
    t0 = time.time()
    for k, (name, r, p) in enumerate(steps, 1):
        try:
            rec = run_one(name, r, p)
            records.append(rec)
            print(f"[{k:>3}/{len(steps)}] OK    {name:<60s} "
                  f"spec={rec['spectator']}/{rec['natoms']} "
                  f"broken={rec['broken']} formed={rec['formed']}")
        except Exception as e:
            tb = traceback.format_exc()
            records.append({"step": name, "status": "fail", "error": str(e)})
            print(f"[{k:>3}/{len(steps)}] FAIL  {name:<60s} {str(e)[:120]}")
        if k % 10 == 0:
            (OUT / "progress.json").write_text(json.dumps(records, indent=2))
            write_index(records)
    (OUT / "results.json").write_text(json.dumps(records, indent=2))
    write_index(records)
    print(f"[build_all] done in {time.time()-t0:.0f}s. Index: {OUT}/index.html")


if __name__ == "__main__":
    main()
