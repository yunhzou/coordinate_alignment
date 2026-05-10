"""
Regenerate just out/bgcp_views/index.html with extra oracle columns.
Reads from out/bgcp_alignment_eval.json (no alignment re-run).
"""
from __future__ import annotations
import json
from pathlib import Path

PROJECT = Path(__file__).resolve().parent
JSON = PROJECT / "out" / "bgcp_alignment_eval.json"
OUT_ROOT = PROJECT / "out" / "bgcp_views"

W_RXN, W_CORE, IMAG_PEN = 1.0, 0.2, 0.3
def full_score(r):
    n = max(r.get("n_imag", 0), 1)
    return (r.get("beta", 0.0) * (1 + W_RXN * r.get("rho", 0.0))
            * (1 + W_CORE * r.get("kappa", 0.0)) / n ** IMAG_PEN)


def main():
    data = json.loads(JSON.read_text())
    rows = []
    for rec in data:
        if rec.get("error"): continue
        igs = [g for g in rec.get("igs", []) if "beta" in g]
        if not igs: continue
        # Verifier's top-1 = highest S
        sorted_by_S = sorted(igs, key=lambda r: -full_score(r))
        top1 = sorted_by_S[0]
        # Oracle = highest mwc
        sorted_by_mwc = sorted(igs, key=lambda r: -r["mwc_to_GT"])
        oracle = sorted_by_mwc[0]
        rows.append({
            "step": rec["step"],
            "n_ig": len(igs),
            "top1_label": top1["label"],
            "top1_S":     full_score(top1),
            "top1_mwc":   top1["mwc_to_GT"],
            "oracle_label": oracle["label"],
            "oracle_S":   full_score(oracle),
            "oracle_mwc": oracle["mwc_to_GT"],
        })
    rows.sort(key=lambda r: r["step"])

    rows_html = "".join(
        f"<tr><td><a href='{r['step']}/view.html'>{r['step']}</a></td>"
        f"<td>{r['n_ig']}</td>"
        f"<td>{r['top1_label']}</td>"
        f"<td>{r['top1_S']:.3f}</td>"
        f"<td>{r['top1_mwc']:.3f}</td>"
        f"<td>{r['oracle_label']}</td>"
        f"<td>{r['oracle_S']:.3f}</td>"
        f"<td>{r['oracle_mwc']:.3f}</td>"
        f"<td>{r['oracle_mwc'] - r['top1_mwc']:+.3f}</td></tr>"
        for r in rows
    )
    n = len(rows)
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    (OUT_ROOT / "index.html").write_text(f"""<!doctype html><html><head><meta charset="utf-8">
<title>BGCP ranked views</title>
<style>
body{{font-family:-apple-system,sans-serif;margin:20px;max-width:1400px}}
table{{border-collapse:collapse}} th,td{{border:1px solid #ccc;padding:6px 10px;font-size:13px}}
caption{{caption-side:top;text-align:left;font-size:14px;padding:6px 0;font-weight:600}}
.note{{color:#666;font-size:12px;margin:8px 0 14px}}
</style></head><body>
<h2>BGCP ranked views ({n} steps)</h2>
<p class="note">
  <b>top1</b>: the verifier's rank-1 IG (highest S).
  <b>oracle</b>: the IG with the highest mwc-to-GT (best alignment with GT, regardless of ranker).
  <b>S</b>: ranker score &beta;(1+&rho;)(1+0.2&kappa;)/n_imag^0.3.
  <b>mwc</b>: mass-weighted cosine of picked imag mode vs GT's picked imag mode (informational; not a substitute for chemist judgment).
  <b>&Delta;mwc</b> = oracle_mwc &minus; top1_mwc; large positive means a better-aligned IG exists in the pool that the ranker missed.
</p>
<table>
<tr>
  <th>step</th><th>n_ig</th>
  <th>top1<br>label</th><th>top1<br>S</th><th>top1<br>mwc</th>
  <th>oracle<br>label</th><th>oracle<br>S</th><th>oracle<br>mwc</th>
  <th>&Delta;mwc</th>
</tr>
{rows_html}
</table></body></html>""")
    print(f"wrote {OUT_ROOT / 'index.html'}  ({n} rows)")
    # Quick stats
    miss = sum(1 for r in rows if r["oracle_mwc"] - r["top1_mwc"] > 0.2)
    perfect_oracle = sum(1 for r in rows if r["oracle_mwc"] >= 0.95)
    perfect_top1 = sum(1 for r in rows if r["top1_mwc"] >= 0.95)
    print(f"  oracle mwc >= 0.95: {perfect_oracle}/{n}")
    print(f"  top1 mwc >= 0.95: {perfect_top1}/{n}")
    print(f"  ranker missed (oracle - top1 > 0.2): {miss}/{n}")


if __name__ == "__main__":
    main()
