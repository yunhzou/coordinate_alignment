"""
For each of the 7 PQ-vs-OLD regression cases, generate 10 random-seed
trace HTMLs of the OLD algorithm (so we can compare how seed ordering
shapes the mapping). Outputs grouped per-step:

  out/regressions/<step>/seed_0.html
  out/regressions/<step>/seed_1.html
  ...
  out/regressions/<step>/seed_9.html
  out/regressions/<step>/index.html

  out/regressions/index.html       — top-level link to each step

Usage: python build_regression_traces.py
"""
from __future__ import annotations
import json
import random
import re
from pathlib import Path

from rxn_core_frag import (
    run_xtb, build_graph, write_xyz_str, classify_bonds,
)
from trace_run import find_islands_with_trace, HTML, _resolve_step_inputs


REGRESSIONS = [
    "pr7.V.dodh_ts56-triplet",
    "pr7.V.dodh_ts56-singlet",
    "pr7.V.dodh_ts71",
    "pr7.V.dodh_ts1314",
    "Jackie_TS_10",
    "pr19.heck_ts1",
    "pr16.carbocation_ts5",
]

OUT_ROOT = Path(__file__).parent / "out" / "regressions"
N_SEEDS = 10
RNG_SEED = 42


def run_step(step):
    sanitized = re.sub(r"[^A-Za-z0-9._-]", "_", step)
    out_dir = OUT_ROOT / sanitized
    out_dir.mkdir(parents=True, exist_ok=True)

    rxyz_path, pxyz_path, work, chg, uhf = _resolve_step_inputs(step)
    elR, xyzR_arr, wboR = run_xtb(rxyz_path, work / "R", charge=chg, uhf=uhf)
    elP, xyzP_arr, wboP = run_xtb(pxyz_path, work / "P", charge=chg, uhf=uhf)
    g_R = build_graph(elR, wboR)
    g_P = build_graph(elP, wboP)

    nodes = list(g_R.nodes())
    rng = random.Random(RNG_SEED)
    seed_orders = []
    for _ in range(N_SEEDS):
        perm = list(nodes); rng.shuffle(perm); seed_orders.append(perm)

    summaries = []
    for i, order in enumerate(seed_orders):
        mapping, events = find_islands_with_trace(g_R, g_P, seed_order=order)
        br, fm, _, _ = classify_bonds(mapping, wboR, wboP)
        title = (f"{step}  seed#{i}  br/fm={len(br)}/{len(fm)}  "
                 f"events={len(events)}  mapped={len(mapping)}/{len(elR)}")
        html = HTML.format(
            title=title,
            xyzR_json=json.dumps(write_xyz_str(elR, xyzR_arr, comment="reactant")),
            xyzP_json=json.dumps(write_xyz_str(elP, xyzP_arr, comment="product")),
            events_json=json.dumps(events),
            wboR_json=json.dumps(wboR.tolist()),
            wboP_json=json.dumps(wboP.tolist()),
            elements_R_json=json.dumps(elR),
            elements_P_json=json.dumps(elP),
        )
        out_path = out_dir / f"seed_{i}.html"
        out_path.write_text(html)
        summaries.append({
            "i": i, "br": len(br), "fm": len(fm),
            "events": len(events), "mapped": len(mapping),
            "first": order[:5], "file": out_path.name,
        })

    # per-step index
    rows = ""
    for s in summaries:
        rows += (f"<tr><td>{s['i']}</td>"
                 f"<td>{s['br']}/{s['fm']}</td>"
                 f"<td>{s['events']}</td>"
                 f"<td>{s['mapped']}/{len(elR)}</td>"
                 f"<td>{s['first']}…</td>"
                 f"<td><a href='{s['file']}' target='_blank'>open</a></td></tr>")
    idx = f"""<!doctype html><html><head><meta charset="utf-8">
<title>{step} seed comparison</title>
<style>body{{font-family:-apple-system,sans-serif;margin:20px}}
table{{border-collapse:collapse}}
th,td{{border:1px solid #ccc;padding:6px 10px;font-size:13px}}</style>
</head><body>
<h2>{step}</h2>
<p>N atoms: {len(elR)} | charge: {chg} | uhf: {uhf} | algorithm: OLD (rxn_core_frag) | {N_SEEDS} random seeds</p>
<p><a href="../index.html">↑ all regressions</a></p>
<table><tr><th>seed#</th><th>br/fm</th><th>events</th><th>mapped</th><th>first 5 seeds</th><th>trace</th></tr>{rows}</table>
</body></html>"""
    (out_dir / "index.html").write_text(idx)
    return len(elR), summaries


def main():
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    top_rows = ""
    for step in REGRESSIONS:
        print(f"=== {step} ===")
        try:
            n, summaries = run_step(step)
        except Exception as e:
            print(f"  ERROR: {e}")
            top_rows += (f"<tr><td>{step}</td><td colspan='4'>ERROR: {e}</td></tr>")
            continue
        # per-step summary stats
        br_set = sorted({s['br'] for s in summaries})
        fm_set = sorted({s['fm'] for s in summaries})
        m_set = sorted({s['mapped'] for s in summaries})
        sanitized = re.sub(r"[^A-Za-z0-9._-]", "_", step)
        top_rows += (
            f"<tr><td><a href='{sanitized}/index.html'>{step}</a></td>"
            f"<td>{n}</td>"
            f"<td>{br_set}</td>"
            f"<td>{fm_set}</td>"
            f"<td>{m_set}</td></tr>")
        for s in summaries:
            print(f"  seed#{s['i']:>2}  br/fm={s['br']}/{s['fm']}  "
                  f"events={s['events']:>4}  mapped={s['mapped']:>3}/{n}")

    top = f"""<!doctype html><html><head><meta charset="utf-8">
<title>PQ regression diagnostic traces</title>
<style>body{{font-family:-apple-system,sans-serif;margin:20px;max-width:900px}}
table{{border-collapse:collapse}}
th,td{{border:1px solid #ccc;padding:6px 10px;font-size:13px}}</style>
</head><body>
<h2>PQ-vs-OLD regressions: 10-seed OLD-algorithm traces</h2>
<p>For each step, OLD algorithm (rxn_core_frag) was run with {N_SEEDS} random
seed orderings. Click a step to open per-seed trace HTMLs.</p>
<table>
<tr><th>step</th><th>N</th><th>broken (set across seeds)</th><th>formed (set)</th><th>mapped (set)</th></tr>
{top_rows}
</table>
</body></html>"""
    (OUT_ROOT / "index.html").write_text(top)
    print(f"\ntop index: {OUT_ROOT / 'index.html'}")


if __name__ == "__main__":
    main()
