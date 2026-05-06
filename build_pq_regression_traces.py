"""
Generate NEW-PQ-algorithm trace HTMLs for each of the 7 regressed steps,
10 random seed orderings each. Output: out/regressions/<step>/pq_seed_*.html
plus a per-step pq_index.html linking to all 10 seeds.

Reuses trace_html.HTML — events emitted by grow_island_pq use the same
schema (seed_start, commit, seed_end, island_locked, pass_start, done)
plus an extra 'consumed' event type (silently skipped by the renderer
unless we patch it; it still increments the slider so the user can step
through every algorithm decision).

Usage: python build_pq_regression_traces.py
"""
from __future__ import annotations
import json
import random
import re
from pathlib import Path

from rxn_core_pq import find_islands_pq
from rxn_core_frag import (
    run_xtb, build_graph, write_xyz_str, classify_bonds, expand_mapping,
)
from trace_html import HTML
from bgcp_io import BGCP_ROOT, LOOKUP, WORK as BGCP_WORK, read_xyzs


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


def patch_html_for_pq(html_str):
    """Inject handling for 'consumed' events into the trace HTML's
    event-log JS so they render meaningfully."""
    extra = """} else if (lastEvent.type === 'consumed') {
    txt += `  consumed edge: R[${lastEvent.frag_atom}] → R[${lastEvent.ext_atom}]  WBO=${lastEvent.wbo}  reason=${lastEvent.reason}`;
"""
    # Inject right before the pass_start clause
    return html_str.replace(
        "}} else if (lastEvent.type === 'pass_start') {{",
        extra.replace("{", "{{").replace("}", "}}") +
        "\n  }} else if (lastEvent.type === 'pass_start') {{"
    )


def run_step(step):
    sanitized = re.sub(r"[^A-Za-z0-9._-]", "_", step)
    out_dir = OUT_ROOT / sanitized
    out_dir.mkdir(parents=True, exist_ok=True)

    # Use BGCP data + LOOKUP charge/multiplicity (not the old Benchmark fallback)
    # so the WBO matches what analyze_pq sees.
    chg, uhf = LOOKUP.get(step, (0, 0))
    work = BGCP_WORK / sanitized
    work.mkdir(parents=True, exist_ok=True)
    rxyz = read_xyzs(BGCP_ROOT / step / "reactants")
    pxyz = read_xyzs(BGCP_ROOT / step / "products")
    if rxyz is None or pxyz is None:
        raise RuntimeError(f"missing R or P for {step}")
    rxyz_path = work / "reactant.xyz"
    pxyz_path = work / "product.xyz"
    rxyz_path.write_text(rxyz)
    pxyz_path.write_text(pxyz)

    elR, xyzR_arr, wboR = run_xtb(rxyz_path, work / "R", charge=chg, uhf=uhf)
    elP, xyzP_arr, wboP = run_xtb(pxyz_path, work / "P", charge=chg, uhf=uhf)
    g_R = build_graph(elR, wboR, bond_cut=0.2)
    g_P = build_graph(elP, wboP, bond_cut=0.2)

    nodes = list(g_R.nodes())
    rng = random.Random(RNG_SEED)
    seed_orders = []
    for _ in range(N_SEEDS):
        perm = list(nodes); rng.shuffle(perm); seed_orders.append(perm)

    summaries = []
    for i, order in enumerate(seed_orders):
        events = []
        branches = find_islands_pq(g_R, g_P, order, events=events)
        if not branches:
            continue
        b = branches[0]
        mapping = expand_mapping(dict(b.mapping), g_R, g_P)
        br, fm, _, _ = classify_bonds(mapping, wboR, wboP)
        title = (f"PQ {step}  seed#{i}  br/fm={len(br)}/{len(fm)}  "
                 f"events={len(events)}  mapped={len(mapping)}/{len(elR)}  "
                 f"branches={len(branches)}")
        html = HTML.format(
            title=title,
            xyzR_json=json.dumps(write_xyz_str(elR, xyzR_arr, comment="R")),
            xyzP_json=json.dumps(write_xyz_str(elP, xyzP_arr, comment="P")),
            events_json=json.dumps(events),
            wboR_json=json.dumps(wboR.tolist()),
            wboP_json=json.dumps(wboP.tolist()),
            elements_R_json=json.dumps(elR),
            elements_P_json=json.dumps(elP),
        )
        html = patch_html_for_pq(html)
        out_path = out_dir / f"pq_seed_{i}.html"
        out_path.write_text(html)
        summaries.append({
            "i": i, "br": len(br), "fm": len(fm),
            "events": len(events), "mapped": len(mapping),
            "branches": len(branches), "first": order[:5],
            "file": out_path.name,
        })

    rows = ""
    for s in summaries:
        rows += (f"<tr><td>{s['i']}</td>"
                 f"<td>{s['br']}/{s['fm']}</td>"
                 f"<td>{s['events']}</td>"
                 f"<td>{s['mapped']}/{len(elR)}</td>"
                 f"<td>{s['branches']}</td>"
                 f"<td>{s['first']}…</td>"
                 f"<td><a href='{s['file']}' target='_blank'>open</a></td></tr>")
    idx = f"""<!doctype html><html><head><meta charset="utf-8">
<title>PQ {step} traces</title>
<style>body{{font-family:-apple-system,sans-serif;margin:20px}}
table{{border-collapse:collapse}}
th,td{{border:1px solid #ccc;padding:6px 10px;font-size:13px}}</style>
</head><body>
<h2>{step} — NEW PQ algorithm, {N_SEEDS} random seeds</h2>
<p>N atoms: {len(elR)} | charge: {chg} | uhf: {uhf}</p>
<p><a href="index.html">↻ OLD 10-seed traces</a> |
<a href="pq_result.html">→ PQ static result</a> |
<a href="../index.html">↑↑ all regressions</a></p>
<table><tr><th>seed#</th><th>br/fm</th><th>events</th>
<th>mapped</th><th>#branches</th><th>first 5 seeds</th><th>trace</th></tr>{rows}</table>
</body></html>"""
    (out_dir / "pq_index.html").write_text(idx)
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
            import traceback; traceback.print_exc()
            continue
        sanitized = re.sub(r"[^A-Za-z0-9._-]", "_", step)
        br_set = sorted({s['br'] for s in summaries})
        fm_set = sorted({s['fm'] for s in summaries})
        m_set = sorted({s['mapped'] for s in summaries})
        for s in summaries:
            print(f"  seed#{s['i']:>2}  br/fm={s['br']}/{s['fm']}  "
                  f"events={s['events']:>4}  mapped={s['mapped']:>3}/{n}  "
                  f"branches={s['branches']}")
        top_rows += (
            f"<tr><td>{step}</td>"
            f"<td>{n}</td>"
            f"<td>{br_set}</td><td>{fm_set}</td><td>{m_set}</td>"
            f"<td><a href='{sanitized}/pq_index.html'>PQ traces</a></td>"
            f"<td><a href='{sanitized}/pq_result.html'>PQ result</a></td>"
            f"<td><a href='{sanitized}/index.html'>OLD traces</a></td></tr>"
        )

    top = f"""<!doctype html><html><head><meta charset="utf-8">
<title>PQ regressions diagnostic</title>
<style>body{{font-family:-apple-system,sans-serif;margin:20px;max-width:1200px}}
table{{border-collapse:collapse}}
th,td{{border:1px solid #ccc;padding:6px 10px;font-size:13px}}</style>
</head><body>
<h2>PQ-vs-OLD regressions: side-by-side traces</h2>
<table>
<tr><th>step</th><th>N</th><th>PQ broken (set)</th><th>PQ formed (set)</th>
<th>PQ mapped (set)</th>
<th>NEW PQ traces</th><th>NEW PQ result</th><th>OLD traces</th></tr>
{top_rows}
</table>
</body></html>"""
    (OUT_ROOT / "index.html").write_text(top)
    print(f"\ntop index: {OUT_ROOT / 'index.html'}")


if __name__ == "__main__":
    main()
