"""
Generate per-seed-order trace HTMLs for a single tsdisco step, so the
user can compare how different seed orderings explore the alignment.

Usage:
  python build_seed_compare.py <step_id>     # e.g. pr1.tempo_ts8
"""
from __future__ import annotations
import json
import re
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from rxn_core_frag import (
    run_xtb, build_graph, write_xyz_str,
    classify_bonds,
)
from trace_run import find_islands_with_trace, HTML
from build_tsdisco_viewer import step_inputs


TSDISCO_INDEX = Path(
    "/Users/yunhengz/empty_for_claude/tsdisco_benchmark_visualization_plain_portable/index.html"
)
ROOT = Path("/Users/yunhengz/empty_for_claude/rxn_core")
OUT = ROOT / "out"
WORK = ROOT / "work_seedcompare"
OUT.mkdir(parents=True, exist_ok=True)
WORK.mkdir(parents=True, exist_ok=True)


def load_data():
    html = TSDISCO_INDEX.read_text()
    m = re.search(r"const TSDISCO_DATA = (\{.*?\});\s*\n", html, re.DOTALL)
    return json.loads(m.group(1))


def main():
    step_id = sys.argv[1] if len(sys.argv) > 1 else "pr1.tempo_ts8"
    n_seeds = int(sys.argv[2]) if len(sys.argv) > 2 else 10
    rng_seed = int(sys.argv[3]) if len(sys.argv) > 3 else 42

    data = load_data()
    step = next((s for s in data["steps"] if s["step_id"] == step_id), None)
    if step is None:
        print(f"step '{step_id}' not found")
        sys.exit(1)

    sanitized = re.sub(r"[^A-Za-z0-9._-]", "_", step_id)
    wd = WORK / sanitized
    wd.mkdir(parents=True, exist_ok=True)
    chg = step.get("charge", 0) or 0
    uhf = max(0, (step.get("multiplicity", 1) or 1) - 1)

    rxyz, pxyz, _, _ = step_inputs(step)
    (wd / "reactant.xyz").write_text(rxyz)
    (wd / "product.xyz").write_text(pxyz)

    elR, xyzR_arr, wboR = run_xtb(wd / "reactant.xyz", wd / "R", charge=chg, uhf=uhf)
    elP, xyzP_arr, wboP = run_xtb(wd / "product.xyz", wd / "P", charge=chg, uhf=uhf)
    g_R = build_graph(elR, wboR)
    g_P = build_graph(elP, wboP)

    nodes = list(g_R.nodes())
    rng = random.Random(rng_seed)
    seed_orders = []
    for i in range(n_seeds):
        perm = list(nodes); rng.shuffle(perm); seed_orders.append(perm)

    print(f"Running {n_seeds} random seedings on {step_id}  (N={len(elR)})")
    summaries = []
    for i, order in enumerate(seed_orders):
        mapping, events = find_islands_with_trace(g_R, g_P, seed_order=order)
        br, fm, _, _ = classify_bonds(mapping, wboR, wboP)
        title = f"{step_id}  seed#{i}  br/fm={len(br)}/{len(fm)}  events={len(events)}  mapped={len(mapping)}"
        xyzR_str = write_xyz_str(elR, xyzR_arr, comment="reactant")
        xyzP_str = write_xyz_str(elP, xyzP_arr, comment="product")
        html = HTML.format(
            title=title,
            xyzR_json=json.dumps(xyzR_str),
            xyzP_json=json.dumps(xyzP_str),
            events_json=json.dumps(events),
            wboR_json=json.dumps(wboR.tolist()),
            wboP_json=json.dumps(wboP.tolist()),
            elements_R_json=json.dumps(elR),
            elements_P_json=json.dumps(elP),
        )
        out_path = OUT / f"seed_{sanitized}_{i}.html"
        out_path.write_text(html)
        summaries.append((i, len(br), len(fm), len(events), len(mapping), order[:5], str(out_path)))
        print(f"  seed#{i:>2}  br/fm={len(br)}/{len(fm)}  "
              f"events={len(events):>4}  mapped={len(mapping):>3}  "
              f"first-seeds={order[:5]}  -> {out_path.name}")

    # Index page linking to all seeds
    rows = ""
    for i, br, fm, ev, m, head, path in summaries:
        rows += f"<tr><td>{i}</td><td>{br}/{fm}</td><td>{ev}</td><td>{m}</td>"
        rows += f"<td>{head}…</td><td><a href='{Path(path).name}' target='_blank'>open</a></td></tr>"
    index = f"""<!doctype html><html><head><meta charset="utf-8"><title>{step_id} seed comparison</title>
<style>body{{font-family:-apple-system,sans-serif;margin:20px}} table{{border-collapse:collapse}} th,td{{border:1px solid #ccc;padding:6px 10px;font-size:13px}}</style>
</head><body>
<h2>{step_id} — {n_seeds} random seed orderings</h2>
<p>N atoms: {len(elR)}; charge={chg}; uhf={uhf}</p>
<table><tr><th>seed#</th><th>br/fm</th><th>events</th><th>mapped</th><th>first 5 seeds</th><th>trace</th></tr>{rows}</table>
</body></html>"""
    index_path = OUT / f"seed_{sanitized}_index.html"
    index_path.write_text(index)
    print(f"\nindex: {index_path}")


if __name__ == "__main__":
    main()
