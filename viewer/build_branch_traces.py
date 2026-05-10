"""
Render slider-driven HTML alignment traces for every UNIQUE alignment
branch of one IG (e.g. ru_nh3_1a iter20).

The PQ aligner returns a per-seed branches list; calling find_islands_pq
with events= records ONLY branch 0's trajectory per seed (per the design
in src/rxn_core_pq.py). So to get traces for distinct branches we run
many random seed orderings and look at branches[0] from each — different
seeds tend to land on different "best" mappings.

Pipeline:
  1. Load R and IG xtb output from out/ranked_views/<step>/xtb/.
  2. Run alignment with return_all=True to enumerate the K unique
     mappings we want traces for.
  3. Try up to N_SEEDS random seed orderings; for each, record
     (resulting branches[0] mapping, events trace, br/fm counts).
  4. Group by mapping equivalence; keep one representative seed
     per unique mapping. Stop when we've hit K.
  5. Render one trace HTML per unique branch + an index.

Output:
  out/ranked_views/<step>/branches/<ig_label>/
    branch_<i>.html
    index.html

Usage:
  python viewer/build_branch_traces.py <step> <ig_label> [n_seeds=200]
  e.g.  python viewer/build_branch_traces.py ru_nh3_1a iter20
"""
from __future__ import annotations
import sys as _sys
from pathlib import Path as _Path
_HERE = _Path(__file__).resolve().parent
_sys.path.insert(0, str(_HERE.parent / "src"))
_sys.path.insert(0, str(_HERE))
PROJECT_ROOT = _HERE.parent

import json
import random
import re
import sys
from pathlib import Path

import numpy as np

from rxn_core_pq import find_islands_pq, align_from_arrays, build_graph
from rxn_core_frag import (
    parse_xyz, write_xyz_str, classify_bonds, expand_mapping,
)
from analyze_core_modes import (
    parse_g98_modes, core_atoms_in_R_frame, reindex_modes_to_R,
    bond_reaction_vector, bond_overlap_per_mode,
    rxn_overlap_per_mode, reaction_coord_delta,
)
from align_bgcp_coords import load_cached_xtb, fill_unmapped_greedy
from trace_html import HTML


W_RXN, W_CORE, IMAG_PEN = 1.0, 0.2, 0.3
RNG_SEED = 42


def patch_html_for_pq(html_str):
    extra = """} else if (lastEvent.type === 'consumed') {
    txt += `  consumed edge: R[${lastEvent.frag_atom}] → R[${lastEvent.ext_atom}]  WBO=${lastEvent.wbo}  reason=${lastEvent.reason}`;
"""
    return html_str.replace(
        "}} else if (lastEvent.type === 'pass_start') {{",
        extra.replace("{", "{{").replace("}", "}}") +
        "\n  }} else if (lastEvent.type === 'pass_start') {{"
    )


def load_xtb_for_hess(hess_dir: Path):
    """Hess dirs hold both the IG input xyz and xtbhess.xyz; pick the input."""
    ig_xyz = next(p for p in hess_dir.glob("*.xyz")
                  if "xtbhess" not in p.name)
    el, xyz = parse_xyz(ig_xyz)
    n = len(el)
    wbo = np.zeros((n, n))
    for ln in (hess_dir / "wbo").read_text().splitlines():
        parts = ln.split()
        if len(parts) < 3: continue
        i, j = int(parts[0]) - 1, int(parts[1]) - 1
        v = float(parts[2])
        wbo[i, j] = v; wbo[j, i] = v
    return el, np.asarray(xyz, float), wbo


def score_mapping(elR, xyzR, xyzT, mapping_RT, modes_TS, freqs,
                  broken_R, formed_R, core_R, delta_RP):
    """Compute (score, picked_k, b, r, c) for a given IG<->R mapping."""
    modes_R = reindex_modes_to_R(modes_TS, mapping_RT, len(elR))
    sq = (modes_R ** 2).sum(axis=2)
    total = sq.sum(axis=1)
    core_e = (sq[:, core_R].sum(axis=1) if core_R
              else np.zeros(modes_R.shape[0]))
    kappa = np.where(total > 1e-12, core_e / total, 0.0)
    rho = rxn_overlap_per_mode(modes_R, delta_RP, core_R)
    ts_xyz_in_R = np.zeros_like(np.asarray(xyzR))
    for r, t in mapping_RT.items():
        ts_xyz_in_R[r] = xyzT[t]
    V = bond_reaction_vector(ts_xyz_in_R, broken_R, formed_R)
    beta = bond_overlap_per_mode(modes_R, V)
    imag_idx = list(np.where(freqs < 0)[0])
    if not imag_idx:
        return 0.0, None, 0.0, 0.0, 0.0
    n_imag = len(imag_idx)
    picked_k = max(imag_idx, key=lambda k: beta[k])
    b = float(beta[picked_k])
    r_ = float(rho[picked_k])
    c = float(kappa[picked_k])
    score = (b * (1 + W_RXN * r_) * (1 + W_CORE * c)
             / max(n_imag, 1) ** IMAG_PEN)
    return score, picked_k, b, r_, c


def main():
    step = sys.argv[1]
    label = sys.argv[2]
    n_seeds = int(sys.argv[3]) if len(sys.argv) > 3 else 200

    cache = PROJECT_ROOT / "out" / "ranked_views" / step / "xtb"
    elR, xyzR, wboR, _ = load_cached_xtb(cache / "R")
    elP, xyzP, wboP, _ = load_cached_xtb(cache / "P")
    elT, xyzT, wboT = load_xtb_for_hess(cache / f"hess_{label}")
    freqs, modes_TS = parse_g98_modes(cache / f"hess_{label}" / "g98.out")

    # Step-level info
    rp = align_from_arrays(elR, xyzR, wboR, elP, xyzP, wboP)
    mapping_RP = dict(rp["mapping"])
    inv_RP = {v: k for k, v in mapping_RP.items()}
    core_R = core_atoms_in_R_frame(mapping_RP, rp["broken"], rp["formed"])
    full_RP = fill_unmapped_greedy(elR, xyzR, elP, xyzP, mapping_RP)
    delta_RP = reaction_coord_delta(np.asarray(xyzR, float),
                                     np.asarray(xyzP, float), full_RP)
    broken_R = [(int(a), int(b)) for (a, b, _, _) in rp["broken"]]
    formed_R = [(int(inv_RP[a]), int(inv_RP[b]))
                for (a, b, _, _) in rp["formed"]
                if a in inv_RP and b in inv_RP]

    # Enumerate the unique branches we want to trace
    it = align_from_arrays(elR, xyzR, wboR, elT, xyzT, wboT, return_all=True)
    raw_branches = it.get("all_scored", [])
    target_keys = []
    seen_keys = set()
    for (_, mapping, _, _, _) in raw_branches:
        key = tuple(sorted(dict(mapping).items()))
        if key not in seen_keys:
            seen_keys.add(key)
            target_keys.append(key)
    print(f"unique target branches: {len(target_keys)}")

    # Spin random seed orderings; for each, capture branches[0] mapping
    # + its event trace, group by mapping key
    g_R = build_graph(elR, wboR, bond_cut=0.2)
    g_T = build_graph(elT, wboT, bond_cut=0.2)
    nodes = list(g_R.nodes())
    rng = random.Random(RNG_SEED)
    captured = {}  # expanded mapping_key -> (events, branch, seed_idx, order)
    for seed_i in range(n_seeds):
        order = list(nodes); rng.shuffle(order)
        events = []
        try:
            branches = find_islands_pq(g_R, g_T, order, events=events)
        except Exception as e:
            continue
        if not branches:
            continue
        b0 = branches[0]
        # Match the same expand_mapping call that align_from_arrays uses
        # internally, so the key here lines up with target_keys.
        expanded = expand_mapping(dict(b0.mapping), g_R, g_T)
        mkey = tuple(sorted(expanded.items()))
        if mkey in captured:
            continue
        captured[mkey] = (events, b0, seed_i, order)
        if len(captured) >= len(target_keys):
            break

    print(f"covered {len(captured)} / {len(target_keys)} target branches "
          f"after {seed_i+1} seeds")

    # Render
    out_dir = PROJECT_ROOT / "out" / "ranked_views" / step / "branches" / label
    out_dir.mkdir(parents=True, exist_ok=True)

    # First pass: compute scores for every captured mapping. Then sort
    # by score descending and assign sequential branch IDs 1..N. This
    # avoids the problem that align_from_arrays' 10-seed search and
    # this tracer's wider seed search cover overlapping-but-not-equal
    # subsets of the actually-distinct mapping space.
    scored = []
    for mkey, (events, b0, seed_i, order) in captured.items():
        mapping_RT = fill_unmapped_greedy(elR, xyzR, elT, xyzT, dict(b0.mapping))
        score, picked_k, b, r, c = score_mapping(
            elR, xyzR, xyzT, mapping_RT, modes_TS, freqs,
            broken_R, formed_R, core_R, delta_RP,
        )
        scored.append({
            "mkey": mkey, "events": events, "b0": b0,
            "seed_i": seed_i, "score": score, "picked_k": picked_k,
            "beta": b, "rho": r, "kappa": c,
        })
    scored.sort(key=lambda r: -r["score"])

    summaries = []
    for target_id, item in enumerate(scored, 1):
        mkey = item["mkey"]; events = item["events"]; b0 = item["b0"]
        seed_i = item["seed_i"]; score = item["score"]
        b = item["beta"]; r = item["rho"]; c = item["kappa"]
        mapping_RT = fill_unmapped_greedy(elR, xyzR, elT, xyzT, dict(b0.mapping))
        mapping_full = expand_mapping(mapping_RT, g_R, g_T)
        br, fm, _, _ = classify_bonds(mapping_full, wboR, wboT)
        score, picked_k, b, r, c = score_mapping(
            elR, xyzR, xyzT, mapping_RT, modes_TS, freqs,
            broken_R, formed_R, core_R, delta_RP,
        )
        n_imag = int(np.sum(freqs < 0))
        title = (f"{step}  {label}  branch #{target_id}  "
                 f"mapped={len(mapping_full)}/{len(elR)}  "
                 f"events={len(events)}  S={score:.4f}  "
                 f"β={b:.3f} ρ={r:.3f} κ={c:.3f} n_imag={n_imag}  "
                 f"(seed#{seed_i})")
        html = HTML.format(
            title=title,
            xyzR_json=json.dumps(write_xyz_str(elR, xyzR, comment="R")),
            xyzP_json=json.dumps(write_xyz_str(elT, xyzT, comment=label)),
            events_json=json.dumps(events),
            wboR_json=json.dumps(wboR.tolist()),
            wboP_json=json.dumps(wboT.tolist()),
            elements_R_json=json.dumps(elR),
            elements_P_json=json.dumps(elT),
        )
        html = patch_html_for_pq(html)
        out_path = out_dir / f"branch_{target_id}.html"
        out_path.write_text(html)
        summaries.append({
            "id": target_id, "score": score, "beta": b, "rho": r,
            "kappa": c, "n_imag": n_imag, "br": len(br), "fm": len(fm),
            "mapped": len(mapping_full), "events": len(events),
            "seed_i": seed_i, "file": out_path.name,
        })
        print(f"  branch #{target_id}: S={score:.4f}  β={b:.3f}  "
              f"ρ={r:.3f}  κ={c:.3f}  events={len(events)}  "
              f"-> {out_path.name}")

    summaries.sort(key=lambda s: -s["score"])
    rows = ""
    for s in summaries:
        rows += (f"<tr><td>{s['id']}</td>"
                 f"<td>{s['score']:.4f}</td>"
                 f"<td>{s['beta']:.3f}</td>"
                 f"<td>{s['rho']:.3f}</td>"
                 f"<td>{s['kappa']:.3f}</td>"
                 f"<td>{s['n_imag']}</td>"
                 f"<td>{s['br']}/{s['fm']}</td>"
                 f"<td>{s['mapped']}/{len(elR)}</td>"
                 f"<td>{s['events']}</td>"
                 f"<td>{s['seed_i']}</td>"
                 f"<td><a href='{s['file']}' target='_blank'>open</a></td></tr>")
    idx = f"""<!doctype html><html><head><meta charset="utf-8">
<title>{step} / {label} branch traces</title>
<style>body{{font-family:-apple-system,sans-serif;margin:20px;max-width:1200px}}
table{{border-collapse:collapse}}
th,td{{border:1px solid #ccc;padding:6px 10px;font-size:13px}}
caption{{caption-side:top;text-align:left;font-size:14px;padding:6px 0;font-weight:600}}</style>
</head><body>
<h2>{step} — {label} alignment branches</h2>
<p>Each branch is one unique IG&hArr;R atom mapping returned by the
priority-queue aligner. Click a row to open the slider-driven trace
of how that mapping was built. The ranker score column applies the
verifier's S formula to the picked imaginary mode for each mapping;
build_ranked_view_external.py picks the row with the highest S as
the displayed result.</p>
<table>
<tr><th>branch&nbsp;#</th><th>S</th><th>β</th><th>ρ</th><th>κ</th>
<th>n_imag</th><th>br/fm</th><th>mapped</th><th>events</th>
<th>seed&nbsp;#</th><th>trace</th></tr>
{rows}
</table>
</body></html>"""
    (out_dir / "index.html").write_text(idx)
    print(f"\nindex: {out_dir / 'index.html'}")


if __name__ == "__main__":
    main()
