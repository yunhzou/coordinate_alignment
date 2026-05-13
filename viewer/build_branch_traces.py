"""
Render ONE slider-driven HTML trace for an alignment (R <-> GT or
R <-> IG), with all unique mappings tabulated in a panel below.

The WBO aligner returns a per-seed branches list; find_islands with
events= records the trajectory of branches[0] per seed. We spin many
random seed orderings to discover all unique mappings the algorithm can
land on, score them all, and render a single HTML containing:

  - Slider trace from the highest-S seed (the one the verifier picks)
  - Static table of every unique mapping found (S, beta, rho, kappa,
    br/fm count, witness mapping R->T)
  - In the event log, seed_end events now list ALL final candidate
    mappings (via the all_isos field in seed_end events)

Cache layouts supported:
  - BGCP:    appendix_perparation/xtb_frequency_calculations/<step>/
             {R, P, sp_<label>, hess_<label>}/
  - Legacy:  out/ranked_views/<step>/xtb/{R, P, hess_<label>}/

Output: out/bgcp_views/<step>/branches/<label>.html (BGCP)
        out/ranked_views/<step>/branches/<label>.html (legacy)

Usage:
  python viewer/build_branch_traces.py <step> <label> [n_seeds=200]
  e.g.  python viewer/build_branch_traces.py pr7.V.dodh_ts910 groundtruth
        python viewer/build_branch_traces.py pr7.V.dodh_ts910 iter11
"""
from __future__ import annotations
import sys as _sys
from pathlib import Path as _Path
_HERE = _Path(__file__).resolve().parent
_sys.path.insert(0, str(_HERE.parent / "src"))
PROJECT_ROOT = _HERE.parent

import json
import random
import re
import sys
from pathlib import Path

import numpy as np

from rxn_core import (
    align_from_arrays, find_islands, build_graph,
    bond_overlap_per_mode, bond_reaction_vector,
    classify_bonds, core_atoms_in_R_frame, expand_mapping,
    load_cached_xtb,
    parse_g98_modes, parse_xyz, reaction_coord_delta,
    reindex_modes_to_R, rxn_overlap_per_mode, write_xyz_str,
)
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
    mode_norms = np.linalg.norm(modes_TS.reshape(modes_TS.shape[0], -1), axis=1)
    sq = (modes_R ** 2).sum(axis=2)
    total = mode_norms ** 2
    core_e = (sq[:, core_R].sum(axis=1) if core_R
              else np.zeros(modes_R.shape[0]))
    kappa = np.where(total > 1e-12, core_e / total, 0.0)
    rho = rxn_overlap_per_mode(modes_R, delta_RP, core_R,
                                mode_norms=mode_norms)
    ts_xyz_in_R = np.asarray(xyzR, float).copy()
    for r, t in mapping_RT.items():
        ts_xyz_in_R[r] = xyzT[t]
    V = bond_reaction_vector(ts_xyz_in_R, broken_R, formed_R)
    beta = bond_overlap_per_mode(modes_R, V, mode_norms=mode_norms)
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

    # BGCP cache layout. Legacy ranked_views layout falls through.
    bgcp_cache = PROJECT_ROOT / "appendix_perparation" / "xtb_frequency_calculations" / step
    legacy_cache = PROJECT_ROOT / "out" / "ranked_views" / step / "xtb"
    if bgcp_cache.exists():
        cache = bgcp_cache
        sp_dir = cache / ("sp_groundtruth" if label == "groundtruth" else f"sp_{label}")
        hess_dir = cache / ("hess_groundtruth" if label == "groundtruth" else f"hess_{label}")
    else:
        cache = legacy_cache
        sp_dir = cache / f"hess_{label}"  # legacy: xyz lives next to hess
        hess_dir = cache / f"hess_{label}"
    elR, xyzR, wboR, _ = load_cached_xtb(cache / "R")
    elP, xyzP, wboP, _ = load_cached_xtb(cache / "P")
    elT, xyzT, wboT = load_xtb_for_hess(sp_dir)
    freqs, modes_TS = parse_g98_modes(hess_dir / "g98.out")

    # Enumerate ALL unique R<->P branches, then DEDUP by chemistry class
    # (broken, formed, core). Within a class, branches differ only in
    # spectator-atom permutations that produce 3rd-decimal noise on S.
    rp_all = align_from_arrays(elR, xyzR, wboR, elP, xyzP, wboP,
                               return_all=True, max_branches=1_000_000)
    n_rp_raw = 0
    seen_witness = set()
    seen_chem = {}  # (broken, formed, core) -> rp_branch dict
    for (_, mapping, broken, formed, _) in rp_all.get("all_scored", []):
        witness_key = tuple(sorted(dict(mapping).items()))
        if witness_key in seen_witness: continue
        seen_witness.add(witness_key)
        n_rp_raw += 1
        mapping_RP = dict(mapping)
        inv_RP = {v: k for k, v in mapping_RP.items()}
        delta_RP = reaction_coord_delta(np.asarray(xyzR, float),
                                         np.asarray(xyzP, float), mapping_RP)
        broken_R = [(int(a), int(b)) for (a, b, _, _) in broken]
        formed_R = [(int(inv_RP[a]), int(inv_RP[b]))
                    for (a, b, _, _) in formed if a in inv_RP and b in inv_RP]
        core_R = list(core_atoms_in_R_frame(mapping_RP, broken, formed))
        chem_key = (
            tuple(sorted((min(a, b), max(a, b)) for (a, b) in broken_R)),
            tuple(sorted((min(a, b), max(a, b)) for (a, b) in formed_R)),
            tuple(sorted(core_R)),
        )
        if chem_key in seen_chem:
            seen_chem[chem_key]["class_size"] += 1
            continue
        seen_chem[chem_key] = {
            "mapping_RP": mapping_RP, "core_R": core_R,
            "broken_R": broken_R, "formed_R": formed_R, "delta_RP": delta_RP,
            "n_broken": len(broken), "n_formed": len(formed),
            "class_size": 1,
        }
    rp_branches = list(seen_chem.values())
    print(f"R<->P: {n_rp_raw} raw branches -> {len(rp_branches)} chemistry classes")

    # Enumerate ALL unique R<->GT branches
    it = align_from_arrays(elR, xyzR, wboR, elT, xyzT, wboT,
                           return_all=True, max_branches=1_000_000)
    raw_branches = it.get("all_scored", [])
    target_keys = []
    seen_keys = set()
    for (_, mapping, _, _, _) in raw_branches:
        key = tuple(sorted(dict(mapping).items()))
        if key not in seen_keys:
            seen_keys.add(key)
            target_keys.append(key)
    print(f"unique R<->GT branches: {len(target_keys)}")

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
            branches = find_islands(g_R, g_T, order, events=events)
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

    # Output: ONE combined HTML, not N branch files. Slider trace from a
    # representative seed + a static "all unique mappings" table below.
    out_root = (PROJECT_ROOT / "out" / "bgcp_views"
                if bgcp_cache.exists()
                else PROJECT_ROOT / "out" / "ranked_views")
    out_dir = out_root / step / "branches"
    out_dir.mkdir(parents=True, exist_ok=True)

    # Cross-product scoring: every (R<->P branch) x (R<->GT branch).
    # The R<->P branch determines the broken/formed bonds, core_R, and
    # delta_RP that go into the verifier score. Pick the highest-S pair.
    captured_by_target_key = {}
    for mkey, (events, b0, seed_i, order) in captured.items():
        captured_by_target_key[mkey] = (events, b0, seed_i)

    scored = []
    for rp_idx, rp_b in enumerate(rp_branches, 1):
        for tkey in target_keys:
            target_mapping = dict(tkey)
            mapping_RT = dict(target_mapping)
            mapping_full = expand_mapping(dict(mapping_RT), g_R, g_T)
            witness_key = tuple(sorted(mapping_full.items()))
            cap = captured_by_target_key.get(witness_key) or captured_by_target_key.get(tkey)
            events_t = cap[0] if cap else None
            seed_i = cap[2] if cap else None
            score, picked_k, b, r, c = score_mapping(
                elR, xyzR, xyzT, mapping_RT, modes_TS, freqs,
                rp_b["broken_R"], rp_b["formed_R"], rp_b["core_R"], rp_b["delta_RP"],
            )
            br, fm, _, _ = classify_bonds(mapping_full, wboR, wboT)
            scored.append({
                "rp_idx": rp_idx,
                "rp_branch": rp_b,
                "events": events_t, "seed_i": seed_i,
                "mapping_RT": mapping_RT, "mapping_full": mapping_full,
                "tkey": tkey,
                "score": score, "picked_k": picked_k,
                "beta": b, "rho": r, "kappa": c,
                "n_imag": int(np.sum(freqs < 0)),
                "n_broken_rp": rp_b["n_broken"], "n_formed_rp": rp_b["n_formed"],
                "n_broken": len(br), "n_formed": len(fm),
                "n_mapped": len(mapping_full),
                "has_trace": cap is not None,
            })
    scored.sort(key=lambda r: -r["score"])
    print(f"scored {len(scored)} (R<->P) x (R<->GT) combinations")

    # Slider trace: use the highest-S R<->GT mapping that was captured by
    # a random seed. (R<->P branch only affects scoring, not the R<->GT
    # event log.) Pick the first captured one in score order.
    best = scored[0]
    captured_scored = [s for s in scored if s["has_trace"]]
    slider_choice = captured_scored[0] if captured_scored else best
    events = slider_choice["events"] or []
    print(f"  best combo: R<->P #{best['rp_idx']} + R<->GT  -->  "
          f"S={best['score']:.4f}  β={best['beta']:.3f}  "
          f"ρ={best['rho']:.3f}  κ={best['kappa']:.3f}")
    print(f"  best combo R27 -> T{best['mapping_RT'].get(27)}    "
          f"R22 -> T{best['mapping_RT'].get(22)}")
    print(f"  unique R<->GT keys with trace: "
          f"{len(set(s['tkey'] for s in captured_scored))}")
    # Show top 20 in console
    print(f"\n  Top 20 by S (out of {len(scored)} combos):")
    for i, s in enumerate(scored[:20], 1):
        ev_n = len(s['events']) if s['events'] else 0
        tag = "T" if s["has_trace"] else "-"
        print(f"   #{i:>3} RP{s['rp_idx']:<2} [{tag}]  S={s['score']:.4f}  "
              f"β={s['beta']:.3f}  ρ={s['rho']:.3f}  κ={s['kappa']:.3f}  "
              f"br/fm(RP)={s['n_broken_rp']}/{s['n_formed_rp']}  "
              f"R27->T{s['mapping_RT'].get(27)}  R22->T{s['mapping_RT'].get(22)}  ev={ev_n}")

    # Build the static table. Limit to top-40 rows to keep the HTML small;
    # the rest is summarized.
    MAX_ROWS = 60

    def fmt_mapping_inline(m_full):
        items = sorted(m_full.items())
        return ", ".join(f"R[{r}]({elR[r]})&rarr;T[{t}]({elT[t]})"
                          for r, t in items)
    mappings_rows = ""
    for i, s in enumerate(scored[:MAX_ROWS], 1):
        is_slider = (s is slider_choice)
        markers = []
        if i == 1: markers.append("&#9733;")
        if is_slider: markers.append("&#9654; trace")
        marker = " ".join(markers)
        ev_n = len(s['events']) if s['events'] else 0
        seed_str = str(s['seed_i']) if s['seed_i'] is not None else '-'
        trace_status = "yes" if s["has_trace"] else "no"
        r27 = s["mapping_RT"].get(27)
        r22 = s["mapping_RT"].get(22)
        mappings_rows += (
            f"<tr class='{'picked' if i == 1 else ''}'>"
            f"<td>{marker} #{i}</td>"
            f"<td>RP#{s['rp_idx']}</td>"
            f"<td><b>{s['score']:.4f}</b></td>"
            f"<td>{s['beta']:.3f}</td>"
            f"<td>{s['rho']:.3f}</td>"
            f"<td>{s['kappa']:.3f}</td>"
            f"<td>{s['n_broken_rp']}/{s['n_formed_rp']}</td>"
            f"<td>R27&rarr;T{r27}({elT[r27] if r27 is not None else '?'})</td>"
            f"<td>R22&rarr;T{r22}({elT[r22] if r22 is not None else '?'})</td>"
            f"<td>{trace_status}</td>"
            f"<td>{ev_n}</td>"
            f"<td style='font-family:ui-monospace,monospace; font-size:11px; "
            f"max-width:600px; word-break:break-all'>"
            f"{fmt_mapping_inline(s['mapping_full'])}</td></tr>"
        )
    truncation_note = (
        f"<div class='note'>showing top {MAX_ROWS} of "
        f"{len(scored)} combos.</div>"
        if len(scored) > MAX_ROWS else ""
    )
    # Friendly side titles. Defined here so we can reference right_short
    # in the summary text.
    if label == "groundtruth":
        right_title = "Ground-truth TS"
        right_short = "GT"
    elif label.startswith("iter"):
        right_title = f"IG {label} TS"
        right_short = label
    else:
        right_title = label
        right_short = label[:3].upper()
    n_with_trace = sum(1 for s in scored if s["has_trace"])
    mappings_summary_html = (
        "<style>"
        "table.mt {border-collapse:collapse; width:100%}"
        "table.mt th, table.mt td {border:1px solid #ccc; padding:3px 6px; "
        "font-size:11px; text-align:left}"
        "table.mt th {background:#eef}"
        "table.mt tr.picked td {background:#fff7d6}"
        ".note {color:#666; font-size:12px; margin:4px 0}"
        "</style>"
        f"<div class='note'>"
        f"Cross-product: {len(rp_branches)} unique R&hArr;P branches "
        f"&times; {len(target_keys)} unique R&hArr;{right_short} branches "
        f"= {len(scored)} total scorings. Each R&hArr;P branch yields its own "
        f"(broken, formed, core_R, &delta;<sub>RP</sub>), which drives the verifier S "
        f"for every R&hArr;{right_short} mapping under it.<br>"
        f"&#9733; = highest-S combination. "
        f"&#9654; trace = which R&hArr;{right_short} mapping the slider above "
        f"is showing. {n_with_trace}/{len(target_keys)} R&hArr;{right_short} "
        f"mappings have an event trace recovered from {n_seeds} seeds.</div>"
        f"{truncation_note}"
        "<table class='mt'>"
        "<tr><th>rank</th><th>R&hArr;P</th><th>S</th><th>&beta;</th><th>&rho;</th>"
        "<th>&kappa;</th><th>br/fm(RP)</th><th>R27&rarr;T</th><th>R22&rarr;T</th>"
        "<th>trace?</th><th>events</th>"
        f"<th>witness mapping (R(el) &rarr; {label}(el))</th></tr>"
        f"{mappings_rows}"
        "</table>"
    )

    n_imag = int(np.sum(freqs < 0))
    title = (
        f"{step}  {label}  &mdash;  {len(scored)} unique mappings  "
        f"(highest S = {best['score']:.4f}: "
        f"&beta;={best['beta']:.3f}  &rho;={best['rho']:.3f}  "
        f"&kappa;={best['kappa']:.3f}  n_imag={n_imag})"
    )
    html = HTML.format(
        title=title,
        left_title="Reactant",
        right_title=right_title,
        left_short="R",
        right_short=right_short,
        xyzR_json=json.dumps(write_xyz_str(elR, xyzR, comment="R")),
        xyzP_json=json.dumps(write_xyz_str(elT, xyzT, comment=label)),
        events_json=json.dumps(events),
        wboR_json=json.dumps(wboR.tolist()),
        wboP_json=json.dumps(wboT.tolist()),
        elements_R_json=json.dumps(elR),
        elements_P_json=json.dumps(elT),
        mappings_summary_html=mappings_summary_html,
    )
    html = patch_html_for_pq(html)
    out_path = out_dir / f"{label}.html"
    out_path.write_text(html)
    print(f"\nwrote {out_path}")


if __name__ == "__main__":
    main()
