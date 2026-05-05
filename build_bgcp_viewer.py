"""
Build viewers for Benchmark_Guesses_Collective_Package (BGCP).

Each step dir under BGCP_ROOT has:
  reactants/      — one or more reactant xyz fragments
  products/       — one or more product xyz fragments
  groundtruth/    — ground-truth TS xyz
  initial_guess/  — N candidate TS xyz files (typically 20)

Charge / multiplicity are not in the xyz files themselves; we cross-
reference the tsdisco index.html data by step_id to recover them
(falling back to 0 / 1 if not found).

Outputs:
  out/bgcp_viewer.html    — 2-panel R / P with broken/formed dashes
  out/bgcp_rtsp_viewer.html — 3-panel R / TS-rank-0 / P
  out/bgcp_bonds.csv      — per-step bond summary
"""
from __future__ import annotations
import argparse
import json
import re
import time
import traceback
from pathlib import Path

import numpy as np

from rxn_core_frag import (
    run_xtb, build_graph, find_islands, expand_mapping,
    classify_bonds, write_xyz_str, _generate_seed_orders,
)
from build_tsdisco_viewer import _parse_xyz_text, concat_xyz, HTML as TSDISCO_HTML
from build_rtsp_viewer import (
    HTML as RTSP_HTML, all_bonds, changed_bonds, reindex_xyz_to_target,
    best_mapping,
)


BGCP_ROOT = Path("/Users/yunhengz/empty_for_claude/rxn_core/Benchmark_Guesses_Collective_Package")
TSDISCO_INDEX = Path(
    "/Users/yunhengz/empty_for_claude/tsdisco_benchmark_visualization_plain_portable/index.html"
)
ROOT = Path("/Users/yunhengz/empty_for_claude/rxn_core")
OUT = ROOT / "out"
WORK = ROOT / "work_bgcp"
OUT.mkdir(parents=True, exist_ok=True)
WORK.mkdir(parents=True, exist_ok=True)


def load_tsdisco_chg_uhf_lookup():
    """step_id -> (charge, uhf) by parsing tsdisco index.html"""
    if not TSDISCO_INDEX.exists():
        return {}
    html = TSDISCO_INDEX.read_text()
    m = re.search(r"const TSDISCO_DATA = (\{.*?\});\s*\n", html, re.DOTALL)
    if m is None:
        return {}
    data = json.loads(m.group(1))
    out = {}
    for s in data["steps"]:
        chg = s.get("charge", 0) or 0
        mult = s.get("multiplicity", 1) or 1
        out[s["step_id"]] = (int(chg), max(0, int(mult) - 1))
    return out


def list_step_dirs():
    return sorted(d for d in BGCP_ROOT.iterdir()
                  if d.is_dir() and not d.name.startswith("."))


def read_xyzs(dir_path):
    """Return concatenated xyz text for all xyz files in dir_path
    (translated apart by 10 Å gap if multiple)."""
    files = sorted(dir_path.glob("*.xyz"))
    if not files:
        return None
    texts = [f.read_text() for f in files]
    return concat_xyz(texts)


def list_initial_guesses(step_dir):
    g = sorted((step_dir / "initial_guess").glob("*.xyz"))
    return g


def iter_num(path):
    """Extract iter number from filename like '...iter12_xxxx.xyz'."""
    m = re.search(r"_iter(\d+)_", path.name)
    return int(m.group(1)) if m else 9999


def bond_table_html(rows):
    if not rows: return "<i>none</i>"
    head = "<tr><th>i</th><th>j</th><th>WBO_R</th><th>WBO_P</th></tr>"
    body = []
    for i, j, wR, wP in rows:
        wR_s = "—" if wR is None else f"{wR:.2f}"
        wP_s = "—" if wP is None else f"{wP:.2f}"
        body.append(f"<tr><td>{i}</td><td>{j}</td><td>{wR_s}</td><td>{wP_s}</td></tr>")
    return f"<table class='bondtab'>{head}{''.join(body)}</table>"


def analyze_rp(step_dir, chg, uhf, name):
    sanitized = re.sub(r"[^A-Za-z0-9._-]", "_", name)
    wd = WORK / sanitized
    wd.mkdir(parents=True, exist_ok=True)
    rxyz = read_xyzs(step_dir / "reactants")
    pxyz = read_xyzs(step_dir / "products")
    if rxyz is None or pxyz is None:
        raise RuntimeError("missing reactant or product")
    (wd / "reactant.xyz").write_text(rxyz)
    (wd / "product.xyz").write_text(pxyz)
    elR, xyzR, wboR = run_xtb(wd / "reactant.xyz", wd / "R", charge=chg, uhf=uhf)
    elP, xyzP, wboP = run_xtb(wd / "product.xyz", wd / "P", charge=chg, uhf=uhf)
    g_R = build_graph(elR, wboR); g_P = build_graph(elP, wboP)
    mapping = best_mapping(g_R, g_P, wboR, wboP)
    broken, formed, _, _ = classify_bonds(mapping, wboR, wboP)
    return dict(
        wd=wd, elR=elR, elP=elP, xyzR=xyzR, xyzP=xyzP, wboR=wboR, wboP=wboP,
        g_R=g_R, g_P=g_P, mapping=mapping, broken=broken, formed=formed,
    )


def analyze_step_rp(step_dir):
    name = step_dir.name
    chg, uhf = LOOKUP.get(name, (0, 0))
    res = analyze_rp(step_dir, chg, uhf, name)
    elR = res["elR"]; elP = res["elP"]; m = res["mapping"]
    inv = {v: k for k, v in m.items()}
    elP_in_R, xyzP_in_R = reindex_xyz_to_target(m, elP, res["xyzP"], len(elR), elR, res["xyzR"])
    return {
        "name": name,
        "xyzR": write_xyz_str(elR, res["xyzR"], comment="reactant"),
        "xyzP": write_xyz_str(elP_in_R, xyzP_in_R, comment="product (R-frame)"),
        "broken_idx": [[i, j] for (i, j, _, _) in res["broken"]],
        "formed_idx": [[inv[i] if i in inv else i, inv[j] if j in inv else j]
                        for (i, j, _, _) in res["formed"]],
        "natoms": len(elR),
        "n_broken": len(res["broken"]),
        "n_formed": len(res["formed"]),
        "n_mapped": len(m),
        "broken_table": bond_table_html(res["broken"]),
        "formed_table": bond_table_html(res["formed"]),
        "charge": chg, "uhf": uhf,
        "mechanism": "", "step": "",
    }


def analyze_step_rtsp(step_dir, ts_path):
    """3-panel: R / TS / P, each in its OWN native coordinate frame.
    Bond-change cylinders are mapped through R↔P and R↔TS alignments
    so they render at the correct atoms on each panel without needing
    to reindex the geometries. This avoids Frankenstein displays when
    some atoms don't map (their fallback coords would be wrong)."""
    name = step_dir.name
    chg, uhf = LOOKUP.get(name, (0, 0))
    base = analyze_rp(step_dir, chg, uhf, name)
    wd = base["wd"]
    elR = base["elR"]; xyzR = base["xyzR"]; wboR = base["wboR"]
    elP = base["elP"]; xyzP = base["xyzP"]; wboP = base["wboP"]
    map_R_to_P = base["mapping"]

    # TS xtb in its own native frame
    sanitized_ts = re.sub(r"[^A-Za-z0-9._-]", "_", ts_path.stem)[:80]
    ts_local = wd / f"{sanitized_ts}.xyz"
    ts_local.write_text(ts_path.read_text())
    ts_workdir = wd / f"TS_{sanitized_ts}"
    elTS, xyzTS, wboTS = run_xtb(ts_local, ts_workdir, charge=chg, uhf=uhf)
    g_TS = build_graph(elTS, wboTS)
    map_R_to_TS = best_mapping(base["g_R"], g_TS, wboR, wboTS)

    # Bond-change classification still done in R-frame (R↔P semantics).
    # We project P-side WBO into an n_R-by-n_R matrix using the R→P map
    # ONLY for the chemistry classification; we do NOT use it for
    # rendering coordinates.
    n_R = len(elR)
    wboP_r = np.zeros_like(wboR)
    for ri, pi in map_R_to_P.items():
        for rj, pj in map_R_to_P.items():
            if ri < rj:
                wboP_r[ri, rj] = wboP[pi, pj]; wboP_r[rj, ri] = wboP[pi, pj]
    wboTS_r = np.zeros_like(wboR)
    for ri, ti in map_R_to_TS.items():
        for rj, tj in map_R_to_TS.items():
            if ri < rj:
                wboTS_r[ri, rj] = wboTS[ti, tj]; wboTS_r[rj, ri] = wboTS[ti, tj]

    rt = changed_bonds(wboR, wboTS_r)
    tp = changed_bonds(wboTS_r, wboP_r)
    rp = changed_bonds(wboR, wboP_r)
    inflight = []
    for (i, j, wR_, wP_, d) in rp:
        wT = float(wboTS_r[i, j])
        lo = min(wR_, wP_); hi = max(wR_, wP_)
        if lo + 0.15 < wT < hi - 0.15:
            inflight.append((i, j, wR_, wT, wP_))

    # Translate R-indexed bond-change pairs to P-indexed and TS-indexed
    # for rendering on each panel. Bonds whose endpoints aren't both
    # mapped get dropped (they'd have nowhere to render).
    def to_p(rp_list):
        out = []
        for (i, j, wA, wB, d) in rp_list:
            if i in map_R_to_P and j in map_R_to_P:
                out.append([map_R_to_P[i], map_R_to_P[j], wA, wB, d])
        return out
    def to_ts(rp_list):
        out = []
        for (i, j, wA, wB, d) in rp_list:
            if i in map_R_to_TS and j in map_R_to_TS:
                out.append([map_R_to_TS[i], map_R_to_TS[j], wA, wB, d])
        return out

    return {
        "name": name,
        "natoms": len(elR),
        "elements": elR,
        # Native-frame coordinates for each panel (no reindexing)
        "xyzR": write_xyz_str(elR, xyzR, comment="reactant"),
        "xyzTS": write_xyz_str(elTS, xyzTS, comment=f"TS: {ts_path.name}"),
        "xyzP": write_xyz_str(elP, xyzP, comment="product"),
        # Bond-change lists per panel (in panel-native indexing)
        "rt_changes": rt,                # R-indexed (renders on R panel)
        "tp_changes": tp,                # R-indexed
        "rp_changes": rp,                # R-indexed (renders on R panel)
        "rp_changes_P": to_p(rp),        # P-indexed (renders on P panel)
        "rp_changes_TS": to_ts(rp),      # TS-indexed (renders on TS panel)
        "inflight": inflight,            # R-indexed
        "inflight_TS": [[map_R_to_TS[i], map_R_to_TS[j], wR_, wT, wP_]
                         for (i, j, wR_, wT, wP_) in inflight
                         if i in map_R_to_TS and j in map_R_to_TS],
        "n_R_to_TS": len(rt),
        "n_TS_to_P": len(tp),
        "n_R_to_P": len(rp),
        "n_inflight": len(inflight),
        "n_mapped_R_to_P": len(map_R_to_P),
        "n_mapped_R_to_TS": len(map_R_to_TS),
        "charge": chg, "uhf": uhf,
        "ts_source": ts_path.name,
    }


# Cache the tsdisco lookup at module load
LOOKUP = load_tsdisco_chg_uhf_lookup()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--start", type=int, default=0)
    ap.add_argument("--mode", choices=["rp", "rtsp", "both"], default="both")
    ap.add_argument("--ts-pick", choices=["iter1", "groundtruth"], default="iter1",
                    help="rtsp mode: which TS to use (lowest-iter initial guess, or groundtruth)")
    args = ap.parse_args()

    step_dirs = list_step_dirs()
    if args.limit is not None:
        step_dirs = step_dirs[args.start:args.start + args.limit]
    else:
        step_dirs = step_dirs[args.start:]
    print(f"[bgcp] {len(step_dirs)} steps; mode={args.mode}; lookup_size={len(LOOKUP)}")

    rp_data = {}
    rtsp_data = {}
    for k, sd in enumerate(step_dirs, 1):
        name = sd.name
        t = time.time()
        try:
            if args.mode in ("rp", "both"):
                d = analyze_step_rp(sd)
                rp_data[name] = d
            if args.mode in ("rtsp", "both"):
                if args.ts_pick == "groundtruth":
                    gt_files = sorted((sd / "groundtruth").glob("*.xyz"))
                    if not gt_files: raise RuntimeError("no groundtruth ts")
                    ts_path = gt_files[0]
                else:
                    guesses = sorted(list_initial_guesses(sd), key=iter_num)
                    if not guesses: raise RuntimeError("no initial guesses")
                    ts_path = guesses[0]
                d3 = analyze_step_rtsp(sd, ts_path)
                rtsp_data[name] = d3
            br, fm = (rp_data.get(name, {}).get("n_broken", 0),
                      rp_data.get(name, {}).get("n_formed", 0))
            print(f"[{k:>3}/{len(step_dirs)}]  {time.time()-t:5.1f}s  OK   "
                  f"{name:<70s}  br/fm={br}/{fm}")
        except Exception as e:
            print(f"[{k:>3}/{len(step_dirs)}]  FAIL {name}: {e}")
            traceback.print_exc()

    if args.mode in ("rp", "both") and rp_data:
        html = TSDISCO_HTML.replace("__DATA__", json.dumps(rp_data))
        out = OUT / "bgcp_viewer.html"
        out.write_text(html)
        print(f"[bgcp] wrote {out}  ({len(rp_data)} steps)")

        # CSV
        import csv
        csv_path = OUT / "bgcp_bonds.csv"
        with csv_path.open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=[
                "step_id", "n_atoms", "n_mapped", "n_broken", "n_formed",
                "charge", "multiplicity", "broken_bonds", "formed_bonds"
            ])
            w.writeheader()
            for name, d in sorted(rp_data.items()):
                els = d["xyzR"].strip().splitlines()
                n = int(els[0])
                rE = [ln.split()[0] for ln in els[2:2+n]]
                pE_lines = d["xyzP"].strip().splitlines()
                pE = [ln.split()[0] for ln in pE_lines[2:2+n]]
                br_str = ";".join(f"{i}({rE[i]})-{j}({rE[j]})" for (i, j) in d["broken_idx"])
                fm_str = ";".join(f"{i}({pE[i]})-{j}({pE[j]})" for (i, j) in d["formed_idx"])
                w.writerow({
                    "step_id": name,
                    "n_atoms": d["natoms"], "n_mapped": d["n_mapped"],
                    "n_broken": d["n_broken"], "n_formed": d["n_formed"],
                    "charge": d["charge"], "multiplicity": d["uhf"] + 1,
                    "broken_bonds": br_str, "formed_bonds": fm_str,
                })
        print(f"[bgcp] wrote {csv_path}")

    if args.mode in ("rtsp", "both") and rtsp_data:
        html = RTSP_HTML.replace("__DATA__", json.dumps(rtsp_data))
        out = OUT / "bgcp_rtsp_viewer.html"
        out.write_text(html)
        print(f"[bgcp] wrote {out}  ({len(rtsp_data)} steps)")


if __name__ == "__main__":
    main()
