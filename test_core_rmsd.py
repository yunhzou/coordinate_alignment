"""
Take the first 5 elementary steps and, for each, align the proposed
candidate TS to the reactant via the multi-seed mapping algorithm.
Then compare its core atoms (atoms touching R->P broken/formed bonds)
against the ground-truth TS's core atoms via Kabsch + RMSD.

Pipeline per step:
  1. xtb on R, P, GT-TS, proposed TS (all cached when possible)
  2. align P -> R, GT-TS -> R, proposed-TS -> R (multi-seed)
  3. classify R->P bonds; core atoms = atoms touching any flagged bond
  4. extract core-atom coords from each in R-frame indexing
  5. Kabsch align proposed-TS core to GT-TS core; report RMSD
"""
from __future__ import annotations
import json
import re
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from rxn_core_frag import (
    run_xtb, build_graph, find_islands, expand_mapping,
    classify_bonds, kabsch, _generate_seed_orders,
)
from build_tsdisco_viewer import step_inputs


TSDISCO_INDEX = Path(
    "/Users/yunhengz/empty_for_claude/tsdisco_benchmark_visualization_plain_portable/index.html"
)
WORK = Path("/Users/yunhengz/empty_for_claude/rxn_core/work_core_rmsd")
WORK.mkdir(parents=True, exist_ok=True)


def load_data():
    html = TSDISCO_INDEX.read_text()
    m = re.search(r"const TSDISCO_DATA = (\{.*?\});\s*\n", html, re.DOTALL)
    return json.loads(m.group(1))


def best_mapping(g_R, g_P, wboR, wboP, n_seeds=10):
    orders = _generate_seed_orders(g_R, n_seeds)
    best = None
    for order in orders:
        m, _ = find_islands(g_R, g_P, seed_order=order)
        m = expand_mapping(m, g_R, g_P)
        br, fm, _, _ = classify_bonds(m, wboR, wboP)
        score = (len(br) + len(fm), -len(m))
        if best is None or score < best[0]:
            best = (score, m)
    return best[1]


def reindex_coords(map_target_to_src, src_coords, n_target, fallback):
    """Return target-indexed coords; unmapped fall back."""
    out = np.array(fallback, dtype=float).copy()
    for tgt, src in map_target_to_src.items():
        out[tgt] = src_coords[src]
    return out


def core_atoms_R_frame(mapping_R_to_P, broken, formed):
    """Indices in R-frame of atoms touching any flagged bond.
    `broken` items are R-frame already; `formed` items are P-frame -- map back."""
    inv = {v: k for k, v in mapping_R_to_P.items()}
    core = set()
    for (i, j, _, _) in broken:
        core.add(i); core.add(j)
    for (ip, jp, _, _) in formed:
        if ip in inv: core.add(inv[ip])
        if jp in inv: core.add(inv[jp])
    return sorted(core)


def kabsch_rmsd_full(coords_a_full, coords_b_full, core_indices):
    """Kabsch using FULL molecule, then RMSD on core atoms only."""
    R, t = kabsch(coords_b_full, coords_a_full)
    aligned_b = (R @ coords_b_full.T).T + t
    a_core = coords_a_full[core_indices]
    b_core = aligned_b[core_indices]
    diff = a_core - b_core
    return float(np.sqrt(np.mean(np.sum(diff * diff, axis=1))))


def kabsch_rmsd_spectator(coords_a_full, coords_b_full, core_indices):
    """Kabsch using SPECTATOR atoms only (atoms not in core), then
    measure RMSD on the core atoms after applying that alignment.
    This anchors the alignment to the unchanging scaffold and reports
    the geometric error in the reactive region."""
    n = coords_a_full.shape[0]
    core_set = set(core_indices)
    spec = [i for i in range(n) if i not in core_set]
    if len(spec) < 3:
        return None  # Kabsch underdetermined
    R, t = kabsch(coords_b_full[spec], coords_a_full[spec])
    aligned_b = (R @ coords_b_full.T).T + t
    a_core = coords_a_full[core_indices]
    b_core = aligned_b[core_indices]
    diff = a_core - b_core
    return float(np.sqrt(np.mean(np.sum(diff * diff, axis=1))))


def direct_rmsd(coords_a, coords_b):
    diff = coords_a - coords_b
    return float(np.sqrt(np.mean(np.sum(diff * diff, axis=1))))


def get_role_xyz(step, role):
    for ist in step["input_structures"]:
        if ist["role"] == role:
            return ist["xyz"]
    return None


def write_xyz(path, text):
    Path(path).write_text(text)


def analyze_step(step, n_cands=3):
    name = f"{step['dataset']}/{step['step_id']}"
    sanitized = re.sub(r"[^A-Za-z0-9._-]", "_", name)
    wd = WORK / sanitized
    wd.mkdir(parents=True, exist_ok=True)
    chg = step.get("charge", 0) or 0
    uhf = max(0, (step.get("multiplicity", 1) or 1) - 1)

    # R, P (concatenated multi-fragment)
    rxyz, pxyz, _, _ = step_inputs(step)
    write_xyz(wd / "reactant.xyz", rxyz)
    write_xyz(wd / "product.xyz", pxyz)

    # GT TS (reference_ts; falls back to reference_ts_xtb_optimized)
    gt_ts_xyz = get_role_xyz(step, "reference_ts") or get_role_xyz(step, "reference_ts_xtb_optimized")
    if gt_ts_xyz is None:
        return f"{name}: no reference_ts in input_structures"
    write_xyz(wd / "gt_ts.xyz", gt_ts_xyz)

    # xtb on R, P, GT-TS
    elR, xyzR, wboR = run_xtb(wd / "reactant.xyz", wd / "R", charge=chg, uhf=uhf)
    elP, xyzP, wboP = run_xtb(wd / "product.xyz", wd / "P", charge=chg, uhf=uhf)
    elGT, xyzGT, wboGT = run_xtb(wd / "gt_ts.xyz", wd / "GT_TS", charge=chg, uhf=uhf)

    g_R = build_graph(elR, wboR); g_P = build_graph(elP, wboP); g_GT = build_graph(elGT, wboGT)

    # Align all three to R
    map_R_to_P = best_mapping(g_R, g_P, wboR, wboP)
    map_R_to_GT = best_mapping(g_R, g_GT, wboR, wboGT)

    # Identify core atoms in R-frame from R->P chemistry
    br, fm, _, _ = classify_bonds(map_R_to_P, wboR, wboP)
    core = core_atoms_R_frame(map_R_to_P, br, fm)

    # Reindex GT-TS coords into R-frame
    n_R = len(elR)
    gt_coords_in_R = reindex_coords(map_R_to_GT, xyzGT, n_R, xyzR)

    # Now process top-K candidates
    cands = step.get("candidates", [])[:n_cands]
    out_rows = []
    for ci, cand in enumerate(cands):
        ts_xyz = cand.get("xyz", "")
        if not ts_xyz:
            out_rows.append((cand.get("rank", ci), "no xyz", None, None, None))
            continue
        write_xyz(wd / f"cand_{ci}.xyz", ts_xyz)
        try:
            elTS, xyzTS, wboTS = run_xtb(wd / f"cand_{ci}.xyz", wd / f"CAND_{ci}",
                                          charge=chg, uhf=uhf)
        except Exception as e:
            out_rows.append((cand.get("rank", ci), f"xtb fail: {e}", None, None, None))
            continue
        g_TS = build_graph(elTS, wboTS)
        map_R_to_TS = best_mapping(g_R, g_TS, wboR, wboTS)
        ts_coords_in_R = reindex_coords(map_R_to_TS, xyzTS, n_R, xyzR)

        # Three metrics:
        #   full→GT: Kabsch on full molecule, RMSD on core (compromise alignment)
        #   spec→GT: Kabsch on spectator atoms only, RMSD on core (chemistry metric)
        #   dir →GT: no re-alignment, direct RMSD on core (pre-aligned input)
        if core:
            full_GT = kabsch_rmsd_full(gt_coords_in_R, ts_coords_in_R, core)
            spec_GT = kabsch_rmsd_spectator(gt_coords_in_R, ts_coords_in_R, core)
            dir_GT  = direct_rmsd(gt_coords_in_R[core], ts_coords_in_R[core])
            spec_R  = kabsch_rmsd_spectator(xyzR, ts_coords_in_R, core)
        else:
            full_GT = spec_GT = dir_GT = spec_R = None
        out_rows.append((cand.get("rank", ci), len(map_R_to_TS),
                         spec_R, full_GT, spec_GT, dir_GT))
    return name, len(elR), len(core), out_rows


def main():
    data = load_data()
    steps = data["steps"][:5]

    for step in steps:
        result = analyze_step(step, n_cands=3)
        if isinstance(result, str):
            print(result); continue
        name, natoms, n_core, rows = result
        print(f"\n=== {name}  N={natoms}  core={n_core} ===")
        print(f"  {'rank':>4}  {'mapped':>6}  {'spec→R':>7}  {'full→GT':>8}  {'spec→GT':>8}  {'dir→GT':>7}")
        for rank, mapped, spec_R, full_GT, spec_GT, dir_GT in rows:
            f = lambda x: f"{x:.3f}" if isinstance(x, float) else str(x) if x is not None else "—"
            print(f"  {rank:>4}  {mapped:>6}  {f(spec_R):>7}  {f(full_GT):>8}  {f(spec_GT):>8}  {f(dir_GT):>7}")


if __name__ == "__main__":
    main()
