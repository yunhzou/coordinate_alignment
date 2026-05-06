"""
Per-step core-atom mode analysis. For each BGCP step:
  1. Run R/P alignment via cached xtb to get core atoms (atoms touching
     broken or formed bonds, expressed in R-frame).
  2. For each TS structure (groundtruth + 20 initial guesses):
       a. Parse vibration modes from cached xtb Hessian (g98.out).
       b. Re-derive R↔TS atom mapping (cache hit) and reindex modes
          into R-frame so atom i in mode array corresponds to R atom i.
       c. For every imaginary mode, compute core_fraction =
          ||disp on core atoms||² / ||disp on all atoms||²
  3. Write CSV `out/mode_analysis/<step>.csv` with one row per
     (TS_label, mode_idx) for imaginary modes, ranked by core_fraction.

Multiprocessing.Pool parallelizes across steps.

Output dir: out/mode_analysis/
"""
from __future__ import annotations
import argparse
import csv
import multiprocessing as mp
import os
import re
import sys
import time
import traceback
from pathlib import Path

import numpy as np

from rxn_core_pq import align_from_arrays
from rxn_core_frag import parse_xyz, build_graph, classify_bonds
from bgcp_io import BGCP_ROOT, LOOKUP, list_step_dirs
from align_bgcp_coords import load_cached_xtb, fill_unmapped_greedy


WORK_MODES = Path("/Users/yunhengz/empty_for_claude/rxn_core/work_modes")
OUT_CSV_DIR = Path("/Users/yunhengz/empty_for_claude/rxn_core/out/mode_analysis")
OUT_CSV_DIR.mkdir(parents=True, exist_ok=True)


def parse_g98_modes(path):
    """Parse xtb's g98.out (Gaussian-style normal-mode block).
    Returns (freqs_array, modes_array) where
      freqs_array shape (n_modes,)
      modes_array shape (n_modes, n_atoms, 3) — Cartesian displacements."""
    text = Path(path).read_text().splitlines()
    # Determine n_atoms from the Standard orientation block
    n_atoms = 0
    for i, line in enumerate(text):
        if "Standard orientation" in line:
            j = i + 5
            while j < len(text) and text[j].strip() and not text[j].strip().startswith("---"):
                if text[j].split()[0].isdigit():
                    n_atoms += 1
                j += 1
            break
    freqs = []
    modes = []
    i = 0
    while i < len(text):
        line = text[i]
        if line.lstrip().startswith("Frequencies --"):
            parts = line.split("--", 1)[1].split()
            block_freqs = [float(x) for x in parts]
            n_block = len(block_freqs)
            j = i + 1
            while j < len(text) and not text[j].lstrip().startswith("Atom"):
                j += 1
            block_modes = [[] for _ in range(n_block)]
            for a in range(n_atoms):
                row = text[j + 1 + a].split()
                vals = [float(x) for x in row[2:2 + 3 * n_block]]
                for b in range(n_block):
                    block_modes[b].extend(vals[3*b:3*b + 3])
            freqs.extend(block_freqs)
            modes.extend(block_modes)
            i = j + 1 + n_atoms
        else:
            i += 1
    freqs = np.array(freqs)
    if modes:
        modes_arr = np.array(modes).reshape(len(modes), n_atoms, 3)
    else:
        modes_arr = np.zeros((0, n_atoms, 3))
    return freqs, modes_arr


def core_atoms_in_R_frame(mapping_R_to_P, broken, formed):
    """Atoms touching broken (R indices) or formed (P, mapped back via inv).
    Returns sorted set of R-frame indices."""
    inv = {v: k for k, v in mapping_R_to_P.items()}
    core = set()
    for (i, j, _, _) in broken:
        core.add(int(i)); core.add(int(j))
    for (ip, jp, _, _) in formed:
        if ip in inv: core.add(int(inv[ip]))
        if jp in inv: core.add(int(inv[jp]))
    return sorted(core)


def reindex_modes_to_R(modes_TS, mapping_R_to_TS, n_R):
    """modes_TS shape (n_modes, n_atoms_TS, 3). Returns shape
    (n_modes, n_R, 3) where row i is the displacement of R-atom i
    (taken from modes_TS at the mapped TS index, zero if unmapped)."""
    n_modes = modes_TS.shape[0]
    out = np.zeros((n_modes, n_R, 3))
    for r, t in mapping_R_to_TS.items():
        out[:, r, :] = modes_TS[:, t, :]
    return out


def list_ts_caches(modes_dir):
    """List (label, hess_dir, sp_dir) tuples for every TS we have cached."""
    out = []
    if (modes_dir / "hess_groundtruth").exists():
        out.append(("groundtruth", modes_dir / "hess_groundtruth",
                    modes_dir / "sp_groundtruth"))
    for hess in sorted(modes_dir.glob("hess_iter*")):
        if not hess.is_dir(): continue
        m = re.match(r"hess_iter(\d+)$", hess.name)
        if not m: continue
        iter_n = int(m.group(1))
        sp = modes_dir / f"sp_iter{iter_n}"
        if not sp.exists(): continue
        out.append((f"iter{iter_n}", hess, sp))
    return out


def process_step(name):
    """Worker: compute core-fraction rows for every TS / imag mode of
    one step. Returns (name, list_of_dicts) on success or
    (name, exception_str) on failure."""
    step_modes_dir = WORK_MODES / name
    if not (step_modes_dir / "R" / "wbo").exists():
        raise RuntimeError(f"missing R cache for {name}")
    if not (step_modes_dir / "P" / "wbo").exists():
        raise RuntimeError(f"missing P cache for {name}")

    # R/P alignment → core atoms in R-frame
    elR, xyzR, wboR, _ = load_cached_xtb(step_modes_dir / "R")
    elP, xyzP, wboP, _ = load_cached_xtb(step_modes_dir / "P")
    rp_res = align_from_arrays(elR, xyzR, wboR, elP, xyzP, wboP)
    mapping_RP = dict(rp_res['mapping'])
    core_R = core_atoms_in_R_frame(mapping_RP, rp_res['broken'], rp_res['formed'])
    n_R = len(elR)

    rows = []
    for label, hess_dir, sp_dir in list_ts_caches(step_modes_dir):
        g98 = hess_dir / "g98.out"
        if not g98.exists():
            continue
        try:
            freqs, modes_TS = parse_g98_modes(g98)
        except Exception as e:
            rows.append({
                'step': name, 'ts_label': label, 'mode_idx': -1,
                'freq': None, 'core_fraction': None,
                'mode_rank': None, 'n_imag': None, 'error': f'parse: {e}',
            })
            continue
        if modes_TS.shape[0] == 0:
            continue

        # R↔TS mapping (cache hit on sp_dir's xyz/wbo)
        elT, xyzT, wboT, _ = load_cached_xtb(sp_dir)
        try:
            ts_res = align_from_arrays(elR, xyzR, wboR, elT, xyzT, wboT)
            mapping_RT = dict(ts_res['mapping'])
            mapping_RT = fill_unmapped_greedy(elR, xyzR, elT, xyzT, mapping_RT)
        except Exception as e:
            rows.append({
                'step': name, 'ts_label': label, 'mode_idx': -1,
                'freq': None, 'core_fraction': None,
                'mode_rank': None, 'n_imag': None, 'error': f'align: {e}',
            })
            continue

        modes_R = reindex_modes_to_R(modes_TS, mapping_RT, n_R)
        # Per-mode total energy (sum of squared displacements)
        sq = (modes_R ** 2).sum(axis=2)         # (n_modes, n_atoms)
        total = sq.sum(axis=1)                  # (n_modes,)
        core_e = sq[:, core_R].sum(axis=1) if core_R else np.zeros(modes_R.shape[0])
        fraction = np.where(total > 1e-12, core_e / total, 0.0)

        imag_idx = np.where(freqs < 0)[0]
        n_imag = len(imag_idx)
        # Rank imag modes by core_fraction desc
        order = sorted(imag_idx, key=lambda k: -fraction[k])
        for rank, k in enumerate(order):
            rows.append({
                'step': name,
                'ts_label': label,
                'mode_idx': int(k),
                'freq': float(freqs[k]),
                'core_fraction': float(fraction[k]),
                'mode_rank': rank,
                'n_imag': n_imag,
                'n_modes_total': int(modes_R.shape[0]),
                'n_core_atoms': len(core_R),
                'core_atoms': ','.join(map(str, core_R)),
                'error': '',
            })

    # Per-step CSV
    sanitized = re.sub(r"[^A-Za-z0-9._-]", "_", name)
    csv_path = OUT_CSV_DIR / f"{sanitized}.csv"
    with csv_path.open('w', newline='') as f:
        fieldnames = ['step', 'ts_label', 'mode_idx', 'freq', 'core_fraction',
                      'mode_rank', 'n_imag', 'n_modes_total', 'n_core_atoms',
                      'core_atoms', 'error']
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)
    return rows


def _safe(name):
    try:
        rows = process_step(name)
        return (name, True, len(rows))
    except Exception as e:
        return (name, False, f"{type(e).__name__}: {e}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--workers', type=int, default=os.cpu_count())
    ap.add_argument('--steps', nargs='+', default=None)
    ap.add_argument('--limit', type=int, default=None)
    args = ap.parse_args()

    all_steps = [d.name for d in list_step_dirs()]
    if args.steps:
        steps = [s for s in all_steps if s in set(args.steps)]
    elif args.limit:
        steps = all_steps[:args.limit]
    else:
        steps = all_steps
    steps = [s for s in steps if (WORK_MODES / s / "R" / "wbo").exists()]

    print(f"Mode analysis on {len(steps)} steps with {args.workers} workers")
    t0 = time.time()
    n_done = n_err = 0
    with mp.Pool(args.workers) as pool:
        for i, (name, ok, payload) in enumerate(pool.imap_unordered(_safe, steps), 1):
            if ok:
                print(f"[{i:3d}/{len(steps)}] {name[:55]:55s}  rows={payload}")
                n_done += 1
            else:
                print(f"[{i:3d}/{len(steps)}] {name[:55]:55s}  ERROR: {payload}")
                n_err += 1
            sys.stdout.flush()
    print(f"\nDone: {n_done} ok, {n_err} errors in {time.time()-t0:.1f}s")
    print(f"CSVs: {OUT_CSV_DIR}")


if __name__ == "__main__":
    main()
