"""
Align all BGCP structures (reactants/products/groundtruth/initial_guess)
to a single shared atom indexing per step. The reactant ordering is the
canonical frame; products and every TS (groundtruth + initial guesses)
are reindexed so atom i refers to the same chemical atom across R, TS, P.

Uses pre-computed xtb output cached in `work_modes/<step>/`:
  R/{reactant.xyz, wbo}                              — reactant
  P/{product.xyz, wbo}                               — product
  sp_groundtruth/{groundtruth__*.xyz, wbo}           — GT TS
  sp_iter1/{iter1__*.xyz, wbo} ... sp_iter20/...     — initial guesses

No xtb is invoked. Multiprocessing.Pool parallelizes across steps.

Output mirrors BGCP folder layout under
  /Users/yunhengz/empty_for_claude/rxn_core/Benchmark_Guesses_Coordinate_Aligned_Version

Usage:
  python align_bgcp_coords.py
  python align_bgcp_coords.py --workers 14
  python align_bgcp_coords.py --steps Jackie_TS_10 pr1.tempo_ts2
"""
from __future__ import annotations
import argparse
import multiprocessing as mp
import os
import re
import sys
import time
import traceback
from pathlib import Path

import numpy as np

from rxn_core_pq import align_from_arrays
from rxn_core_frag import write_xyz_str, parse_xyz
from bgcp_io import BGCP_ROOT, LOOKUP, list_step_dirs


OUT_ROOT = Path("/Users/yunhengz/empty_for_claude/rxn_core/"
                "Benchmark_Guesses_Coordinate_Aligned_Version")
WORK_MODES = Path("/Users/yunhengz/empty_for_claude/rxn_core/work_modes")


def load_cached_xtb(workdir):
    """Load (elements, coords, wbo) from a pre-computed xtb workdir.
    The workdir must contain exactly one *.xyz and a `wbo` file."""
    workdir = Path(workdir)
    xyz_files = [f for f in workdir.iterdir() if f.suffix == '.xyz']
    if len(xyz_files) != 1:
        raise RuntimeError(f"expected 1 xyz in {workdir}, found {len(xyz_files)}")
    elements, coords = parse_xyz(xyz_files[0])
    n = len(elements)
    wbo = np.zeros((n, n))
    wf = workdir / "wbo"
    if not wf.exists():
        raise RuntimeError(f"no wbo file in {workdir}")
    for ln in wf.read_text().splitlines():
        parts = ln.split()
        if len(parts) < 3:
            continue
        i, j = int(parts[0]) - 1, int(parts[1]) - 1
        v = float(parts[2])
        wbo[i, j] = v; wbo[j, i] = v
    return elements, coords, wbo, xyz_files[0]


def fill_unmapped_greedy(elR, xyzR, elT, xyzT, mapping):
    """Greedy nearest-element-match for atoms PQ left unmapped. Pairs
    each unmapped R-atom with the closest still-free same-element T-atom."""
    used_T = set(mapping.values())
    unmapped_R = [i for i in range(len(elR)) if i not in mapping]
    if not unmapped_R:
        return dict(mapping)
    free_by_el = {}
    for j in range(len(elT)):
        if j in used_T:
            continue
        free_by_el.setdefault(elT[j], []).append(j)
    out = dict(mapping)
    for i in unmapped_R:
        cands = free_by_el.get(elR[i], [])
        if not cands:
            continue
        d = [np.linalg.norm(xyzT[j] - xyzR[i]) for j in cands]
        best = cands[int(np.argmin(d))]
        out[i] = best
        cands.remove(best)
    return out


def reindex_to_R_frame(elR, xyzR, elT, xyzT, mapping):
    """Build (elements, coords) of length len(elR) in R atom order."""
    n = len(elR)
    out_el = list(elR)
    out_xyz = np.array(xyzR, dtype=float).copy()
    missing = []
    for i in range(n):
        j = mapping.get(i)
        if j is None:
            missing.append(i)
            continue
        out_el[i] = elT[j]
        out_xyz[i] = xyzT[j]
    return out_el, out_xyz, missing


def align_target_to_R(elR, xyzR, wboR, T_dir):
    """Align a target structure (loaded from T_dir cache) to R-frame.
    Returns dict with elements, coords, missing list, pq_mapped count."""
    elT, xyzT, wboT, _ = load_cached_xtb(T_dir)
    res = align_from_arrays(elR, xyzR, wboR, elT, xyzT, wboT, n_seeds=10)
    pq_mapping = dict(res['mapping'])
    full = fill_unmapped_greedy(elR, xyzR, elT, xyzT, pq_mapping)
    aligned_el, aligned_xyz, missing = reindex_to_R_frame(
        elR, xyzR, elT, xyzT, full)
    return {
        'elements': aligned_el,
        'coords': aligned_xyz,
        'missing': missing,
        'pq_mapped': len(pq_mapping),
        'fallback_added': len(full) - len(pq_mapping),
        'n_atoms': len(elR),
    }


def list_ts_caches(modes_dir):
    """List (label, cache_subdir, source_xyz) tuples for every TS in a
    work_modes step directory: sp_groundtruth + sp_iter1..sp_iter20.
    Source xyz refers to the TS's original filename in the BGCP step
    (used to preserve naming in the output)."""
    out = []
    for sp in sorted(modes_dir.glob("sp_*")):
        if not sp.is_dir():
            continue
        # The xyz file inside sp_* has the original BGCP filename
        xyzs = [f for f in sp.iterdir() if f.suffix == '.xyz']
        if not xyzs:
            continue
        out.append((sp.name, sp, xyzs[0].name))
    return out


def process_step(name):
    """Top-level worker (must be picklable for multiprocessing)."""
    step_dir = BGCP_ROOT / name
    modes_dir = WORK_MODES / name
    out_dir = OUT_ROOT / name

    # Caches must exist
    if not (modes_dir / "R" / "wbo").exists():
        raise RuntimeError(f"missing cached R wbo at {modes_dir}/R")
    if not (modes_dir / "P" / "wbo").exists():
        raise RuntimeError(f"missing cached P wbo at {modes_dir}/P")

    elR, xyzR, wboR, R_xyz_path = load_cached_xtb(modes_dir / "R")

    # Mirror BGCP layout
    for sub in ('reactants', 'products', 'groundtruth', 'initial_guess'):
        (out_dir / sub).mkdir(parents=True, exist_ok=True)

    # Reactant: copy as canonical frame
    (out_dir / "reactants" / "reactant_aligned.xyz").write_text(
        R_xyz_path.read_text())

    # Product
    p_res = align_target_to_R(elR, xyzR, wboR, modes_dir / "P")
    (out_dir / "products" / "product_aligned.xyz").write_text(
        write_xyz_str(p_res['elements'], p_res['coords'],
                      comment=f"product_aligned_to_R  pq_mapped={p_res['pq_mapped']}/{p_res['n_atoms']}  fallback={p_res['fallback_added']}  missing={len(p_res['missing'])}"))

    # TS: groundtruth + initial guesses
    gt_res = []; ig_res = []
    for label, sp_dir, src_name in list_ts_caches(modes_dir):
        try:
            res = align_target_to_R(elR, xyzR, wboR, sp_dir)
            xyz_str = write_xyz_str(
                res['elements'], res['coords'],
                comment=f"ts_aligned_to_R [{label}]  pq_mapped={res['pq_mapped']}/{res['n_atoms']}  fallback={res['fallback_added']}  missing={len(res['missing'])}")
            if label == "sp_groundtruth":
                # Use the source name without the sp_ prefix transformation
                # The original BGCP groundtruth file lives at step_dir/groundtruth/
                gt_files = sorted((step_dir / "groundtruth").glob("*.xyz"))
                target_name = gt_files[0].name if gt_files else src_name
                (out_dir / "groundtruth" / target_name).write_text(xyz_str)
                gt_res.append((target_name, res))
            else:  # sp_iterN
                # Match to original BGCP initial_guess filename
                # The src_name in the cache is e.g. "iter1__Jackie_TS_10_..."
                # Original is e.g. "Jackie_TS_10_benchmark_plain_iter1_..."
                ig_files = sorted((step_dir / "initial_guess").glob("*.xyz"))
                # Derive target name by matching iter number
                m = re.match(r"sp_iter(\d+)$", label)
                target = src_name  # fallback
                if m:
                    iter_n = int(m.group(1))
                    for f in ig_files:
                        mm = re.search(r"_iter(\d+)_", f.name)
                        if mm and int(mm.group(1)) == iter_n:
                            target = f.name; break
                (out_dir / "initial_guess" / target).write_text(xyz_str)
                ig_res.append((target, res))
        except Exception as e:
            (out_dir / ("groundtruth" if label == "sp_groundtruth" else "initial_guess")
             / f"{label}.ERROR").write_text(f"{e}\n{traceback.format_exc()}")

    return {
        'name': name,
        'p_missing': len(p_res['missing']),
        'p_fallback': p_res['fallback_added'],
        'n_atoms': p_res['n_atoms'],
        'gt_count': len(gt_res),
        'ig_count': len(ig_res),
        'ig_max_missing': max((len(x['missing']) for _, x in ig_res), default=0),
        'ig_max_fallback': max((x['fallback_added'] for _, x in ig_res), default=0),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--workers', type=int, default=os.cpu_count())
    ap.add_argument('--steps', nargs='+', default=None)
    ap.add_argument('--limit', type=int, default=None)
    args = ap.parse_args()

    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    all_steps = [d.name for d in list_step_dirs()]
    if args.steps:
        steps = [s for s in all_steps if s in set(args.steps)]
    elif args.limit:
        steps = all_steps[:args.limit]
    else:
        steps = all_steps

    # Filter to those that have a work_modes cache
    steps = [s for s in steps if (WORK_MODES / s / "R" / "wbo").exists()]

    print(f"Aligning {len(steps)} BGCP steps using {args.workers} workers")
    print(f"  cache:  {WORK_MODES}")
    print(f"  output: {OUT_ROOT}")
    print()

    t0 = time.time()
    n_done = n_err = 0
    with mp.Pool(args.workers) as pool:
        for i, result in enumerate(pool.imap_unordered(_safe_process, steps), 1):
            name, ok, payload = result
            if ok:
                r = payload
                print(f"[{i:3d}/{len(steps)}] {name[:55]:55s}  "
                      f"N={r['n_atoms']}  P_miss={r['p_missing']} fb={r['p_fallback']}  "
                      f"GT={r['gt_count']}  IG={r['ig_count']}  "
                      f"ig_miss_max={r['ig_max_missing']}  ig_fb_max={r['ig_max_fallback']}")
                n_done += 1
            else:
                print(f"[{i:3d}/{len(steps)}] {name[:55]:55s}  ERROR: {payload}")
                n_err += 1
            sys.stdout.flush()

    print()
    print(f"Done: {n_done} ok, {n_err} errors in {time.time()-t0:.1f}s")
    print(f"Aligned dataset: {OUT_ROOT}")


def _safe_process(name):
    try:
        return (name, True, process_step(name))
    except Exception as e:
        return (name, False, f"{type(e).__name__}: {e}")


if __name__ == "__main__":
    main()
