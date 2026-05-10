"""
External-dataset version of build_ranked_view_one_step.py.

Runs the full pipeline from scratch (no cached payload, no groundtruth
required) on a directory shaped like the El Agente Pathways "plain"
mode output:

    <step_dir>/
      source/reactant*.xyz       # one (or more) reactant xyz
      source/product*.xyz        # product xyz
      initial_ts_guesses/        # 20 LLM-generated IG-TS xyz files
      generation_report.json     # used to read charge, multiplicity

Pipeline:
  1. xtb GFN2 single-point on R, P -> WBO matrices
  2. PQ alignment R <-> P -> atom mapping + broken / formed bonds
  3. core_atoms in R-frame; bond-reaction vector V at TS coords
  4. xtb hess on each IG (parallel) -> g98.out -> normal modes
  5. align each IG <-> R (every branch), reindex modes to R-frame,
     keep the branch with the highest ranker score
  6. per-mode features: bond_overlap (beta), rxn_overlap (rho),
     core_fraction (kappa); pick max-bond_overlap imag mode per IG
  7. score = beta * (1 + w_r * rho) * (1 + w_c * kappa) / n_imag^p,
     sort descending. Every IG with >=1 imaginary mode gets animated
     on its picked mode; IGs with n_imag == 0 render as static.
     No n_imag-count or rho filter on the viewer side.

Output: every artifact for one step lives under a single folder so the
step is fully self-contained (you can move / share / delete a step
without touching anything else):

  out/ranked_views/<workflow_name>/
    view.html              the ranked viewer (open in any browser)
    alignment.json         mapping_RP, broken / formed bonds, core atoms,
                           score-formula constants
    scores.csv             per-IG ranked scores + features
    aligned/
      R.xyz                R geometry (R-frame)
      P_in_R_frame.xyz     P reindexed to R atom order
      iter<N>_in_R_frame.xyz  IG reindexed to R atom order
    modes/
      <label>_picked.xyz   picked imaginary mode in extended xyz format
                           (element x y z dx dy dz)
    xtb/
      R/, P/               xtb single-point output (xyz, wbo, ...)
      hess_iter<N>/        xtb --hess output (g98.out, hessian, ...)

The xtb subtree is the cache: re-running the script on the same step
hits the cache and rebuilds the artifacts in seconds.

Usage:
  python viewer/build_ranked_view_external.py <step_dir> [step_dir ...]

  e.g.
  python viewer/build_ranked_view_external.py \\
      /Users/.../run_20260508_202120_plain_parallel/{1a,1b,1c}

Output:
  out/ranked_views/<workflow_name>.html
"""
from __future__ import annotations
import sys as _sys
from pathlib import Path as _Path
_HERE = _Path(__file__).resolve().parent
_sys.path.insert(0, str(_HERE.parent / "src"))
_sys.path.insert(0, str(_HERE))
PROJECT_ROOT = _HERE.parent

import json
import os
import re
import shutil
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

from rxn_core_pq import align_from_arrays
from rxn_core_frag import run_xtb, parse_xyz
from analyze_core_modes import (
    parse_g98_modes, core_atoms_in_R_frame, reindex_modes_to_R,
    bond_reaction_vector, bond_overlap_per_mode,
    rxn_overlap_per_mode, reaction_coord_delta,
)
from align_bgcp_coords import fill_unmapped_greedy


OUT_DIR = PROJECT_ROOT / "out" / "ranked_views"

# Score-formula hyperparameters (same shape as rk_clean_v2 score, but
# applied without the n_imag<=2 / rho>=0.10 filters — the viewer shows
# every IG that has any imaginary mode, ranked by score).
W_RXN, W_CORE, IMAG_PEN = 1.0, 0.2, 0.3

# xtb parallelism: each worker uses OMP_NUM_THREADS, and we run N_WORKERS
# in parallel. 4 workers x 4 threads = 16, fine on a 14-core box.
N_WORKERS = 4
OMP_THREADS = 4


def parse_xyz_file(path: Path):
    el, xyz = parse_xyz(path)
    return el, xyz


def read_first_xyz(d: Path):
    files = sorted(d.glob("*.xyz"))
    if not files:
        raise FileNotFoundError(f"no xyz in {d}")
    return parse_xyz_file(files[0])


def run_xtb_hess(xyz_path: Path, workdir: Path, charge: int = 0, uhf: int = 0):
    """Run `xtb input.xyz --gfn 2 --hess`. Cached: skip if g98.out + wbo
    are already there with a matching xyz copy. Returns (elements,
    coords, wbo, freqs, modes_TS) where freqs is (n_modes,) and
    modes_TS is (n_modes, n_atoms, 3)."""
    workdir.mkdir(parents=True, exist_ok=True)
    local = workdir / Path(xyz_path).name
    src_text = Path(xyz_path).read_text()
    cached_text = local.read_text() if local.exists() else None
    g98 = workdir / "g98.out"
    wbo = workdir / "wbo"
    cached = (cached_text == src_text) and g98.exists() and wbo.exists()
    if not cached:
        shutil.copy(xyz_path, local)
        cmd = ["xtb", local.name, "--gfn", "2", "--hess"]
        if charge: cmd += ["--chrg", str(charge)]
        if uhf:    cmd += ["--uhf",  str(uhf)]
        env = os.environ.copy()
        env["OMP_NUM_THREADS"] = str(OMP_THREADS)
        res = subprocess.run(cmd, cwd=workdir, capture_output=True,
                             text=True, env=env)
        if res.returncode != 0:
            raise RuntimeError(f"xtb hess failed in {workdir}: "
                               f"{res.stderr[-500:]}")
        if not g98.exists() or not wbo.exists():
            raise RuntimeError(f"missing g98.out or wbo in {workdir}")
    elements, coords = parse_xyz(local)
    n = len(elements)
    wbo_arr = np.zeros((n, n))
    for ln in wbo.read_text().splitlines():
        parts = ln.split()
        if len(parts) < 3: continue
        i, j = int(parts[0]) - 1, int(parts[1]) - 1
        v = float(parts[2])
        wbo_arr[i, j] = v; wbo_arr[j, i] = v
    freqs, modes_TS = parse_g98_modes(g98)
    return elements, coords, wbo_arr, freqs, modes_TS


def _hess_worker(args):
    """Standalone for ProcessPoolExecutor."""
    xyz_path, workdir, charge, uhf = args
    t0 = time.time()
    try:
        elT, xyzT, wboT, freqs, modes_TS = run_xtb_hess(
            Path(xyz_path), Path(workdir), charge, uhf)
        return {"ok": True, "xyz_path": str(xyz_path),
                "secs": time.time() - t0,
                "elT": elT, "xyzT": xyzT.tolist() if hasattr(xyzT, 'tolist') else xyzT,
                "wboT": wboT.tolist(), "freqs": freqs.tolist(),
                "modes_TS": modes_TS.tolist()}
    except Exception as e:
        return {"ok": False, "xyz_path": str(xyz_path),
                "secs": time.time() - t0, "error": str(e)}


def label_for_ig(path: Path):
    m = re.search(r"_iter(\d+)_", path.name)
    return f"iter{m.group(1)}" if m else path.stem


def process_step(step_dir: Path):
    rep = json.loads((step_dir / "generation_report.json").read_text())
    spec = rep["generation_spec"]
    workflow_name = spec["workflow_name"]
    charge = int(spec.get("charge", 0))
    mult = int(spec.get("multiplicity", 1))
    uhf = max(0, mult - 1)
    print(f"\n=== {workflow_name}  charge={charge} mult={mult} ===", flush=True)

    # All artifacts (xtb cache, aligned coords, modes, viewer, JSON, CSV)
    # live under one self-contained directory.
    run_dir = OUT_DIR / workflow_name
    cache = run_dir / "xtb"
    cache.mkdir(parents=True, exist_ok=True)

    # 1. R, P single-points
    rxyz_path = sorted((step_dir / "source").glob("reactant*.xyz"))[0]
    pxyz_path = sorted((step_dir / "source").glob("product*.xyz"))[0]
    print(f"  xtb sp on R, P ...", flush=True)
    t0 = time.time()
    elR, xyzR, wboR = run_xtb(rxyz_path, cache / "R", charge=charge, uhf=uhf)
    elP, xyzP, wboP = run_xtb(pxyz_path, cache / "P", charge=charge, uhf=uhf)
    print(f"    sp done in {time.time()-t0:.1f}s ({len(elR)} atoms each)", flush=True)

    # 2. R<->P alignment + core atoms + bond-reaction vector at R-coords
    print(f"  PQ alignment R<->P ...", flush=True)
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
    print(f"    mapping={len(mapping_RP)}/{len(elR)}, broken={len(broken_R)}, "
          f"formed={len(formed_R)}, core_atoms={len(core_R)}", flush=True)

    # 3. xtb hess on each IG, in parallel
    ig_paths = sorted((step_dir / "initial_ts_guesses").glob("*.xyz"))
    print(f"  xtb hess on {len(ig_paths)} IGs (workers={N_WORKERS}, omp={OMP_THREADS}) ...",
          flush=True)
    t0 = time.time()
    args = [(str(p), str(cache / f"hess_{label_for_ig(p)}"), charge, uhf)
            for p in ig_paths]
    results = []
    with ProcessPoolExecutor(max_workers=N_WORKERS) as ex:
        futures = {ex.submit(_hess_worker, a): a for a in args}
        for done in as_completed(futures):
            r = done.result()
            results.append(r)
            tag = "ok" if r["ok"] else f"FAIL({r.get('error','?')[:60]})"
            print(f"    [{r['secs']:>5.1f}s]  {Path(r['xyz_path']).name}  {tag}",
                  flush=True)
    print(f"  all hess done in {time.time()-t0:.1f}s", flush=True)

    # 4. Per-IG: align IG<->R, reindex modes, compute features, score
    n_R = len(elR)
    bond_V = None  # built per-IG using IG coords in R-frame
    ig_records = []
    for r in sorted(results, key=lambda r: r["xyz_path"]):
        label = label_for_ig(Path(r["xyz_path"]))
        if not r["ok"]:
            ig_records.append({
                "label": label, "score": 0.0,
                "beta": 0.0, "rho": 0.0, "kappa": 0.0, "n_imag": 0,
                "picked_freq": None, "picked_disp": None,
                "xyz_elements": elR, "xyz_coords": xyzR.tolist() if hasattr(xyzR,'tolist') else xyzR,
                "error": r.get("error", ""),
            })
            continue
        elT = r["elT"]
        xyzT = np.asarray(r["xyzT"], float)
        wboT = np.asarray(r["wboT"], float)
        freqs = np.asarray(r["freqs"], float)
        modes_TS = np.asarray(r["modes_TS"], float)

        # Align IG<->R, asking for ALL equally-good branches so we can
        # pick the mapping that maximizes the ranker score (rather than
        # using whichever branch happened to sort first under the
        # alignment-only score, which doesn't see chemistry).
        try:
            it = align_from_arrays(elR, xyzR, wboR, elT, xyzT, wboT,
                                   return_all=True)
        except Exception as e:
            print(f"    align fail for {label}: {e}")
            ig_records.append({
                "label": label, "score": 0.0,
                "beta": 0.0, "rho": 0.0, "kappa": 0.0, "n_imag": 0,
                "picked_freq": None, "picked_disp": None,
                "xyz_elements": elT, "xyz_coords": xyzT.tolist(),
                "error": f"align: {e}",
            })
            continue

        # Precompute the imag-mode list once (independent of mapping).
        imag_idx = list(np.where(freqs < 0)[0])
        n_imag = len(imag_idx)
        if n_imag == 0:
            ig_records.append({
                "label": label, "score": 0.0,
                "beta": 0.0, "rho": 0.0, "kappa": 0.0, "n_imag": 0,
                "picked_freq": None, "picked_disp": None,
                "xyz_elements": elT, "xyz_coords": xyzT.tolist(),
            })
            continue

        # Try every alignment branch; keep the one whose picked-mode
        # score is highest. Two-level dedup:
        #   1. full mapping dedup -- collapse exact duplicates from the
        #      multi-seed sweep
        #   2. core-only mapping dedup -- branches that differ only in
        #      non-core (spectator) atom assignments give identical
        #      beta/rho/kappa, so they're score-equivalent and we only
        #      need to evaluate one. This is the score-relevant set.
        all_branches = it.get("all_scored", [])
        if not all_branches:
            all_branches = [(None, dict(it["mapping"]), None, None, None)]
        seen_full = set()
        seen_core = set()
        best = None  # (score, b, r, c, picked_k, modes_R, ts_xyz_in_R)
        for (_, br_mapping, _, _, _) in all_branches:
            br_d = dict(br_mapping)
            full_key = tuple(sorted(br_d.items()))
            if full_key in seen_full:
                continue
            seen_full.add(full_key)

            # Score-equivalent core-only signature
            core_key = tuple(sorted(
                (c, br_d[c]) for c in core_R if c in br_d
            ))
            if core_key in seen_core:
                continue
            seen_core.add(core_key)

            mapping_RT = fill_unmapped_greedy(elR, xyzR, elT, xyzT, br_d)
            modes_R = reindex_modes_to_R(modes_TS, mapping_RT, n_R)
            sq = (modes_R ** 2).sum(axis=2)
            total = sq.sum(axis=1)
            core_e = (sq[:, core_R].sum(axis=1)
                      if core_R else np.zeros(modes_R.shape[0]))
            kappa = np.where(total > 1e-12, core_e / total, 0.0)
            rho = rxn_overlap_per_mode(modes_R, delta_RP, core_R)
            ts_xyz_in_R = np.zeros_like(np.asarray(xyzR))
            for r_idx, t_idx in mapping_RT.items():
                ts_xyz_in_R[r_idx] = xyzT[t_idx]
            V = bond_reaction_vector(ts_xyz_in_R, broken_R, formed_R)
            beta = bond_overlap_per_mode(modes_R, V)

            picked_k = max(imag_idx, key=lambda k: beta[k])
            b  = float(beta[picked_k])
            r_ = float(rho[picked_k])
            c  = float(kappa[picked_k])
            score = (b * (1 + W_RXN * r_) * (1 + W_CORE * c)
                     / max(n_imag, 1) ** IMAG_PEN)
            if best is None or score > best[0]:
                best = (score, b, r_, c, picked_k, modes_R, ts_xyz_in_R)

        score, b, r_, c, picked_k, modes_R, ts_xyz_in_R = best
        if len(seen_full) > 1:
            collapsed = (f" (collapsed from {len(seen_full)})"
                         if len(seen_core) < len(seen_full) else "")
            print(f"    {label}: {len(seen_core)} core-unique branches"
                  f"{collapsed}; best score={score:.3f}",
                  flush=True)
        ig_records.append({
            "label": label,
            "score": float(score),
            "beta": b, "rho": r_, "kappa": c,
            "n_imag": int(n_imag),
            "picked_freq": float(freqs[picked_k]),
            "picked_disp": modes_R[picked_k].tolist(),
            "xyz_elements": elR,
            "xyz_coords":   ts_xyz_in_R.tolist(),
            "n_branches": len(seen_core),
            "n_branches_full": len(seen_full),
        })

    # 5. Sort by score descending. IGs with no imag mode (score=0) sink
    # to the bottom naturally.
    ig_records.sort(key=lambda x: -x["score"])

    # 6. Reindex P xyz to R-frame so bond pairs (which live in R-frame)
    # draw on the right atoms in the P panel. full_RP maps every R-index
    # to a P-index (alignment + greedy fill).
    xyzP_arr = np.asarray(xyzP, float)
    xyzP_in_R = np.zeros_like(np.asarray(xyzR, float))
    for i_R, i_P in full_RP.items():
        xyzP_in_R[i_R] = xyzP_arr[i_P]

    data = {
        "step": workflow_name,
        "n_atoms": n_R,
        "core_atoms": list(core_R),
        "broken_bonds":  [list(b) for b in broken_R],
        "formed_bonds_R":[list(b) for b in formed_R],
        "reactant":  {"xyz_elements": elR,
                      "xyz_coords":   xyzR.tolist() if hasattr(xyzR,'tolist') else xyzR},
        # P shown in R-frame so atom indices line up with broken/formed pairs
        "product":   {"xyz_elements": elR,
                      "xyz_coords":   xyzP_in_R.tolist()},
        "igs": ig_records,
    }
    # 7. Write all artifacts under run_dir (one self-contained folder).
    out_path = run_dir / "view.html"
    out_path.write_text(HTML.format(
        title=f"Ranked view — {workflow_name}",
        n_atoms=n_R,
        n_broken=len(broken_R),
        n_formed=len(formed_R),
        n_core=len(core_R),
        data_json=json.dumps(data),
    ))
    write_artifacts(
        run_dir=run_dir,
        workflow_name=workflow_name, charge=charge, mult=mult,
        elR=elR, xyzR=np.asarray(xyzR, float),
        xyzP_in_R=xyzP_in_R,
        mapping_RP=mapping_RP,
        broken_R=broken_R, formed_R=formed_R, core_R=core_R,
        ig_records=ig_records,
    )

    has_mode = sum(1 for ig in ig_records if ig["n_imag"] > 0)
    print(f"  IGs with at least one imag mode: {has_mode}/{len(ig_records)}", flush=True)
    print(f"  top-3 by score:")
    for ig in ig_records[:3]:
        print(f"    {ig['label']:>8s}  S={ig['score']:.3f}  "
              f"beta={ig['beta']:.3f}  rho={ig['rho']:.3f}  "
              f"kappa={ig['kappa']:.3f}  n_imag={ig['n_imag']}", flush=True)
    print(f"  wrote {run_dir}/", flush=True)
    return run_dir


def _xyz_block(elements, xyz, comment=""):
    n = len(elements)
    lines = [str(n), comment]
    for el, c in zip(elements, xyz):
        lines.append(f"{el:<3s}  {c[0]:14.8f}  {c[1]:14.8f}  {c[2]:14.8f}")
    return "\n".join(lines) + "\n"


def _xyz_with_disp(elements, xyz, disp, comment=""):
    """Extended xyz: element x y z dx dy dz."""
    n = len(elements)
    lines = [str(n), comment]
    for el, c, d in zip(elements, xyz, disp):
        lines.append(f"{el:<3s}  {c[0]:14.8f}  {c[1]:14.8f}  {c[2]:14.8f}"
                     f"   {d[0]:10.6f}  {d[1]:10.6f}  {d[2]:10.6f}")
    return "\n".join(lines) + "\n"


def write_artifacts(run_dir, workflow_name, charge, mult,
                    elR, xyzR, xyzP_in_R, mapping_RP,
                    broken_R, formed_R, core_R, ig_records):
    aligned = run_dir / "aligned"
    modes_d = run_dir / "modes"
    aligned.mkdir(exist_ok=True)
    modes_d.mkdir(exist_ok=True)

    (aligned / "R.xyz").write_text(
        _xyz_block(elR, xyzR.tolist(),
                   comment=f"R for {workflow_name} (R-frame)"))
    (aligned / "P_in_R_frame.xyz").write_text(
        _xyz_block(elR, xyzP_in_R.tolist(),
                   comment=f"P (reindexed to R-frame) for {workflow_name}"))

    for ig in ig_records:
        label = ig["label"]
        (aligned / f"{label}_in_R_frame.xyz").write_text(
            _xyz_block(elR, ig["xyz_coords"],
                       comment=f"{label} (reindexed to R-frame) for {workflow_name}"))
        if ig.get("picked_disp") is not None:
            freq = ig["picked_freq"]
            comment = (f"{label} picked imag mode  "
                       f"freq={freq:.1f}  "
                       f"beta={ig['beta']:.3f}  rho={ig['rho']:.3f}  "
                       f"kappa={ig['kappa']:.3f}  score={ig['score']:.3f}")
            (modes_d / f"{label}_picked.xyz").write_text(
                _xyz_with_disp(elR, ig["xyz_coords"], ig["picked_disp"],
                               comment=comment))

    alignment = {
        "step": workflow_name,
        "n_atoms": len(elR),
        "charge": charge,
        "multiplicity": mult,
        "elements_R": elR,
        "broken_bonds_R":   [list(b) for b in broken_R],
        "formed_bonds_R":   [list(b) for b in formed_R],
        "core_atoms_R":     list(core_R),
        "mapping_RP":       {str(k): int(v) for k, v in mapping_RP.items()},
        "score_formula": {
            "form":  "S = beta * (1 + w_r * rho) * (1 + w_c * kappa) / n_imag^p",
            "w_r": W_RXN, "w_c": W_CORE, "p": IMAG_PEN,
            "filter": "none (every IG with >=1 imag mode is ranked)",
        },
    }
    (run_dir / "alignment.json").write_text(json.dumps(alignment, indent=2))

    csv_lines = ["rank,label,score,beta,rho,kappa,n_imag,picked_freq,"
                 "n_branches_core,n_branches_full"]
    for rank, ig in enumerate(ig_records, 1):
        freq_str = (f"{ig['picked_freq']:.2f}"
                    if ig.get('picked_freq') is not None else "")
        nb_core = ig.get('n_branches', 1)
        nb_full = ig.get('n_branches_full', nb_core)
        csv_lines.append(
            f"{rank},{ig['label']},{ig['score']:.4f},{ig['beta']:.4f},"
            f"{ig['rho']:.4f},{ig['kappa']:.4f},{ig['n_imag']},"
            f"{freq_str},{nb_core},{nb_full}"
        )
    (run_dir / "scores.csv").write_text("\n".join(csv_lines) + "\n")

    (run_dir / "README.md").write_text(f"""# {workflow_name}

End-to-end ranked-view artifacts for one elementary reaction step.

## Layout

```
view.html              open in any browser; shows R + P + 20 sorted IG panels
alignment.json         mapping_RP, broken/formed/core (R-frame), score-formula constants
scores.csv             per-IG ranked features (rank, score, beta, rho, kappa, ...)
aligned/
  R.xyz                reactant geometry (R-frame indexing)
  P_in_R_frame.xyz     product reindexed to R atom order
  iter<N>_in_R_frame.xyz   IG reindexed to R atom order
modes/
  <label>_picked.xyz   picked imaginary mode in extended xyz
                       (element x y z dx dy dz)
xtb/
  R/, P/               xtb GFN2 single-point output (xyz, wbo, ...)
  hess_iter<N>/        xtb --hess output (g98.out, hessian, vibspectrum, ...)
```

## Score formula

S = beta * (1 + 1.0 * rho) * (1 + 0.2 * kappa) / n_imag^0.3

- beta = picked-mode displacement projected on broken/formed bond axes
- rho  = picked-mode core-atom motion projected on R->P direction
- kappa = fraction of picked-mode energy localized on the reactive core

IGs are ranked by S descending; every IG with at least one imaginary
mode is animated on its picked imag mode (max-beta among imag modes).
For each IG we sweep every alignment branch (when the IG<->R mapping
is set-non-unique) and keep the branch with the highest S, so a
spectator-symmetry tie can't artificially deflate a good IG.

## Reproducing

```
python viewer/build_ranked_view_external.py <step_dir>
```

re-runs the pipeline; the xtb subtree is the cache, so a re-run on the
same step is a sub-second rebuild of the artifacts.
""")


HTML = r"""<!doctype html>
<html><head><meta charset="utf-8">
<title>{title}</title>
<script src="https://3dmol.org/build/3Dmol-min.js"></script>
<style>
 html, body {{ margin:0; padding:0; }}
 body {{ font-family:-apple-system,sans-serif; background:#fafafa; padding:14px; box-sizing:border-box; }}
 h2 {{ margin:0 0 4px; font-size:18px; }}
 .sub {{ color:#444; font-size:13px; margin-bottom:14px; }}
 .legend span {{ display:inline-block; padding:2px 8px; border-radius:4px; font-size:12px; margin-right:6px; }}
 .leg-broken {{ background:#fcd3d3; color:#a00; }}
 .leg-formed {{ background:#cdebd0; color:#070; }}
 .leg-mode   {{ background:#d6e7ff; color:#024; }}
 .leg-static {{ background:#eee;    color:#666; }}
 .ref-row {{ display:grid; grid-template-columns:repeat(2,1fr); gap:12px; margin-bottom:18px; }}
 .ig-grid {{ display:grid; grid-template-columns:repeat(4,1fr); gap:12px; }}
 .panel  {{ background:white; border:1px solid #ddd; border-radius:6px; padding:6px 8px 8px; }}
 .panel.no-imag {{ background:#f7f7f7; }}
 .ph    {{ display:flex; justify-content:space-between; align-items:baseline; font-size:12px; margin-bottom:4px; }}
 .ph .lbl {{ font-weight:600; font-size:13px; }}
 .ph .rk  {{ font-family:ui-monospace,monospace; color:#024; }}
 .panel.no-imag .ph .rk {{ color:#888; }}
 .vw    {{ position:relative; width:100%; height:230px; }}
 .ref-row .vw {{ height:300px; }}
 .vwbox {{ position:absolute; inset:0; }}
 .meta  {{ font-family:ui-monospace,monospace; font-size:11px; color:#444; padding:3px 0 0; line-height:1.4; }}
 .meta b {{ color:#024; }}
 .badge-static {{ display:inline-block; background:#eee; color:#666; padding:1px 6px; border-radius:3px; font-size:11px; margin-left:6px; }}
</style>
</head>
<body>
<h2>{title}</h2>
<div class="sub">
  N atoms: <b>{n_atoms}</b> &nbsp;|&nbsp;
  broken: <b>{n_broken}</b> &nbsp;|&nbsp;
  formed: <b>{n_formed}</b> &nbsp;|&nbsp;
  core atoms: <b>{n_core}</b><br>
  <span class="legend">
    <span class="leg-broken">broken bond</span>
    <span class="leg-formed">formed bond</span>
    <span class="leg-mode">mode arrow (core atoms)</span>
    <span class="leg-static">n_imag = 0 (rendered static)</span>
  </span>
  <br>
  IGs are sorted by score
  <code>S = &beta; (1 + w_r &rho;) (1 + w_c &kappa;) / n_imag^p</code>
  with <code>w_r=1.0, w_c=0.2, p=0.3</code>, descending. Every IG with
  &ge;1 imaginary mode is animated on its picked mode (max-&beta; imag);
  IGs with <code>n_imag = 0</code> render as static structures.
</div>

<div class="ref-row">
  <div class="panel"><div class="ph"><span class="lbl">Reactant</span><span class="rk">static</span></div>
    <div class="vw"><div id="vw_R" class="vwbox"></div></div></div>
  <div class="panel"><div class="ph"><span class="lbl">Product</span><span class="rk">static</span></div>
    <div class="vw"><div id="vw_P" class="vwbox"></div></div></div>
</div>

<div class="ig-grid" id="grid"></div>

<script>
const DATA = {data_json};

function buildBody(elements, xyz) {{
  const n = xyz.length;
  let s = `${{n}}\nframe\n`;
  for (let i = 0; i < n; i++) {{
    s += `${{elements[i]}}  ${{xyz[i][0].toFixed(6)}}  ${{xyz[i][1].toFixed(6)}}  ${{xyz[i][2].toFixed(6)}}\n`;
  }}
  return s;
}}
function buildBodyAt(elements, xyz, disp, scale) {{
  const n = xyz.length;
  let s = `${{n}}\nframe\n`;
  for (let i = 0; i < n; i++) {{
    const x = xyz[i][0] + scale * disp[i][0];
    const y = xyz[i][1] + scale * disp[i][1];
    const z = xyz[i][2] + scale * disp[i][2];
    s += `${{elements[i]}}  ${{x.toFixed(6)}}  ${{y.toFixed(6)}}  ${{z.toFixed(6)}}\n`;
  }}
  return s;
}}
function xyzAt(xyz, disp, scale) {{
  const out = new Array(xyz.length);
  for (let i = 0; i < xyz.length; i++) {{
    out[i] = [
      xyz[i][0] + scale * disp[i][0],
      xyz[i][1] + scale * disp[i][1],
      xyz[i][2] + scale * disp[i][2],
    ];
  }}
  return out;
}}
// which: 'broken' (R panel only) | 'formed' (P panel only) | 'both' (TS-like)
// Solid-color dashed cylinder, radius 0.20 — wider than the default
// 3Dmol stick-bond radius (0.10) so the marker stays visible on R / P
// where 3Dmol auto-draws a stick bond on top of the same atom pair.
// Uses xyz coords directly so cylinders can be added before render().
function decorateEvent(viewer, xyz, pair, color) {{
  const [i, j] = pair;
  if (i >= xyz.length || j >= xyz.length) return;
  const a = xyz[i], b = xyz[j];
  viewer.addCylinder({{
    start:{{x:a[0], y:a[1], z:a[2]}}, end:{{x:b[0], y:b[1], z:b[2]}},
    color:color, radius:0.20, dashed:true,
  }});
}}
function decorateBonds(viewer, which, xyz) {{
  if (which !== 'formed') {{
    for (const pair of DATA.broken_bonds) decorateEvent(viewer, xyz, pair, 'red');
  }}
  if (which !== 'broken') {{
    for (const pair of DATA.formed_bonds_R) decorateEvent(viewer, xyz, pair, 'green');
  }}
}}
function drawArrows(viewer, xyz, disp) {{
  for (const i of DATA.core_atoms) {{
    if (i >= xyz.length || !disp || !disp[i]) continue;
    const d = disp[i];
    const len = Math.hypot(d[0], d[1], d[2]);
    if (len < 1e-3) continue;
    const ax = xyz[i][0], ay = xyz[i][1], az = xyz[i][2];
    viewer.addArrow({{
      start:{{x:ax, y:ay, z:az}},
      end:  {{x:ax + d[0]*1.5, y:ay + d[1]*1.5, z:az + d[2]*1.5}},
      color:'#0066cc', radius:0.06,
    }});
  }}
}}
function makeStatic(divId, ts, which) {{
  const v = $3Dmol.createViewer(divId, {{backgroundColor:'white'}});
  v.addModel(buildBody(ts.xyz_elements, ts.xyz_coords), 'xyz');
  v.setStyle({{}}, {{stick:{{radius:0.10}}, sphere:{{scale:0.20}}}});
  decorateBonds(v, which || 'both', ts.xyz_coords);
  v.zoomTo();
  v.render();
  return v;
}}
function makeAnimated(divId, ts, disp, which) {{
  const w = which || 'both';
  const v = $3Dmol.createViewer(divId, {{backgroundColor:'white'}});
  v.addModel(buildBody(ts.xyz_elements, ts.xyz_coords), 'xyz');
  v.setStyle({{}}, {{stick:{{radius:0.10}}, sphere:{{scale:0.20}}}});
  decorateBonds(v, w, ts.xyz_coords);
  drawArrows(v, ts.xyz_coords, disp);
  v.zoomTo(); v.render();
  let t = 0; const period = 30; const amp = 0.6;
  setInterval(() => {{
    t = (t + 1) % period;
    const scale = amp * Math.sin(2 * Math.PI * t / period);
    const cur = xyzAt(ts.xyz_coords, disp, scale);
    v.removeAllModels(); v.removeAllShapes();
    v.addModel(buildBodyAt(ts.xyz_elements, ts.xyz_coords, disp, scale), 'xyz');
    v.setStyle({{}}, {{stick:{{radius:0.10}}, sphere:{{scale:0.20}}}});
    decorateBonds(v, w, cur);
    drawArrows(v, cur, disp);
    v.render();
  }}, 60);
  return v;
}}

window.addEventListener('load', () => {{
  // R panel: only the bonds that are BREAKING (red dashed) -- they
  //          exist in R right now and are about to disappear.
  // P panel: only the bonds that just FORMED (green dashed) -- they
  //          exist in P. P xyz is reindexed to R-frame above so the
  //          formed_bonds_R pairs land on the right atoms.
  makeStatic('vw_R', DATA.reactant, 'broken');
  makeStatic('vw_P', DATA.product,  'formed');
  const grid = document.getElementById('grid');
  DATA.igs.forEach((ig, i) => {{
    const div = document.createElement('div');
    const hasMode = !!ig.picked_disp;
    div.className = 'panel' + (hasMode ? '' : ' no-imag');
    const staticBadge = hasMode ? '' :
      `<span class="badge-static">n_imag = 0</span>`;
    const freqStr = ig.picked_freq != null ? ig.picked_freq.toFixed(0) + 'i' : '—';
    div.innerHTML = `
      <div class="ph">
        <span class="lbl">${{ig.label}}${{staticBadge}}</span>
        <span class="rk">S = ${{ig.score.toFixed(3)}}</span>
      </div>
      <div class="vw"><div id="vw_ig${{i}}" class="vwbox"></div></div>
      <div class="meta">
        <b>&beta;</b>=${{ig.beta.toFixed(3)}} &nbsp;
        <b>&rho;</b>=${{ig.rho.toFixed(3)}} &nbsp;
        <b>&kappa;</b>=${{ig.kappa.toFixed(3)}} &nbsp;
        <b>n_imag</b>=${{ig.n_imag}} &nbsp;
        <b>freq</b>=${{freqStr}} cm⁻¹
      </div>`;
    grid.appendChild(div);
    if (ig.picked_disp) {{
      makeAnimated(`vw_ig${{i}}`, ig, ig.picked_disp);
    }} else {{
      makeStatic(`vw_ig${{i}}`, ig);
    }}
  }});
}});
</script>
</body></html>
"""


def main():
    if len(sys.argv) < 2:
        print(__doc__); sys.exit(1)
    for arg in sys.argv[1:]:
        step_dir = Path(arg)
        if not step_dir.is_dir():
            print(f"  WARN: skip {step_dir} (not a directory)")
            continue
        try:
            process_step(step_dir)
        except Exception as e:
            import traceback
            traceback.print_exc()
            print(f"  ERROR processing {step_dir}: {e}")


if __name__ == "__main__":
    main()
