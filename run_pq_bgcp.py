"""
Run priority-queue alignment + bond classification on all BGCP steps.

Reuses existing xtb caches at work_bgcp/<sanitized>/{R,P}/wbo when present
(same workdir layout as build_bgcp_viewer.py). For steps without a cache,
writes the inputs and runs xtb fresh — these become caches for next time.

Output:
  out/bgcp_pq_bonds.csv — step_id, n_atoms, n_mapped, n_broken, n_formed,
                          chirality_violations, charge, multiplicity,
                          broken_bonds, formed_bonds, time_s, error

Usage:
  python run_pq_bgcp.py
  python run_pq_bgcp.py --limit 5
  python run_pq_bgcp.py --steps pr1.tempo_ts2 pr12.Co_Silylation_JACS2015_TS_Dstar-Estar
"""
from __future__ import annotations
import argparse
import csv
import re
import sys
import time
import traceback
from pathlib import Path

from rxn_core_pq import analyze_pq
from bgcp_io import (
    BGCP_ROOT, WORK, LOOKUP, list_step_dirs, read_xyzs,
)

OUT_CSV = Path(__file__).parent / "out" / "bgcp_pq_bonds.csv"


def fmt_bonds(bonds, elements):
    return ";".join(f"{i}({elements[i]})-{j}({elements[j]})"
                    for (i, j, _, _) in bonds)


def run_one(step_dir):
    name = step_dir.name
    chg, uhf = LOOKUP.get(name, (0, 0))
    sanitized = re.sub(r"[^A-Za-z0-9._-]", "_", name)
    wd = WORK / sanitized
    wd.mkdir(parents=True, exist_ok=True)
    rxyz = read_xyzs(step_dir / "reactants")
    pxyz = read_xyzs(step_dir / "products")
    if rxyz is None or pxyz is None:
        raise RuntimeError("missing reactant or product")
    (wd / "reactant.xyz").write_text(rxyz)
    (wd / "product.xyz").write_text(pxyz)

    t0 = time.time()
    out = analyze_pq(wd / "reactant.xyz", wd / "product.xyz", wd,
                     charge=chg, uhf=uhf)
    elapsed = time.time() - t0
    return out, chg, uhf, elapsed


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--steps", nargs="+", default=None,
                    help="Specific step dir names to run")
    ap.add_argument("--out", default=str(OUT_CSV))
    args = ap.parse_args()

    all_steps = list_step_dirs()
    if args.steps:
        wanted = set(args.steps)
        steps = [s for s in all_steps if s.name in wanted]
    elif args.limit:
        steps = all_steps[:args.limit]
    else:
        steps = all_steps

    print(f"Running PQ alignment on {len(steps)} BGCP steps")
    print(f"  cache dir: {WORK}")
    print(f"  output:    {args.out}")
    print()

    rows = []
    n_done = n_err = 0
    t_total = 0.0
    for i, sd in enumerate(steps, 1):
        name = sd.name
        try:
            out, chg, uhf, elapsed = run_one(sd)
            t_total += elapsed
            row = {
                "step_id": name,
                "n_atoms": len(out["elements_R"]),
                "n_mapped": out["n_mapped"],
                "n_broken": out["n_broken"],
                "n_formed": out["n_formed"],
                "chirality_violations": out["chirality_violations"],
                "charge": chg,
                "multiplicity": uhf + 1,
                "broken_bonds": fmt_bonds(out["broken"], out["elements_R"]),
                "formed_bonds": fmt_bonds(out["formed"], out["elements_P"]),
                "time_s": f"{elapsed:.2f}",
                "error": "",
            }
            rows.append(row)
            n_done += 1
            print(f"[{i:3d}/{len(steps)}] {name[:60]:60s} "
                  f"N={row['n_atoms']:>4} m={row['n_mapped']:>4} "
                  f"br/fm={row['n_broken']}/{row['n_formed']:<2} "
                  f"chir={row['chirality_violations']:<2} "
                  f"t={elapsed:>5.1f}s")
            sys.stdout.flush()
        except Exception as e:
            n_err += 1
            err_msg = str(e)[:200]
            rows.append({
                "step_id": name,
                "n_atoms": "", "n_mapped": "", "n_broken": "", "n_formed": "",
                "chirality_violations": "",
                "charge": "", "multiplicity": "",
                "broken_bonds": "", "formed_bonds": "",
                "time_s": "",
                "error": err_msg,
            })
            print(f"[{i:3d}/{len(steps)}] {name[:60]:60s} ERROR: {err_msg}")
            traceback.print_exc()
            sys.stdout.flush()

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="") as f:
        fieldnames = ["step_id", "n_atoms", "n_mapped", "n_broken", "n_formed",
                      "chirality_violations", "charge", "multiplicity",
                      "broken_bonds", "formed_bonds", "time_s", "error"]
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)

    print()
    print("=" * 70)
    print(f"Done: {n_done} ok, {n_err} errors")
    print(f"Total wall time: {t_total:.1f}s "
          f"(avg {t_total/max(1, n_done):.1f}s/step)")
    print(f"CSV: {out_path}")


if __name__ == "__main__":
    main()
