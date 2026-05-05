"""
Same as run_pq_bgcp.py, but using rxn_core_frag.analyze (the old algorithm).
For apples-to-apples comparison against bgcp_pq_bonds.csv.
Output: out/bgcp_old_bonds.csv
"""
from __future__ import annotations
import argparse
import csv
import re
import sys
import time
import traceback
from pathlib import Path

from rxn_core_frag import analyze
from build_bgcp_viewer import BGCP_ROOT, WORK, LOOKUP, list_step_dirs, read_xyzs

OUT_CSV = Path(__file__).parent / "out" / "bgcp_old_bonds.csv"


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
    out = analyze(wd / "reactant.xyz", wd / "product.xyz", wd,
                  charge=chg, uhf=uhf)
    return out, chg, uhf, time.time() - t0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--out", default=str(OUT_CSV))
    args = ap.parse_args()
    steps = list_step_dirs()[:args.limit] if args.limit else list_step_dirs()
    print(f"Running OLD algorithm on {len(steps)} BGCP steps")
    rows = []; t_total = 0; n_done = n_err = 0
    for i, sd in enumerate(steps, 1):
        name = sd.name
        try:
            out, chg, uhf, elapsed = run_one(sd)
            t_total += elapsed
            row = {
                "step_id": name,
                "n_atoms": len(out["elements_R"]),
                "n_mapped": len(out["mapping"]),
                "n_broken": len(out["broken"]),
                "n_formed": len(out["formed"]),
                "charge": chg, "multiplicity": uhf + 1,
                "broken_bonds": fmt_bonds(out["broken"], out["elements_R"]),
                "formed_bonds": fmt_bonds(out["formed"], out["elements_P"]),
                "time_s": f"{elapsed:.2f}", "error": "",
            }
            rows.append(row); n_done += 1
            print(f"[{i:3d}/{len(steps)}] {name[:60]:60s} "
                  f"N={row['n_atoms']:>4} m={row['n_mapped']:>4} "
                  f"br/fm={row['n_broken']}/{row['n_formed']:<2} "
                  f"t={elapsed:>5.1f}s")
            sys.stdout.flush()
        except Exception as e:
            n_err += 1
            rows.append({"step_id": name, "n_atoms": "", "n_mapped": "",
                         "n_broken": "", "n_formed": "", "charge": "",
                         "multiplicity": "", "broken_bonds": "", "formed_bonds": "",
                         "time_s": "", "error": str(e)[:200]})
            print(f"[{i:3d}/{len(steps)}] {name[:60]:60s} ERROR: {e}")
            traceback.print_exc(); sys.stdout.flush()
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="") as f:
        fieldnames = ["step_id", "n_atoms", "n_mapped", "n_broken", "n_formed",
                      "charge", "multiplicity", "broken_bonds", "formed_bonds",
                      "time_s", "error"]
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader(); w.writerows(rows)
    print(f"\nDone: {n_done} ok, {n_err} errors. Total {t_total:.1f}s. CSV: {out_path}")


if __name__ == "__main__":
    main()
