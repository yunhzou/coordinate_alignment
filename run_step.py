"""
Driver: pick a benchmark step, run analysis, write HTML.

Usage:
  python run_step.py <step_name>          # e.g. pr1.tempo_ts1
  python run_step.py <step_name> [--charge N] [--uhf N] [--out PATH]
"""

from __future__ import annotations
import argparse
from pathlib import Path
import sys

from rxn_core_wbo import analyze
from viz_3dmol import render_html


BENCHMARK_ROOT = Path("/Users/yunhengz/empty_for_claude/Benchmark")


def find_step(name):
    p = BENCHMARK_ROOT / name / "plain" / "stage0"
    r = p / "reactant.xyz"
    pr = p / "product.xyz"
    if not (r.exists() and pr.exists()):
        raise SystemExit(f"missing reactant.xyz or product.xyz under {p}")
    return r, pr


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("step")
    ap.add_argument("--charge", type=int, default=0)
    ap.add_argument("--uhf", type=int, default=0)
    ap.add_argument("--bond-cut", type=float, default=0.5)
    ap.add_argument("--max-radius", type=int, default=4)
    ap.add_argument("--out", default=None)
    ap.add_argument("--show-map", action="store_true")
    args = ap.parse_args()

    rxn, prd = find_step(args.step)
    workdir = Path(__file__).parent / "work" / args.step
    print(f"[run] reactant={rxn}\n[run] product={prd}")
    result = analyze(
        rxn, prd, workdir,
        charge=args.charge, uhf=args.uhf,
        max_radius=args.max_radius, bond_cut=args.bond_cut,
    )

    print(f"[result] anchors found: {len(result['anchors'])}")
    if result["anchors"]:
        rmin = min(a[2] for a in result["anchors"])
        rmax = max(a[2] for a in result["anchors"])
        print(f"          anchor radii span: {rmin}..{rmax}")
    print(f"[result] mapping size: {len(result['mapping'])} / {len(result['elements_R'])}"
          f" (spectator phase: {result['n_spectator']})")
    if args.show_map:
        for i in sorted(result["mapping"]):
            print(f"           R[{i:>2}]({result['elements_R'][i]}) -> "
                  f"P[{result['mapping'][i]:>2}]({result['elements_P'][result['mapping'][i]]})")
    print(f"[result] broken bonds: {len(result['broken'])}")
    for (i, j, wR, wP) in result["broken"]:
        wP_s = "—" if wP is None else f"{wP:.2f}"
        print(f"           R[{i}]-R[{j}]  WBO_R={wR:.2f}  WBO_P={wP_s}  "
              f"({result['elements_R'][i]}-{result['elements_R'][j]})")
    print(f"[result] formed bonds: {len(result['formed'])}")
    for (i, j, wR, wP) in result["formed"]:
        wR_s = "—" if wR is None else f"{wR:.2f}"
        print(f"           P[{i}]-P[{j}]  WBO_R={wR_s}  WBO_P={wP:.2f}  "
              f"({result['elements_P'][i]}-{result['elements_P'][j]})")
    print(f"[result] core atoms (R): {result['core_R']}")
    print(f"[result] core atoms (P): {result['core_P']}")

    out = args.out or str(Path(__file__).parent / "out" / f"{args.step}.html")
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    render_html(result, title=args.step, out_path=out)
    print(f"[viz] wrote {out}")


if __name__ == "__main__":
    main()
