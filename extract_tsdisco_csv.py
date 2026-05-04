"""Extract per-step bond-change summary from out/tsdisco_viewer.html into a CSV."""
from __future__ import annotations
import csv
import json
import re
from pathlib import Path

VIEWER = Path("/Users/yunhengz/empty_for_claude/rxn_core/out/tsdisco_viewer.html")
OUT_CSV = Path("/Users/yunhengz/empty_for_claude/rxn_core/out/tsdisco_bonds.csv")


def main():
    html = VIEWER.read_text()
    m = re.search(r"const STEPS = (\{.*?\});\s*\nconst stepNames", html, re.DOTALL)
    if m is None:
        raise SystemExit("Could not find STEPS in viewer HTML")
    steps = json.loads(m.group(1))

    rows = []
    for name, d in steps.items():
        # parse element symbols out of the xyz strings
        def els(xyz):
            lines = xyz.strip().splitlines()
            n = int(lines[0])
            return [ln.split()[0] for ln in lines[2:2 + n]]
        elR = els(d["xyzR"])
        elP = els(d["xyzP"])

        broken_str = ";".join(
            f"{i}({elR[i]})-{j}({elR[j]})" for (i, j) in d["broken_idx"]
        )
        formed_str = ";".join(
            f"{i}({elP[i]})-{j}({elP[j]})" for (i, j) in d["formed_idx"]
        )
        dataset, step_id = name.split("/", 1)
        rows.append({
            "dataset": dataset,
            "step_id": step_id,
            "n_atoms": d["natoms"],
            "n_mapped": d["n_mapped"],
            "n_broken": d["n_broken"],
            "n_formed": d["n_formed"],
            "charge": d["charge"],
            "multiplicity": d["uhf"] + 1,
            "broken_bonds": broken_str,
            "formed_bonds": formed_str,
            "mechanism": d.get("mechanism", ""),
            "step": d.get("step", ""),
        })

    rows.sort(key=lambda r: (r["dataset"], r["step_id"]))
    cols = ["dataset", "step_id", "n_atoms", "n_mapped", "n_broken", "n_formed",
            "charge", "multiplicity", "broken_bonds", "formed_bonds",
            "mechanism", "step"]
    with OUT_CSV.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {OUT_CSV}  ({len(rows)} rows)")


if __name__ == "__main__":
    main()
