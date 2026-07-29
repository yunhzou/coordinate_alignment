#!/usr/bin/env python3
"""Compile finalized batch mechanisms into two aligned XYZ endpoints each."""
from __future__ import annotations

import argparse
import json
import re
import shutil
from pathlib import Path

import numpy as np

from rxn_core.alignment.interpolation import proper_align_coordinates
from rxn_core.chemistry_computations.xyz import parse_xyz, write_xyz_str


def _safe_name(value: str) -> str:
    name = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value)).strip("_")
    return name or "case"


def _records(batch_root: Path):
    records = []
    for tier in ("small", "medium", "large"):
        manifest = json.loads(
            (batch_root / "manifests" / f"{tier}.json").read_text())
        records.extend(manifest["cases"])
    return sorted(records, key=lambda item: int(item["source_index"]))


def compile_batch(batch_root: Path, output: Path, overwrite: bool = False):
    batch_root = batch_root.resolve()
    output = output.resolve()
    if output.exists():
        if not overwrite:
            raise FileExistsError(f"output already exists: {output}")
        shutil.rmtree(output)
    output.mkdir(parents=True)

    case_count = 0
    mechanism_count = 0
    file_count = 0
    for case in _records(batch_root):
        step = str(case["step_id"])
        elements_r, coords_r = parse_xyz(case["reactant_xyz"])
        elements_p, coords_p = parse_xyz(case["product_xyz"])
        stage = json.loads(
            (batch_root / "cases" / step / "rp_stage.json").read_text())
        case_dir = output / _safe_name(step)

        for mechanism in stage.get("mechanisms", []):
            mechanism_id = int(mechanism["id"])
            mapping = {
                int(r): int(p)
                for r, p in mechanism["mapping_RP"].items()
            }
            expected = set(range(len(elements_r)))
            if set(mapping) != expected or set(mapping.values()) != expected:
                raise ValueError(
                    f"{step} mechanism {mechanism_id}: mapping is not bijective")
            if any(elements_r[r] != elements_p[p] for r, p in mapping.items()):
                raise ValueError(
                    f"{step} mechanism {mechanism_id}: element mismatch")

            product_in_r = np.asarray(
                mechanism["product_xyz_in_R"], dtype=float)
            mapped_source = np.asarray(
                [coords_p[mapping[r]] for r in range(len(elements_r))])
            if product_in_r.shape != coords_r.shape or not np.allclose(
                    product_in_r, mapped_source, rtol=0.0, atol=1e-10):
                raise ValueError(
                    f"{step} mechanism {mechanism_id}: stored aligned product "
                    "does not match mapping_RP")

            # Match the default viewer: R atom order followed by a proper rigid
            # alignment of P onto R.  This does not alter internal geometry.
            product_aligned = proper_align_coordinates(product_in_r, coords_r)
            mechanism_dir = case_dir / f"mechanism_{mechanism_id:03d}"
            mechanism_dir.mkdir(parents=True)
            (mechanism_dir / "R.xyz").write_text(write_xyz_str(
                elements_r, coords_r,
                f"{step} mechanism {mechanism_id} R endpoint"))
            (mechanism_dir / "P_aligned.xyz").write_text(write_xyz_str(
                elements_r, product_aligned,
                f"{step} mechanism {mechanism_id} P in R order and pose"))
            mechanism_count += 1
            file_count += 2
        case_count += 1

    return {
        "output": str(output),
        "cases": case_count,
        "mechanisms": mechanism_count,
        "xyz_files": file_count,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--batch-root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    print(json.dumps(compile_batch(
        args.batch_root, args.output, overwrite=args.overwrite), indent=2))


if __name__ == "__main__":
    main()
