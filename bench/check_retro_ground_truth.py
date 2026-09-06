#!/usr/bin/env python3
"""Run the production detector/assembler on bank-verified known ingredients.

No atom mappings or target regions are supplied. Detection checkpoints survive
assembly failures. This checks recoverability, not blind full-bank rank.
"""
import argparse
import csv
import gzip
import json
import pickle
from functools import reduce
from itertools import product
from pathlib import Path
import subprocess
import sys

from rdkit import Chem


def audit_saved_results(directory):
    """Check the union after the production index expands allowed occupations.

    A missing atom here proves no combination, even with unlimited repetitions,
    can cover P using these saved detections. A full union alone proves no
    particular stoichiometric combination.
    """
    report = json.loads((directory / "results.json").read_text())
    target = Chem.AddHs(Chem.MolFromSmiles(report["target_smiles"]))
    index = json.loads((directory / "results.occupations.json").read_text())
    support = {}
    for group in index["groups"]:
        for item in group:
            atoms = support.setdefault(item["precursor_id"], set())
            for occupation in item["target_occupations"]:
                atoms.update(occupation["covered_target_atoms"])
    union = set().union(*support.values())
    missing = sorted(set(range(target.GetNumAtoms())) - union)
    sources = []
    for path in sorted((directory / "checkpoints").glob("*.detection.pkl.gz")):
        row, (_, result, seconds) = pickle.load(gzip.open(path, "rb"))
        initial = [graph for graph in result.search_graphs
                   if graph.contexts[0].objective == "seeded_fragment"]
        sources.append({"id": row[2], "detection_seconds": seconds,
            "candidate_count": len(result.candidates), "status": result.status,
            "complete": result.complete, "maximum_branch_count": result.maximum_branch_count,
            "capped_seed_count": result.capped_seed_count,
            "seed_attempt_count": result.seed_attempt_count,
            "supported_target_atoms": sorted(support[row[2]]),
            "initial_searches": [{"seed": graph.contexts[0].seed_order,
                "stops": [stop.reason for stop in graph.stops],
                "placements": [dict(edge.match["symmetry"]["witness"])
                               for edge in graph.transitions if edge.match]}
                for graph in initial]})
    audit = {"scope": report["search_scope"], "target_atom_count": target.GetNumAtoms(),
             "union_covered_atom_count": len(union), "uncovered_target_atoms": missing,
             "uncovered_elements": [target.GetAtomWithIdx(a).GetSymbol() for a in missing],
             "expected_status": report["expected_status"], "sources": sources}
    (directory / "correctness_audit.json").write_text(json.dumps(audit, indent=2) + "\n")
    report["uncovered_target_atoms"] = missing
    (directory / "results.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({key: value for key, value in audit.items() if key != "sources"}), flush=True)


def build_partial_diagnostic(directory):
    """Maximum union for the specified known copies; never a recommendation.

    Equal coverage masks may share one display witness here because the sole
    diagnostic objective is coverage. Production matching/indexing is unchanged.
    """
    from rxn_core.retrosynthesis.compressed_coverage import place_item
    from rxn_core.retrosynthesis.ranking import build_ranked_assembly
    report = json.loads((directory / "results.json").read_text())
    index = json.loads((directory / "results.occupations.json").read_text())
    pools = {source_id: {} for source_id in report["expected_ids"]}
    for group in index["groups"]:
        for item in group:
            if item["precursor_id"] not in pools:
                continue
            for occupation in item["target_occupations"]:
                mask = sum(1 << a for a in occupation["covered_target_atoms"])
                pools[item["precursor_id"]].setdefault(mask, (item, occupation))
    choices = [tuple(pools[source_id]) for source_id in report["expected_ids"]]
    best, selected = None, None
    for masks in product(*choices):
        union = reduce(int.__or__, masks, 0)
        rank = (-union.bit_count(), sum(m.bit_count() for m in masks) - union.bit_count())
        if best is None or rank < best:
            best, selected = rank, masks
    placed = [place_item(*pools[source_id][mask])
              for source_id, mask in zip(report["expected_ids"], selected, strict=True)]
    target = Chem.AddHs(Chem.MolFromSmiles(report["target_smiles"]))
    covered = set().union(*(set(item["covered_target_atoms"]) for item in placed))
    carried = {tuple(sorted(b)) for item in placed for b in item["preserved_target_bonds"]}
    connections = [sorted((b.GetBeginAtomIdx(), b.GetEndAtomIdx())) for b in target.GetBonds()
        if b.GetBeginAtomIdx() in covered and b.GetEndAtomIdx() in covered
        and tuple(sorted((b.GetBeginAtomIdx(), b.GetEndAtomIdx()))) not in carried]
    diagnostic = build_ranked_assembly(placed, connections)
    diagnostic["score"].update(covered_target_atoms=len(covered), target_atom_count=target.GetNumAtoms())
    diagnostic["uncovered_target_atoms"] = sorted(set(range(target.GetNumAtoms())) - covered)
    diagnostic["precursors"] = [{k: v for k, v in item.items() if k != "target_occupations"}
                                for item in diagnostic["precursors"]]
    diagnostic["semantics"] = "Best partial support from one copy of each known ingredient; not a complete recommendation"
    report["diagnostic_assembly"] = diagnostic
    (directory / "results.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"diagnostic_coverage": len(covered),
        "target_atoms": target.GetNumAtoms(), "uncovered": diagnostic["uncovered_target_atoms"]}), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("phase", choices=("prepare", "detect", "assemble", "audit", "diagnostic"))
    parser.add_argument("--case", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--bank", type=Path)
    parser.add_argument("--source-index", type=int)
    parser.add_argument("--workers", type=int, default=1)
    args = parser.parse_args()
    case = json.loads(args.case.read_text())
    directory = args.output_dir.resolve()
    tools = Path(__file__).resolve().parents[1] / "tools"
    catalog = directory / "verified_bank_subset.csv.gz"
    expected = {r["id"]: r for r in case["reactants"]}
    if args.phase == "prepare":
        with gzip.open(args.bank, "rt") as stream:
            rows = [(i, row) for i, row in enumerate(csv.DictReader(stream))
                    if row["Bank ID"] in expected]
        assert len(rows) == len(expected)
        for _, row in rows:
            normalize = lambda s: Chem.MolToSmiles(Chem.MolFromSmiles(s))
            assert normalize(row["SMILES"]) == normalize(expected[row["Bank ID"]]["smiles"])
        directory.mkdir(parents=True, exist_ok=True)
        (directory / "parts").mkdir()
        (directory / "checkpoints").mkdir()
        (directory / "logs").mkdir()
        with gzip.open(catalog, "xt") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0][1]))
            writer.writeheader()
            writer.writerows(row for _, row in rows)
        provenance = {"case": case, "bank": str(args.bank.resolve()),
                      "bank_rows": [{"bank_row_index": i, **row} for i, row in rows]}
        (directory / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
        print(json.dumps(provenance), flush=True)
    elif args.phase == "detect":
        subprocess.run([sys.executable, str(tools / "search_mcule_retro.py"),
            "--target-smiles", case["target_smiles"], "--catalog", str(catalog),
            "--catalog-format", "csv", "--id-column", "Bank ID",
            "--output", str(directory / "parts" / f"part_{args.source_index}.jsonl.gz"),
            "--checkpoint-directory", str(directory / "checkpoints"),
            "--shard-count", str(len(expected)), "--shard-index", str(args.source_index),
            "--workers", str(args.workers), "--scheduling", "adaptive", "--batch-size", "1",
            "--seed-mode", "all", "--branch-limit", "100", "--iso-tolerance", "0.5"], check=True)
    elif args.phase == "audit":
        audit_saved_results(directory)
    elif args.phase == "diagnostic":
        build_partial_diagnostic(directory)
    else:
        output = directory / "results.json"
        expected_args = [word for source_id in expected for word in ("--expected-id", source_id)]
        subprocess.run([sys.executable, str(tools / "merge_retro_catalog.py"),
            "--parts", str(directory / "parts"), "--target-smiles", case["target_smiles"],
            "--output", str(output), "--index-workers", str(args.workers),
            "--pattern-limit", "1", "--recommendations-per-pattern", "1", *expected_args], check=True)
        report = json.loads(output.read_text())
        report["search_scope"] = case["scope"]
        output.write_text(json.dumps(report, indent=2) + "\n")
        audit_saved_results(directory)
        subprocess.run([sys.executable, str(tools / "build_retro_db_viewer.py"),
            "--results", str(output), "--output", str(directory / "viewer.html"),
            "--title", "Example 5: known-ingredient correctness check (3-entry bank)",
            "--ground-truth-status", report["expected_status"],
            "--ground-truth-note", case["scope"]], check=True)


if __name__ == "__main__":
    main()
