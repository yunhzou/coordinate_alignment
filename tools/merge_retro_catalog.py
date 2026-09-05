#!/usr/bin/env python3
"""Exact assembly of saved detections; output limits never prune the search."""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
import gzip
import json
from pathlib import Path
import time

from rdkit import Chem

from rxn_core.retrosynthesis.assembly import AssemblyProblem
from rxn_core.retrosynthesis.catalog_index import (
    CandidateIndexConfig, build_candidate_index, merge_candidate_indexes,
)
from rxn_core.retrosynthesis.ranking import assembly_rank


def _part_records(path):
    with gzip.open(path, "rt", encoding="utf-8") as stream:
        for line in stream:
            yield json.loads(line)


def _build_part_index(task):
    path, smiles, key, config = task
    target = Chem.AddHs(Chem.MolFromSmiles(smiles))
    return build_candidate_index(_part_records(path), target, key, config=config)


def _serializable_assembly(assembly):
    clean = dict(assembly)
    clean["precursors"] = [
        {key: value for key, value in item.items() if key != "target_occupations"}
        for item in assembly["precursors"]]
    return clean


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--parts", required=True)
    parser.add_argument("--target-smiles", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--assembly-limit", type=int, default=20, help="display only")
    parser.add_argument("--pattern-limit", type=int, default=4, help="display only")
    parser.add_argument("--recommendations-per-pattern", type=int, default=4, help="display only")
    parser.add_argument("--expected-id", action="append", default=[])
    parser.add_argument("--index-workers", type=int, default=1)
    parser.add_argument("--exhaustive", action="store_true", help="write every assembly, not just the certified recommendation prefix")
    args = parser.parse_args()
    implicit = Chem.MolFromSmiles(args.target_smiles)
    if implicit is None:
        raise ValueError("invalid target SMILES")
    target = Chem.AddHs(implicit)
    target_key = Chem.MolToSmiles(implicit, isomericSmiles=True)
    paths = sorted(Path(args.parts).glob("part_*.jsonl.gz"))
    if not paths:
        raise ValueError("no detection shards found")
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    checkpoint = output.with_suffix(".checkpoint.json")
    started = time.perf_counter()
    summaries = [json.loads(path.with_suffix(path.suffix + ".summary.json").read_text())
                 for path in paths]
    if any(s["target_smiles"] != args.target_smiles for s in summaries):
        raise ValueError("detection shard target differs from assembly target")
    scan_counts = Counter()
    for summary in summaries:
        scan_counts.update(summary["counts"])
    config = CandidateIndexConfig(expected_ids=tuple(args.expected_id),
                                  exclude_target_identity=False)
    tasks = [(path, args.target_smiles, target_key, config) for path in paths]
    with ProcessPoolExecutor(max_workers=args.index_workers) as executor:
        index = merge_candidate_indexes(executor.map(_build_part_index, tasks))
    index_seconds = time.perf_counter() - started
    # Save the lossless occupation index BEFORE constructing the cover graph.
    index_path = output.with_suffix(".occupations.json")
    index_path.write_text(json.dumps({
        "schema": "rxn_core.occupation_index/v1",
        "parts": [str(p.resolve()) for p in paths],
        "target_smiles": args.target_smiles,
        "groups": list(index.groups.values()),
    }) + "\n")
    checkpoint.write_text(json.dumps({"phase": "index_complete", "complete": False,
                                      "index_seconds": index_seconds}) + "\n")
    problem = AssemblyProblem.from_index(index, target)
    graph_path = output.with_suffix(".decisions.json")
    graph_path.write_text(json.dumps(problem.decisions.to_record()) + "\n")
    best_by_pattern = defaultdict(list)
    pattern_counts = Counter()
    supplier_sets = defaultdict(set)
    expected = Counter(args.expected_id)
    expected_best = None
    count = 0
    evaluations_path = output.with_suffix(".assemblies.jsonl")
    assembly_started = time.perf_counter()
    selected_patterns = []
    exhausted = True
    with evaluations_path.open("w") as stream:
        for assembly in problem.ranked_assemblies():
            count += 1
            serial = _serializable_assembly(assembly)
            stream.write(json.dumps(serial) + "\n")
            key = assembly["pattern_key"]
            if key not in selected_patterns and len(selected_patterns) < args.pattern_limit:
                selected_patterns.append(key)
            supplier_key = tuple(sorted(assembly["precursor_stoichiometry"].items()))
            new_suppliers = supplier_key not in supplier_sets[key]
            supplier_sets[key].add(supplier_key)
            pattern_counts[key] += new_suppliers
            bucket = best_by_pattern[key]
            if new_suppliers:
                bucket.append(serial)
            bucket.sort(key=assembly_rank)
            del bucket[args.recommendations_per_pattern:]
            if expected and Counter(assembly["precursor_stoichiometry"]) == expected:
                if expected_best is None or assembly_rank(serial) < assembly_rank(expected_best):
                    expected_best = serial
            # Time-based progress, not a search budget.
            now = time.perf_counter()
            if now - getattr(main, "_last_progress", assembly_started) >= 30:
                stream.flush()
                progress = {"phase": "assembly", "complete": False, "evaluated": count,
                            "seconds": now - assembly_started}
                checkpoint.write_text(json.dumps(progress) + "\n")
                print(json.dumps(progress), flush=True)
                main._last_progress = now
            if (not args.exhaustive and len(selected_patterns) == args.pattern_limit
                    and all(pattern_counts[p] >= args.recommendations_per_pattern for p in selected_patterns)):
                exhausted = False
                break
    pattern_order = sorted(best_by_pattern, key=lambda p: (assembly_rank(best_by_pattern[p][0]), p))
    assemblies, patterns = [], []
    for number, key in enumerate(pattern_order[:args.pattern_limit], 1):
        ranks = []
        variants = best_by_pattern[key]
        for assembly in variants:
            if len(assemblies) == args.assembly_limit:
                break
            assembly["construction_pattern"] = number
            assemblies.append(assembly)
            ranks.append(len(assemblies))
        if ranks:
            patterns.append({
                "pattern": number,
                "coverage_atom_sets": [item["covered_target_atoms"] for item in variants[0]["precursors"]],
                "matched_fragment_pattern": [item["target_fragment_atoms"] for item in variants[0]["precursors"]],
                "pattern_certificate": key,
                "fragment_sizes": [len(item["covered_target_atoms"]) for item in variants[0]["precursors"]],
                "recommendation_ranks": ranks,
                "assembly_count": pattern_counts[key],
                "best_score": variants[0]["score"],
            })
    expected_rank = None
    if expected_best is not None:
        # Rank in the SAME saved stream, not a separate rematching experiment.
        expected_rank = 1
        with evaluations_path.open() as stream:
            for line in stream:
                if assembly_rank(json.loads(line)) < assembly_rank(expected_best):
                    expected_rank += 1
    report = {
        "schema": "rxn_core.retro_catalog_assembly/v5",
        "target_smiles": args.target_smiles, "target_atom_count": target.GetNumAtoms(),
        "part_count": len(paths), "scan_counts": dict(scan_counts), "merge_counts": index.counts,
        "coverage_domain_families": len(index.groups),
        "assembly_limit": args.assembly_limit, "pattern_limit": args.pattern_limit,
        "recommendations_per_pattern": args.recommendations_per_pattern,
        "allow_repeated_precursors": True, "require_attachment_bonds": False,
        "exclude_target_identity": False, "attachment_trim_variants": False, "chirality_ranking": False,
        "recommendation_search_truncated": False,
        "detection_complete": not (index.counts.get("incomplete_precursors", 0)
            or any(scan_counts[key] for key in
                   ("incomplete", "capped", "parse_errors", "target_coverage_filtered"))),
        "assembly_complete": exhausted,
        "recommendations_certified": True,
        "candidate_choice_evaluations": count, "assembly_mapping_variants": count,
        "construction_patterns": patterns, "assemblies": assemblies,
        "expected_ids": args.expected_id,
        "expected_ids_found": {item: bool(index.expected.get(item)) for item in args.expected_id},
        "expected_recommendation_rank": expected_rank,
        "expected_assembly_recovered": expected_best is not None,
        "expected_status": ("in_ranked_prefix" if expected_best else
                            "no_cover_in_saved_detections" if exhausted else "not_in_ranked_prefix"),
        "expected_assembly": expected_best, "expected_mapping": None,
        "direct_target_matches": index.direct_target_matches,
        "saved_occupations": str(index_path.resolve()),
        "saved_decisions": str(graph_path.resolve()),
        "saved_assemblies": str(evaluations_path.resolve()),
        "timing": {"index_seconds": index_seconds,
                   "assembly_seconds": time.perf_counter() - assembly_started},
        "semantics": "Geometric fragment support, not a chemically validated or atom-owned reaction.",
    }
    output.write_text(json.dumps(report, indent=2) + "\n")
    checkpoint.write_text(json.dumps({"phase": "complete", "complete": True,
                                      "evaluated": count}) + "\n")
    print(json.dumps({"output": str(output), "evaluated": count, "timing": report["timing"]}))


if __name__ == "__main__":
    main()
