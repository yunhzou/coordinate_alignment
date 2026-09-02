#!/usr/bin/env python3
"""Merge fragment detections and rank joint full-target assemblies."""
from __future__ import annotations

import argparse
import gzip
import heapq
import json
from collections import Counter, defaultdict
from fractions import Fraction
from pathlib import Path

import networkx as nx
from rdkit import Chem

from rxn_core.retrosynthesis.catalog_index import (
    CandidateIndexConfig,
    build_candidate_index,
)
from rxn_core.retrosynthesis.compressed_coverage import (
    CoverageRecommendationConfig,
    assign_candidate_items,
    recommend_compressed_coverage_patterns,
)
from rxn_core.retrosynthesis.ranking import (
    assembly_rank as _assembly_rank,
    build_ranked_assembly as _assembly,
    validate_atom_ownership as _formed_bonds,
)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--parts", required=True)
    parser.add_argument("--target-smiles", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--per-domain-limit", type=int, default=200)
    parser.add_argument("--assembly-limit", type=int, default=20)
    parser.add_argument("--maximum-precursors", type=int, default=3)
    parser.add_argument("--frontier-limit", type=int, default=200)
    parser.add_argument("--pattern-limit", type=int, default=4)
    parser.add_argument("--recommendations-per-pattern", type=int, default=4)
    parser.add_argument("--expected-id", action="append", default=[])
    parser.add_argument("--require-attachment-bonds",
                        action=argparse.BooleanOptionalAction, default=False)
    parser.add_argument("--exclude-target-identity",
                        action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--allow-repeated-precursors",
                        action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--attachment-trim-variants",
                        action=argparse.BooleanOptionalAction, default=False)
    parser.add_argument("--chirality-ranking",
                        action=argparse.BooleanOptionalAction, default=False)
    args = parser.parse_args()

    target_implicit = Chem.MolFromSmiles(args.target_smiles)
    if target_implicit is None:
        raise ValueError("invalid target SMILES")
    target = Chem.AddHs(target_implicit)
    Chem.AssignStereochemistry(target, cleanIt=True, force=True)
    target_key = Chem.MolToSmiles(target_implicit, isomericSmiles=True)
    atom_count = target.GetNumAtoms()
    target_edges = [
        (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
        for bond in target.GetBonds()
    ]

    def construction_pattern_signature(atom_sets):
        atom_sets = [set(atoms) for atoms in atom_sets]
        labels = [Chem.MolFragmentToSmiles(
            target, atomsToUse=sorted(atoms), canonical=True,
            isomericSmiles=True) for atoms in atom_sets]
        owner = {
            atom: module
            for module, atoms in enumerate(atom_sets) for atom in atoms
        }
        module_graph = nx.Graph()
        for module, label in enumerate(labels):
            module_graph.add_node(module, label=label)
        edge_labels = defaultdict(list)
        for bond in target.GetBonds():
            left = owner[bond.GetBeginAtomIdx()]
            right = owner[bond.GetEndAtomIdx()]
            if left == right:
                continue
            edge_labels[tuple(sorted((left, right)))].append(
                str(bond.GetBondType()))
        for (left, right), bond_labels in edge_labels.items():
            module_graph.add_edge(
                left, right, label="|".join(sorted(bond_labels)))
        return (
            tuple(sorted(labels)),
            nx.weisfeiler_lehman_graph_hash(
                module_graph, node_attr="label", edge_attr="label"),
        )

    part_paths = sorted(Path(args.parts).glob("part_*.jsonl.gz"))
    if not part_paths:
        raise RuntimeError("no completed part files found")

    def records():
        for path in part_paths:
            with gzip.open(path, "rt", encoding="utf-8") as stream:
                for line in stream:
                    yield json.loads(line)

    index = build_candidate_index(
        records(),
        target,
        target_key,
        config=CandidateIndexConfig(
            per_domain_limit=args.per_domain_limit,
            exclude_target_identity=args.exclude_target_identity,
            attachment_trim_variants=args.attachment_trim_variants,
            chirality_ranking=args.chirality_ranking,
            expected_ids=tuple(args.expected_id),
        ),
    )
    groups = index.groups
    expected = index.expected
    direct_target_matches = index.direct_target_matches
    counts = index.counts

    def item_set_rank(items, covered_atom_count=0):
        structures = {item["structure_key"] for item in items}
        symmetry_retained = sum(
            item["symmetry_retained_atom_count"] for item in items)
        retained = sum(item["retained_atom_count"] for item in items)
        total = sum(item["total_atom_count"] for item in items)
        symmetry_retention = (
            Fraction(symmetry_retained, total) if total else Fraction())
        set_retention = Fraction(retained, total) if total else Fraction()
        return (
            sum(item["chirality_violations"] for item in items),
            len(structures),
            -symmetry_retention,
            -set_retention,
            -covered_atom_count,
            sum(not item["complete"] for item in items),
            sum(len(item["boundary_bonds"]) for item in items),
            sum(item["leftover_atom_count"] for item in items),
            tuple(item["precursor_id"] for item in items),
        )

    def rank_pattern(pattern, covered_atom_count):
        return (
            -covered_atom_count,
            item_set_rank(
                [groups[signature][0] for signature in pattern],
                covered_atom_count),
        )

    coverage_search = recommend_compressed_coverage_patterns(
        groups,
        atom_count,
        rank_pattern,
        result_limit=max(
            args.pattern_limit * args.recommendations_per_pattern,
            args.assembly_limit,
        ),
        config=CoverageRecommendationConfig(
            maximum_precursors=args.maximum_precursors,
            frontier_limit=args.frontier_limit,
        ),
    )
    covers = coverage_search.patterns
    beam_truncated = coverage_search.truncated

    assemblies = []
    assemblies_by_pattern = defaultdict(list)
    representative_pattern = {}
    cover_count = 0
    expansion_count = 0
    assembly_mapping_variants = 0

    def retain_best(items, limit):
        unique = []
        seen = set()
        for assembly in sorted(items, key=_assembly_rank):
            combination = tuple(
                item["precursor_id"] for item in assembly["precursors"])
            if combination in seen:
                continue
            seen.add(combination)
            unique.append(assembly)
            if len(unique) >= limit:
                break
        return unique

    def best_item_combinations(pools, limit):
        """Yield only requested top combinations from sorted pools."""
        if not pools or any(not pool for pool in pools):
            return
        start = tuple(0 for _pool in pools)
        queue = [(item_set_rank(
            tuple(pool[0] for pool in pools)), start)]
        seen = {start}
        yielded = 0
        while queue and yielded < limit:
            _rank, indices = heapq.heappop(queue)
            items = tuple(
                pool[index] for pool, index in zip(pools, indices))
            yield items
            yielded += 1
            for dimension, pool in enumerate(pools):
                next_index = indices[dimension] + 1
                if next_index >= len(pool):
                    continue
                neighbour = list(indices)
                neighbour[dimension] = next_index
                neighbour = tuple(neighbour)
                if neighbour in seen:
                    continue
                seen.add(neighbour)
                neighbour_items = tuple(
                    candidate_pool[index]
                    for candidate_pool, index in zip(pools, neighbour))
                heapq.heappush(queue, (
                    item_set_rank(neighbour_items), neighbour))

    for cover in covers:
        cover_count += 1
        representative_items = assign_candidate_items(
            [groups[signature][0] for signature in cover], atom_count)
        if representative_items is None:
            raise RuntimeError("compressed coverage pattern has no witness")
        representative_atom_sets = tuple(
            tuple(item["covered_target_atoms"])
            for item in representative_items)
        pattern_key = construction_pattern_signature(
            representative_atom_sets)
        representative_pattern.setdefault(
            pattern_key, representative_atom_sets)
        pattern_assemblies = []
        pools = [groups[signature] for signature in cover]
        item_combinations = best_item_combinations(
            pools, args.recommendations_per_pattern)
        for items in item_combinations:
            expansion_count += 1
            ids = [item["precursor_id"] for item in items]
            if (not args.allow_repeated_precursors and
                    len(set(ids)) != len(ids)):
                continue
            placed_items = assign_candidate_items(items, atom_count)
            if placed_items is None:
                continue
            formed = _formed_bonds(
                placed_items, target_edges, args.require_attachment_bonds)
            if formed is not None:
                pattern_assemblies.append(_assembly(placed_items, formed))
                assembly_mapping_variants += 1
        assemblies_by_pattern[pattern_key].extend(retain_best(
            pattern_assemblies, args.recommendations_per_pattern))

    ranked_patterns = sorted(
        ((pattern, retain_best(items, args.recommendations_per_pattern))
         for pattern, items in assemblies_by_pattern.items() if items),
        key=lambda pair: _assembly_rank(pair[1][0]),
    )[:args.pattern_limit]
    construction_patterns = []
    assemblies = []
    for pattern_index, (pattern_key, variants) in enumerate(ranked_patterns, 1):
        pattern = representative_pattern[pattern_key]
        start_rank = len(assemblies) + 1
        for variant in variants:
            variant["construction_pattern"] = pattern_index
            assemblies.append(variant)
        construction_patterns.append({
            "pattern": pattern_index,
            "coverage_atom_sets": [
                list(atoms) for atoms in pattern
            ],
            "fragment_sizes": [len(atoms) for atoms in pattern],
            "recommendation_ranks": list(range(
                start_rank, start_rank + len(variants))),
            "best_score": variants[0]["score"],
        })
    assemblies = assemblies[:args.assembly_limit]

    expected_assembly = None
    if args.expected_id:
        expected_pools = [expected.get(item, ()) for item in args.expected_id]
        for items in best_item_combinations(
                expected_pools, args.per_domain_limit):
            placed_items = assign_candidate_items(items, atom_count)
            if placed_items is None:
                continue
            formed = _formed_bonds(
                placed_items, target_edges, args.require_attachment_bonds)
            if formed is not None:
                expected_assembly = _assembly(placed_items, formed)
                break

    summaries = []
    for path in part_paths:
        summary_path = path.with_suffix(path.suffix + ".summary.json")
        summary = json.loads(summary_path.read_text())
        if summary["schema"] != "rxn_core.retro_catalog_summary/v3":
            raise ValueError(f"unsupported summary schema: {summary['schema']!r}")
        if summary["target_smiles"] != args.target_smiles:
            raise ValueError("fragment shard target does not match merge target")
        summaries.append(summary)
    scan_counts = Counter()
    for summary in summaries:
        scan_counts.update(summary["counts"])

    report = {
        "schema": "rxn_core.retro_catalog_assembly/v5",
        "target_smiles": args.target_smiles,
        "target_atom_count": atom_count,
        "part_count": len(part_paths),
        "scan_counts": {
            key: scan_counts[key]
            for key in (
                "rows",
                "parse_errors",
                "searched",
                "capped",
                "target_coverage_filtered",
                "matched_precursors",
                "fragment_candidates",
            )
        },
        "merge_counts": dict(counts),
        "coverage_domain_families": len(groups),
        "per_domain_limit": args.per_domain_limit,
        "assembly_limit": args.assembly_limit,
        "maximum_precursors": args.maximum_precursors,
        "recommended_cover_count": cover_count,
        "candidate_choice_evaluations": expansion_count,
        "require_attachment_bonds": args.require_attachment_bonds,
        "exclude_target_identity": args.exclude_target_identity,
        "allow_repeated_precursors": args.allow_repeated_precursors,
        "attachment_trim_variants": args.attachment_trim_variants,
        "chirality_ranking": args.chirality_ranking,
        "frontier_limit": args.frontier_limit,
        "recommendation_search_truncated": beam_truncated,
        "pattern_limit": args.pattern_limit,
        "recommendations_per_pattern": args.recommendations_per_pattern,
        "construction_patterns": construction_patterns,
        "direct_target_matches": direct_target_matches,
        "assembly_mapping_variants": assembly_mapping_variants,
        "expected_ids": args.expected_id,
        "expected_ids_found": {
            item: bool(expected.get(item)) for item in args.expected_id
        },
        "expected_assembly_recovered": expected_assembly is not None,
        "expected_assembly": expected_assembly,
        "assemblies": assemblies,
    }
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({
        "output": str(output),
        "scan_counts": dict(scan_counts),
        "assemblies": len(assemblies),
        "expected_ids_found": report["expected_ids_found"],
        "expected_assembly_recovered": report[
            "expected_assembly_recovered"],
    }, indent=2))


if __name__ == "__main__":
    main()
