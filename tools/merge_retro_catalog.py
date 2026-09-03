#!/usr/bin/env python3
"""Merge fragment detections and rank joint full-target assemblies."""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
import gzip
import heapq
import itertools
import json
import time
from collections import Counter, defaultdict
from fractions import Fraction
from pathlib import Path

import networkx as nx
from rdkit import Chem

from rxn_core.fragment_matching import (
    FragmentDetectionConfig,
    progressive_fragment_matching,
)
from rxn_core.fragment_matching.rdkit_adapter import molecule_to_weighted_graph
from rxn_core.retrosynthesis.catalog_index import (
    CandidateIndexConfig,
    build_candidate_index,
    merge_candidate_indexes,
)
from rxn_core.retrosynthesis.compressed_coverage import (
    CoverageRecommendationConfig,
    assign_candidate_items,
    place_candidate_items,
    recommend_compressed_coverage_patterns,
)
from rxn_core.retrosynthesis.ranking import (
    assembly_rank as _assembly_rank,
    build_ranked_assembly as _assembly,
    validate_atom_ownership as _formed_bonds,
)


def _part_records(path):
    with gzip.open(path, "rt", encoding="utf-8") as stream:
        for line in stream:
            yield json.loads(line)


def _build_part_index(task):
    path, target_smiles, target_key, config = task
    target = Chem.AddHs(Chem.MolFromSmiles(target_smiles))
    Chem.AssignStereochemistry(target, cleanIt=True, force=True)
    return build_candidate_index(
        _part_records(path), target, target_key, config=config)


def _serializable_assembly(assembly):
    """Drop exhaustive occupation lists after a concrete witness is selected."""
    if assembly is None:
        return None
    clean = dict(assembly)
    clean["precursors"] = []
    for precursor in assembly["precursors"]:
        item = dict(precursor)
        item.pop("target_occupations", None)
        clean["precursors"].append(item)
    return clean


def _progressive_item(base, placement, molecule):
    graph = molecule_to_weighted_graph(molecule).to_networkx()
    mapping = dict(placement.mapping)
    retained = set(mapping)
    leftover = set(graph) - retained
    retained_fragments = tuple(sorted(
        (tuple(sorted(fragment)) for fragment in placement.retained_fragments),
        key=lambda fragment: (fragment[0], len(fragment)),
    ))
    retained_fragment_by_atom = {
        atom: index
        for index, fragment in enumerate(retained_fragments)
        for atom in fragment
    }
    leftover_fragments = tuple(sorted(
        (tuple(sorted(component)) for component in
         nx.connected_components(graph.subgraph(leftover))),
        key=lambda fragment: (fragment[0], len(fragment)),
    )) if leftover else ()
    boundary = tuple(sorted(
        tuple(sorted((left, right))) for left, right in graph.edges()
        if ((left in retained) != (right in retained)
            or (left in retained and right in retained
                and retained_fragment_by_atom[left]
                != retained_fragment_by_atom[right]))
    ))
    attachment_source = tuple(sorted({
        atom for edge in boundary for atom in edge if atom in retained
    }))
    retained_heavy = sum(
        molecule.GetAtomWithIdx(atom).GetAtomicNum() > 1 for atom in retained)
    transformed = dict(base)
    transformed.update({
        "mapping": [list(pair) for pair in placement.mapping],
        "covered_target_atoms": sorted(mapping.values()),
        "retained_atoms": sorted(retained),
        "retained_fragments": [list(part) for part in retained_fragments],
        "leftover_fragments": [list(part) for part in leftover_fragments],
        "leftover_atom_count": len(leftover),
        "boundary_bonds": [list(edge) for edge in boundary],
        "attachment_atoms_source": list(attachment_source),
        "attachment_atoms_target": sorted(
            mapping[atom] for atom in attachment_source),
        "target_domains": [
            [source, [target]] for source, target in placement.mapping
        ],
        "target_occupations": ({
            "covered_target_atoms": tuple(sorted(mapping.values())),
            "mapping": placement.mapping,
            "attachment_atoms_target": tuple(sorted(
                mapping[atom] for atom in attachment_source)),
        },),
        "retained_heavy_atoms": retained_heavy,
        "retained_atom_count": len(retained),
        "heavy_atom_retention": retained_heavy / molecule.GetNumHeavyAtoms(),
        "atom_retention": len(retained) / molecule.GetNumAtoms(),
        "symmetry_copy_capacity": 1,
        "symmetry_retained_atoms": sorted(retained),
        "symmetry_retained_atom_count": len(retained),
        "symmetry_retained_heavy_atoms": retained_heavy,
        "symmetry_atom_retention": len(retained) / molecule.GetNumAtoms(),
        "symmetry_heavy_atom_retention": (
            retained_heavy / molecule.GetNumHeavyAtoms()),
        "attachment_trimmed_target_atoms": [],
    })
    return transformed


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--parts", required=True)
    parser.add_argument("--target-smiles", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--per-domain-limit", type=int, default=200)
    parser.add_argument("--assembly-limit", type=int, default=20)
    parser.add_argument("--maximum-precursors", type=int)
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
    parser.add_argument("--index-workers", type=int, default=1)
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
        memberships = defaultdict(list)
        for module, atoms in enumerate(atom_sets):
            for atom in atoms:
                memberships[atom].append(module)
        module_graph = nx.Graph()
        for module, label in enumerate(labels):
            module_graph.add_node(module, label=label)
        edge_labels = defaultdict(list)
        for atom, modules in memberships.items():
            for left, right in itertools.combinations(modules, 2):
                symbol = target.GetAtomWithIdx(atom).GetSymbol()
                edge_labels[(left, right)].append(f"overlap:{symbol}")
        for bond in target.GetBonds():
            for left in memberships[bond.GetBeginAtomIdx()]:
                for right in memberships[bond.GetEndAtomIdx()]:
                    if left != right:
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

    index_config = CandidateIndexConfig(
        per_domain_limit=args.per_domain_limit,
        exclude_target_identity=args.exclude_target_identity,
        attachment_trim_variants=args.attachment_trim_variants,
        chirality_ranking=args.chirality_ranking,
        expected_ids=tuple(args.expected_id),
    )
    index_start = time.perf_counter()
    tasks = tuple(
        (path, args.target_smiles, target_key, index_config)
        for path in part_paths)
    with ProcessPoolExecutor(max_workers=args.index_workers) as executor:
        shard_indexes = tuple(executor.map(_build_part_index, tasks))
    index = merge_candidate_indexes(
        shard_indexes, args.per_domain_limit)
    groups = index.groups
    expected = index.expected
    direct_target_matches = index.direct_target_matches
    counts = index.counts
    print(json.dumps({
        "phase": "candidate_index",
        "seconds": time.perf_counter() - index_start,
        "coverage_families": len(groups),
    }), flush=True)

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

    coverage_start = time.perf_counter()
    coverage_search = recommend_compressed_coverage_patterns(
        groups,
        atom_count,
        rank_pattern,
        result_limit=max(args.assembly_limit, atom_count),
        config=CoverageRecommendationConfig(
            maximum_precursors=args.maximum_precursors,
        ),
    )
    covers = coverage_search.patterns
    search_truncated = coverage_search.truncated
    print(json.dumps({
        "phase": "coverage_recommendation",
        "seconds": time.perf_counter() - coverage_start,
        "patterns": len(covers),
        "truncated": search_truncated,
    }), flush=True)

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

    assembly_start = time.perf_counter()
    for cover, representative_atom_sets in zip(
            covers, coverage_search.occupations, strict=True):
        cover_count += 1
        pattern_key = construction_pattern_signature(
            representative_atom_sets)
        representative_pattern.setdefault(
            pattern_key, representative_atom_sets)
        pattern_assemblies = []
        pools = [groups[signature] for signature in cover]
        unique_signatures = tuple(dict.fromkeys(cover))
        signature_positions = {
            signature: tuple(
                index for index, item in enumerate(cover)
                if item == signature)
            for signature in unique_signatures
        }

        def candidate_combinations():
            seen = set()
            unique_pools = [groups[signature]
                            for signature in unique_signatures]
            homogeneous_limit = (
                args.per_domain_limit
                if len(unique_signatures) <= 2
                else max(
                    args.recommendations_per_pattern,
                    args.per_domain_limit // len(unique_signatures),
                )
            )
            for unique_items in best_item_combinations(
                    unique_pools, homogeneous_limit):
                items = [None] * len(cover)
                for signature, item in zip(
                        unique_signatures, unique_items):
                    for position in signature_positions[signature]:
                        items[position] = item
                key = tuple(item["precursor_id"] for item in items)
                if key not in seen:
                    seen.add(key)
                    yield tuple(items)
            for items in best_item_combinations(
                    pools, args.recommendations_per_pattern):
                key = tuple(item["precursor_id"] for item in items)
                if key not in seen:
                    seen.add(key)
                    yield items

        for items in candidate_combinations():
            expansion_count += 1
            ids = [item["precursor_id"] for item in items]
            if (not args.allow_repeated_precursors and
                    len(set(ids)) != len(ids)):
                continue
            placed_items = place_candidate_items(
                items, representative_atom_sets)
            formed = _formed_bonds(
                placed_items, target_edges, args.require_attachment_bonds)
            if formed is not None:
                pattern_assemblies.append(_assembly(placed_items, formed))
                assembly_mapping_variants += 1
        assemblies_by_pattern[pattern_key].extend(retain_best(
            pattern_assemblies, args.per_domain_limit))
        if cover_count % 10 == 0:
            print(json.dumps({
                "phase": "candidate_expansion",
                "seconds": time.perf_counter() - assembly_start,
                "patterns_processed": cover_count,
            }), flush=True)

    print(json.dumps({
        "phase": "candidate_expansion",
        "seconds": time.perf_counter() - assembly_start,
        "patterns_processed": cover_count,
        "complete": True,
    }), flush=True)

    ranked_patterns = sorted(
        ((pattern, retain_best(items, args.per_domain_limit))
         for pattern, items in assemblies_by_pattern.items() if items),
        key=lambda pair: _assembly_rank(pair[1][0]),
    )
    construction_patterns = []
    assemblies = []
    all_pattern_assemblies = []
    for pattern_index, (pattern_key, variants) in enumerate(ranked_patterns, 1):
        pattern = representative_pattern[pattern_key]
        start_rank = len(assemblies) + 1
        for variant in variants:
            variant["construction_pattern"] = pattern_index
        all_pattern_assemblies.extend(variants)
        displayed_variants = variants[:args.recommendations_per_pattern]
        assemblies.extend(displayed_variants)
        construction_patterns.append({
            "pattern": pattern_index,
            "coverage_atom_sets": [
                list(atoms) for atoms in pattern
            ],
            "fragment_sizes": [len(atoms) for atoms in pattern],
            "recommendation_ranks": list(range(
                start_rank, start_rank + len(displayed_variants))),
            "best_score": variants[0]["score"],
        })
    assemblies = assemblies[:args.assembly_limit]

    expected_stoichiometry = Counter(args.expected_id)
    naturally_ranked = retain_best(
        all_pattern_assemblies, len(all_pattern_assemblies))
    expected_recommendation_rank = next((
        rank for rank, assembly in enumerate(naturally_ranked, 1)
        if Counter(
            item["precursor_id"] for item in assembly["precursors"]
        ) == expected_stoichiometry
    ), None) if expected_stoichiometry else None
    if expected_recommendation_rank is not None:
        expected_recommendation = naturally_ranked[
            expected_recommendation_rank - 1]
        if all(
                Counter(item["precursor_id"] for item in assembly["precursors"])
                != expected_stoichiometry
                for assembly in assemblies):
            assemblies.append(expected_recommendation)

    expected_assembly = None
    expected_mapping = None
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
        if all(expected_pools):
            base_items = [pool[0] for pool in expected_pools]
            molecules = [
                Chem.AddHs(Chem.MolFromSmiles(item["smiles"]))
                for item in base_items
            ]
            precise = progressive_fragment_matching(
                tuple(
                    (item["precursor_id"], molecule_to_weighted_graph(molecule))
                    for item, molecule in zip(
                        base_items, molecules, strict=True)
                ),
                molecule_to_weighted_graph(target),
                config=FragmentDetectionConfig(
                    seed_mode="all",
                    branch_limit=100,
                    candidate_limit=512,
                    iso_tolerance=0.5,
                ),
            )
            precise_items = tuple(
                _progressive_item(item, placement, molecule)
                for item, placement, molecule in zip(
                    base_items, precise.placements, molecules, strict=True)
            )
            formed = _formed_bonds(
                precise_items, target_edges, require_attachment_bonds=False)
            expected_mapping = _assembly(
                precise_items, formed if formed is not None else [])
            covered = {
                atom for item in precise_items
                for atom in item["covered_target_atoms"]
            }
            expected_mapping["score"]["covered_target_atoms"] = len(covered)
            expected_mapping["score"]["target_atom_count"] = atom_count
            expected_mapping["uncovered_target_atoms"] = list(
                precise.uncovered_target_atoms)

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
        "recommendation_search_truncated": search_truncated,
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
        "expected_assembly": _serializable_assembly(expected_assembly),
        "expected_mapping": _serializable_assembly(expected_mapping),
        "expected_recommendation_rank": expected_recommendation_rank,
        "assemblies": [
            _serializable_assembly(assembly) for assembly in assemblies
        ],
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
