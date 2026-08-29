#!/usr/bin/env python3
"""Merge catalog shards and rank complementary two-precursor covers."""
from __future__ import annotations

import argparse
import gzip
import itertools
import json
from collections import Counter, defaultdict
from fractions import Fraction
from functools import lru_cache
from pathlib import Path

import networkx as nx
from rdkit import Chem


def _mask(atoms):
    value = 0
    for atom in atoms:
        value |= 1 << int(atom)
    return value


def _cip(atom):
    return atom.GetProp("_CIPCode") if atom.HasProp("_CIPCode") else None


def _chirality_violations(candidate, source, target):
    reaction_center = {int(edge[0]) for edge in candidate["boundary_bonds"]}
    violations = 0
    for source_atom, target_atom in candidate["mapping"]:
        if int(source_atom) in reaction_center:
            continue
        source_cip = _cip(source.GetAtomWithIdx(int(source_atom)))
        target_cip = _cip(target.GetAtomWithIdx(int(target_atom)))
        if ((source_cip is not None or target_cip is not None)
                and source_cip != target_cip):
            violations += 1
    return violations


def _entry(record, candidate, explicit_molecule, total_heavy_atoms,
           structure_key, target, chirality_ranking):
    leftovers = candidate["leftover_fragments"]
    retained_heavy_atoms = sum(
        explicit_molecule.GetAtomWithIdx(int(atom)).GetAtomicNum() > 1
        for atom in candidate["retained_atoms"]
    )
    retained_atom_count = len(candidate["retained_atoms"])
    total_atom_count = explicit_molecule.GetNumAtoms()
    return {
        "precursor_id": record["precursor_id"],
        "smiles": record["smiles"],
        "row_index": record["row_index"],
        "complete": record["complete"],
        "status": record["status"],
        "best_fragment_size": record["best_fragment_size"],
        "covered_target_atoms": candidate["covered_target_atoms"],
        "retained_atoms": candidate["retained_atoms"],
        "leftover_fragments": leftovers,
        "boundary_bonds": candidate["boundary_bonds"],
        "attachment_atoms_P": candidate["attachment_atoms_P"],
        "mapping": candidate["mapping"],
        "retained_fragments": candidate.get(
            "retained_fragments", [candidate["retained_atoms"]]),
        "leftover_atom_count": sum(map(len, leftovers)),
        "structure_key": structure_key,
        "retained_heavy_atoms": retained_heavy_atoms,
        "total_heavy_atoms": total_heavy_atoms,
        "heavy_atom_retention": retained_heavy_atoms / total_heavy_atoms,
        "retained_atom_count": retained_atom_count,
        "total_atom_count": total_atom_count,
        "atom_retention": retained_atom_count / total_atom_count,
        "attachment_trimmed_target_atoms": candidate.get(
            "attachment_trimmed_target_atoms", []),
        "chirality_violations": (
            _chirality_violations(candidate, explicit_molecule, target)
            if chirality_ranking else 0),
    }


def _attachment_trim_variants(candidate, molecule):
    """Yield connected variants that relinquish one mapped attachment atom.

    This resolves atom-source competition at a coupling junction without
    permitting arbitrary overlap between precursor fragments.
    """
    yield candidate
    mapping = {int(source): int(target)
               for source, target in candidate["mapping"]}
    retained = {int(atom) for atom in candidate["retained_atoms"]}
    attachment_targets = set(map(int, candidate["attachment_atoms_P"]))
    inverse = {target: source for source, target in mapping.items()}
    seen = set()
    for target_atom in sorted(attachment_targets):
        source_atom = inverse.get(target_atom)
        if source_atom is None or source_atom not in retained:
            continue
        new_retained = retained - {source_atom}
        if not new_retained:
            continue
        graph = nx.Graph()
        graph.add_nodes_from(new_retained)
        for bond in molecule.GetBonds():
            left, right = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            if left in new_retained and right in new_retained:
                graph.add_edge(left, right)
        if not nx.is_connected(graph):
            continue
        new_mapping = {source: target for source, target in mapping.items()
                       if source in new_retained}
        leftover = set(range(molecule.GetNumAtoms())) - new_retained
        leftover_graph = nx.Graph()
        leftover_graph.add_nodes_from(leftover)
        boundary = []
        attachment_r = set()
        for bond in molecule.GetBonds():
            left, right = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            if left in leftover and right in leftover:
                leftover_graph.add_edge(left, right)
            elif left in new_retained and right in leftover:
                boundary.append([left, right])
                attachment_r.add(left)
            elif right in new_retained and left in leftover:
                boundary.append([right, left])
                attachment_r.add(right)
        if not attachment_r:
            continue
        key = tuple(sorted(new_mapping.items()))
        if key in seen:
            continue
        seen.add(key)
        variant = dict(candidate)
        variant.update({
            "mapping": [list(item) for item in sorted(new_mapping.items())],
            "retained_atoms": sorted(new_retained),
            "covered_target_atoms": sorted(new_mapping.values()),
            "leftover_fragments": [
                sorted(component)
                for component in nx.connected_components(leftover_graph)
            ],
            "boundary_bonds": sorted(boundary),
            "attachment_atoms_R": sorted(attachment_r),
            "attachment_atoms_P": sorted(
                new_mapping[atom] for atom in attachment_r),
            "attachment_trimmed_target_atoms": [target_atom],
        })
        yield variant


def _entry_rank(entry):
    return (
        entry["chirality_violations"],
        not entry["complete"],
        -entry["atom_retention"],
        len(entry["boundary_bonds"]),
        entry["leftover_atom_count"],
        len(entry["leftover_fragments"]),
        entry["precursor_id"],
        entry["mapping"],
    )


def _formed_bonds(precursors, target_edges, require_attachment_bonds=True):
    owner = {}
    attachments = []
    for index, precursor in enumerate(precursors):
        attachments.append(set(precursor["attachment_atoms_P"]))
        for atom in precursor["covered_target_atoms"]:
            if atom in owner:
                return None
            owner[atom] = index
    formed = []
    for atom_a, atom_b in target_edges:
        owner_a, owner_b = owner.get(atom_a), owner.get(atom_b)
        if owner_a is None or owner_b is None:
            return None
        if owner_a == owner_b:
            continue
        if (require_attachment_bonds and
                (atom_a not in attachments[owner_a] or
                 atom_b not in attachments[owner_b])):
            return None
        formed.append([atom_a, atom_b])
    return formed


def _assembly(items, formed):
    precursors = sorted(items, key=lambda item: item["precursor_id"])
    stoichiometry = Counter(item["precursor_id"] for item in precursors)
    retention_by_structure = defaultdict(list)
    for item in precursors:
        retention_by_structure[item["structure_key"]].append(Fraction(
            item["retained_heavy_atoms"], item["total_heavy_atoms"]))
    unique_retentions = [min(values)
                         for values in retention_by_structure.values()]
    worst_retention = min(unique_retentions)
    mean_retention = sum(unique_retentions, Fraction()) / len(unique_retentions)
    set_retention = Fraction(
        sum(item["retained_heavy_atoms"] for item in precursors),
        sum(item["total_heavy_atoms"] for item in precursors),
    )
    set_atom_retention = Fraction(
        sum(item["retained_atom_count"] for item in precursors),
        sum(item["total_atom_count"] for item in precursors),
    )
    return {
        "precursors": precursors,
        "precursor_stoichiometry": dict(sorted(stoichiometry.items())),
        "formed_bonds": formed,
        "score": {
            "chirality_violations": sum(
                item.get("chirality_violations", 0) for item in precursors),
            "unique_precursor_structures": len(retention_by_structure),
            "set_atom_retention": float(set_atom_retention),
            "set_heavy_atom_retention": float(set_retention),
            "worst_heavy_atom_retention": float(worst_retention),
            "mean_heavy_atom_retention": float(mean_retention),
            "capped_precursors": sum(not item["complete"] for item in precursors),
            "broken_bonds": sum(len(item["boundary_bonds"])
                                for item in precursors),
            "leftover_atoms": sum(item["leftover_atom_count"]
                                  for item in precursors),
            "formed_bonds": len(formed),
        },
    }


def _assembly_rank(assembly):
    score = assembly["score"]
    return (
        score["chirality_violations"],
        score["unique_precursor_structures"],
        -score["set_atom_retention"],
        score["capped_precursors"], score["broken_bonds"],
        score["leftover_atoms"], score["formed_bonds"],
        tuple(item["precursor_id"] for item in assembly["precursors"]),
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--parts", required=True)
    parser.add_argument("--target-smiles", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--per-mask-limit", type=int, default=200)
    parser.add_argument("--assembly-limit", type=int, default=20)
    parser.add_argument("--maximum-precursors", type=int, default=3)
    parser.add_argument("--per-cover-expansion-limit", type=int, default=5000)
    parser.add_argument("--search-mode",
                        choices=("modular", "recommendation", "exhaustive"),
                        default="modular")
    parser.add_argument("--beam-width", type=int, default=200)
    parser.add_argument("--combination-beam-width", type=int, default=200)
    parser.add_argument("--pattern-limit", type=int, default=4)
    parser.add_argument("--recommendations-per-pattern", type=int, default=4)
    parser.add_argument("--patterns-per-coverage", type=int, default=4)
    parser.add_argument("--modular-state-limit", type=int, default=100000)
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
    full_mask = (1 << atom_count) - 1
    target_edges = [
        (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
        for bond in target.GetBonds()
    ]

    def construction_pattern_signature(pattern):
        masks = tuple(sorted(pattern))
        atom_sets = [
            {atom for atom in range(atom_count) if mask & (1 << atom)}
            for mask in masks
        ]
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

    groups = defaultdict(list)
    expected = defaultdict(list)
    direct_target_matches = []
    counts = Counter()
    part_paths = sorted(Path(args.parts).glob("part_*.jsonl.gz"))
    if not part_paths:
        raise RuntimeError("no completed part files found")
    for path in part_paths:
        with gzip.open(path, "rt", encoding="utf-8") as stream:
            for line in stream:
                record = json.loads(line)
                molecule = Chem.MolFromSmiles(record["smiles"])
                molecule_key = (Chem.MolToSmiles(
                    molecule, isomericSmiles=True)
                                 if molecule is not None
                                 else record["precursor_id"])
                structure_key = molecule_key
                if (args.exclude_target_identity and molecule is not None and
                        molecule_key == target_key):
                    direct_target_matches.append({
                        "precursor_id": record["precursor_id"],
                        "smiles": record["smiles"],
                        "row_index": record["row_index"],
                    })
                    continue
                counts["matched_precursors"] += 1
                counts["capped_precursors"] += record["status"] == "capped"
                explicit_molecule = Chem.AddHs(molecule)
                if args.chirality_ranking:
                    Chem.AssignStereochemistry(
                        explicit_molecule, cleanIt=True, force=True)
                total_heavy_atoms = molecule.GetNumHeavyAtoms()
                candidates = record["candidates"]
                if args.attachment_trim_variants:
                    candidates = [
                        variant
                        for candidate in candidates
                        for variant in _attachment_trim_variants(
                            candidate, explicit_molecule)
                    ]
                for candidate in candidates:
                    item = _entry(
                        record, candidate, explicit_molecule,
                        total_heavy_atoms, structure_key, target,
                        args.chirality_ranking)
                    mask = _mask(item["covered_target_atoms"])
                    bucket = groups[mask]
                    bucket.append(item)
                    if len(bucket) >= args.per_mask_limit * 2:
                        bucket.sort(key=_entry_rank)
                        del bucket[args.per_mask_limit:]
                    counts["fragment_candidates"] += 1
                    if item["precursor_id"] in args.expected_id:
                        expected[item["precursor_id"]].append(item)

    for mask, items in groups.items():
        groups[mask] = sorted(items, key=_entry_rank)[:args.per_mask_limit]

    masks_by_atom = defaultdict(list)
    for mask in groups:
        for atom in range(atom_count):
            if mask & (1 << atom):
                masks_by_atom[atom].append(mask)
    for atom in masks_by_atom:
        masks_by_atom[atom].sort(key=lambda mask: (-mask.bit_count(), mask))
    all_masks = tuple(sorted(groups, key=lambda mask: (-mask.bit_count(), mask)))
    largest_mask_size = max((mask.bit_count() for mask in all_masks), default=0)

    @lru_cache(maxsize=None)
    def can_complete(covered, slots):
        uncovered = full_mask ^ covered
        if not uncovered:
            return True
        compatible = [mask for mask in all_masks if not (mask & covered)]
        if not compatible or slots == 0:
            return False
        possible = 0
        sizes = []
        for mask in compatible:
            possible |= mask
            sizes.append(mask.bit_count())
        if uncovered & ~possible:
            return False
        sizes.sort(reverse=True)
        return uncovered.bit_count() <= sum(sizes[:slots])

    def exact_covers(covered=0, selected=()):
        if covered == full_mask:
            yield selected
            return
        if len(selected) >= args.maximum_precursors:
            return
        if not can_complete(
                covered, args.maximum_precursors - len(selected)):
            return
        uncovered = full_mask ^ covered
        pivot = (uncovered & -uncovered).bit_length() - 1
        for mask in masks_by_atom[pivot]:
            if mask & covered:
                continue
            yield from exact_covers(covered | mask, selected + (mask,))

    def item_set_rank(items, covered_atom_count=0):
        structures = {item["structure_key"] for item in items}
        retained = sum(item["retained_atom_count"] for item in items)
        total = sum(item["total_atom_count"] for item in items)
        set_retention = Fraction(retained, total) if total else Fraction()
        return (
            sum(item.get("chirality_violations", 0) for item in items),
            len(structures),
            -set_retention,
            -covered_atom_count,
            sum(not item["complete"] for item in items),
            sum(len(item["boundary_bonds"]) for item in items),
            sum(item["leftover_atom_count"] for item in items),
            tuple(item["precursor_id"] for item in items),
        )

    beam_truncated = False

    def recommendation_covers():
        nonlocal beam_truncated
        states = [(0, ())]
        completed = []
        for _depth in range(args.maximum_precursors):
            next_states = []
            for covered, selected in states:
                uncovered = full_mask ^ covered
                pivot = (uncovered & -uncovered).bit_length() - 1
                for mask in masks_by_atom[pivot]:
                    if mask & covered:
                        continue
                    new_covered = covered | mask
                    new_selected = selected + (mask,)
                    if new_covered == full_mask:
                        completed.append(new_selected)
                    elif ((full_mask ^ new_covered).bit_count()
                          <= (args.maximum_precursors - len(new_selected))
                          * largest_mask_size):
                        next_states.append((new_covered, new_selected))
            if len(next_states) > args.beam_width:
                beam_truncated = True
            next_states.sort(key=lambda state: item_set_rank(
                [groups[mask][0] for mask in state[1]],
                state[0].bit_count()))
            states = next_states[:args.beam_width]
            if not states:
                break
        if len(completed) > args.beam_width:
            beam_truncated = True
        completed.sort(key=lambda cover: item_set_rank(
            [groups[mask][0] for mask in cover], atom_count))
        return completed[:args.beam_width]

    def modular_covers():
        nonlocal beam_truncated
        states = {0: [()]}
        completed = []
        for _depth in range(args.maximum_precursors):
            next_states = defaultdict(list)
            for covered, paths in states.items():
                uncovered = full_mask ^ covered
                pivot = (uncovered & -uncovered).bit_length() - 1
                for selected in paths:
                    for mask in masks_by_atom[pivot]:
                        if mask & covered:
                            continue
                        new_covered = covered | mask
                        new_selected = selected + (mask,)
                        if new_covered == full_mask:
                            completed.append(new_selected)
                            continue
                        remaining_slots = (
                            args.maximum_precursors - len(new_selected))
                        if ((full_mask ^ new_covered).bit_count()
                                > remaining_slots * largest_mask_size):
                            continue
                        bucket = next_states[new_covered]
                        bucket.append(new_selected)
                        if len(bucket) >= args.patterns_per_coverage * 2:
                            bucket.sort(key=lambda path: item_set_rank(
                                [groups[item][0] for item in path],
                                new_covered.bit_count()))
                            del bucket[args.patterns_per_coverage:]
            for covered, paths in next_states.items():
                paths.sort(key=lambda path: item_set_rank(
                    [groups[item][0] for item in path],
                    covered.bit_count()))
                del paths[args.patterns_per_coverage:]
            state_count = sum(map(len, next_states.values()))
            if state_count > args.modular_state_limit:
                beam_truncated = True
                ranked = sorted(
                    ((covered, path)
                     for covered, paths in next_states.items()
                     for path in paths),
                    key=lambda pair: item_set_rank(
                        [groups[item][0] for item in pair[1]],
                        pair[0].bit_count()),
                )[:args.modular_state_limit]
                next_states = defaultdict(list)
                for covered, path in ranked:
                    next_states[covered].append(path)
            states = dict(next_states)
            if not states:
                break
        return completed

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

    prune_at = max(10_000, args.assembly_limit * 4)
    retain_during_search = max(args.assembly_limit * 2, 5_000)
    if args.search_mode == "modular":
        covers = modular_covers()
    elif args.search_mode == "recommendation":
        covers = recommendation_covers()
    else:
        covers = exact_covers()
    for cover in covers:
        cover_count += 1
        pattern_key = construction_pattern_signature(cover)
        representative_pattern.setdefault(pattern_key, tuple(sorted(cover)))
        pattern_assemblies = []
        pools = [groups[mask] for mask in cover]
        if args.search_mode == "recommendation":
            partials = [()]
            for pool in pools:
                expanded = [partial + (item,)
                            for partial in partials for item in pool]
                if len(expanded) > args.combination_beam_width:
                    beam_truncated = True
                expanded.sort(key=lambda items: item_set_rank(items))
                partials = expanded[:args.combination_beam_width]
            item_combinations = partials
        else:
            item_combinations = itertools.islice(
                itertools.product(*pools), args.per_cover_expansion_limit)
        for items in item_combinations:
            expansion_count += 1
            ids = [item["precursor_id"] for item in items]
            if (not args.allow_repeated_precursors and
                    len(set(ids)) != len(ids)):
                continue
            formed = _formed_bonds(
                items, target_edges, args.require_attachment_bonds)
            if formed is not None:
                pattern_assemblies.append(_assembly(items, formed))
                assembly_mapping_variants += 1
            if len(pattern_assemblies) >= prune_at:
                pattern_assemblies = retain_best(
                    pattern_assemblies, retain_during_search)
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
                [atom for atom in range(atom_count) if mask & (1 << atom)]
                for mask in pattern
            ],
            "fragment_sizes": [mask.bit_count() for mask in pattern],
            "recommendation_ranks": list(range(
                start_rank, start_rank + len(variants))),
            "best_score": variants[0]["score"],
        })
    assemblies = assemblies[:args.assembly_limit]

    expected_assembly = None
    if len(args.expected_id) == 2:
        for left in expected.get(args.expected_id[0], ()):
            left_mask = _mask(left["covered_target_atoms"])
            for right in expected.get(args.expected_id[1], ()):
                if (left_mask | _mask(right["covered_target_atoms"])) != full_mask:
                    continue
                if left_mask & _mask(right["covered_target_atoms"]):
                    continue
                formed = _formed_bonds(
                    (left, right), target_edges,
                    args.require_attachment_bonds)
                if formed is not None:
                    expected_assembly = _assembly((left, right), formed)
                    break
            if expected_assembly:
                break

    summaries = []
    for path in part_paths:
        summary_path = path.with_suffix(path.suffix + ".summary.json")
        if summary_path.exists():
            summaries.append(json.loads(summary_path.read_text()))
    scan_counts = Counter()
    for summary in summaries:
        scan_counts.update(summary["counts"])

    report = {
        "schema": "rxn_core.retro_catalog_assembly/v1",
        "target_smiles": args.target_smiles,
        "part_count": len(part_paths),
        "scan_counts": dict(scan_counts),
        "merge_counts": dict(counts),
        "coverage_masks_retained": len(groups),
        "per_mask_limit": args.per_mask_limit,
        "assembly_limit": args.assembly_limit,
        "maximum_precursors": args.maximum_precursors,
        "exact_cover_count": cover_count,
        "cover_expansion_count": expansion_count,
        "require_attachment_bonds": args.require_attachment_bonds,
        "exclude_target_identity": args.exclude_target_identity,
        "allow_repeated_precursors": args.allow_repeated_precursors,
        "attachment_trim_variants": args.attachment_trim_variants,
        "chirality_ranking": args.chirality_ranking,
        "search_mode": args.search_mode,
        "beam_width": args.beam_width,
        "combination_beam_width": args.combination_beam_width,
        "recommendation_search_truncated": beam_truncated,
        "pattern_limit": args.pattern_limit,
        "recommendations_per_pattern": args.recommendations_per_pattern,
        "patterns_per_coverage": args.patterns_per_coverage,
        "modular_state_limit": args.modular_state_limit,
        "construction_patterns": construction_patterns,
        "direct_target_matches": direct_target_matches,
        "assembly_mapping_variants": assembly_mapping_variants,
        "expected_ids": args.expected_id,
        "expected_ids_found": {
            item: bool(expected.get(item)) for item in args.expected_id
        },
        "expected_pair_recovered": expected_assembly is not None,
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
        "expected_pair_recovered": report["expected_pair_recovered"],
    }, indent=2))


if __name__ == "__main__":
    main()
