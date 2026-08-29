#!/usr/bin/env python3
"""Blind one-precursor rearrangement search using explicit-H no-cut AAM."""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import multiprocessing as mp
import time
from collections import Counter
from pathlib import Path

import numpy as np
from rdkit import Chem, RDLogger

from rxn_core.alignment.api import match_wbo_graphs


_TARGET_ELEMENTS = None
_TARGET_WBO = None
_TARGET_COUNTS = None
_MAX_EXCESS = None
_MIN_EXCESS = None
_MAX_EVENTS = None
_BRANCH_LIMIT = None
_TARGET_MOLECULE = None
_USE_CHIRALITY = None
_REQUIRE_ZERO_CHIRALITY = None


def _mol_graph(molecule):
    elements = [atom.GetSymbol() for atom in molecule.GetAtoms()]
    wbo = np.zeros((len(elements), len(elements)), dtype=float)
    for bond in molecule.GetBonds():
        left, right = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        value = float(bond.GetBondTypeAsDouble()) or 1.0
        wbo[left, right] = wbo[right, left] = value
    return elements, wbo


def _worker_init(target_smiles, min_excess, max_excess,
                 max_events, branch_limit, use_chirality,
                 require_zero_chirality):
    global _TARGET_ELEMENTS, _TARGET_WBO, _TARGET_COUNTS
    global _MIN_EXCESS, _MAX_EXCESS, _MAX_EVENTS, _BRANCH_LIMIT
    global _TARGET_MOLECULE, _USE_CHIRALITY, _REQUIRE_ZERO_CHIRALITY
    RDLogger.DisableLog("rdApp.*")
    target_implicit = Chem.MolFromSmiles(target_smiles)
    if target_implicit is None:
        raise ValueError("invalid target SMILES")
    target = Chem.AddHs(target_implicit)
    _TARGET_ELEMENTS, _TARGET_WBO = _mol_graph(target)
    Chem.AssignStereochemistry(target, cleanIt=True, force=True)
    _TARGET_MOLECULE = target
    _TARGET_COUNTS = Counter(_TARGET_ELEMENTS)
    _MIN_EXCESS = int(min_excess)
    _MAX_EXCESS = int(max_excess)
    _MAX_EVENTS = int(max_events)
    _BRANCH_LIMIT = int(branch_limit)
    _USE_CHIRALITY = bool(use_chirality)
    _REQUIRE_ZERO_CHIRALITY = bool(require_zero_chirality)


def _cip(atom):
    return atom.GetProp("_CIPCode") if atom.HasProp("_CIPCode") else None


def _chirality_violations(source, mapping, broken, formed):
    """Count changed/undefined CIP labels at spectator stereocentres."""
    bad_source = {int(atom) for event in broken for atom in event[:2]}
    inverse = {int(image): int(atom) for atom, image in mapping.items()}
    for event in formed:
        for image in event[:2]:
            if int(image) in inverse:
                bad_source.add(inverse[int(image)])
    violations = 0
    for source_index, target_index in mapping.items():
        if int(source_index) in bad_source:
            continue
        source_cip = _cip(source.GetAtomWithIdx(int(source_index)))
        target_cip = _cip(_TARGET_MOLECULE.GetAtomWithIdx(int(target_index)))
        if (source_cip is not None or target_cip is not None) and source_cip != target_cip:
            violations += 1
    return violations


def _event_record(event):
    return [int(event[0]), int(event[1]), float(event[2]), float(event[3])]


def _search_batch(batch):
    counts = Counter(rows=len(batch))
    records = []
    for row_index, smiles, precursor_id in batch:
        molecule_implicit = Chem.MolFromSmiles(smiles)
        if molecule_implicit is None:
            counts["parse_errors"] += 1
            continue
        molecule = Chem.AddHs(molecule_implicit)
        elements, wbo_R = _mol_graph(molecule)
        if _USE_CHIRALITY:
            Chem.AssignStereochemistry(molecule, cleanIt=True, force=True)
        source_counts = Counter(elements)
        if any(source_counts[element] < count
               for element, count in _TARGET_COUNTS.items()):
            counts["composition_filtered"] += 1
            continue
        excess = source_counts - _TARGET_COUNTS
        excess_elements = sorted(excess.elements())
        if not _MIN_EXCESS <= len(excess_elements) <= _MAX_EXCESS:
            counts["composition_filtered"] += 1
            continue
        if len(elements) != len(_TARGET_ELEMENTS) + len(excess_elements):
            counts["composition_filtered"] += 1
            continue

        augmented_elements = list(_TARGET_ELEMENTS) + excess_elements
        wbo_P = np.zeros((len(augmented_elements), len(augmented_elements)))
        n_target = len(_TARGET_ELEMENTS)
        wbo_P[:n_target, :n_target] = _TARGET_WBO
        profile = []
        counts["searched"] += 1
        result = match_wbo_graphs(
            elements, wbo_R, augmented_elements, wbo_P,
            graph_floor=0.2,
            iso_tol=1.0,
            dwbo_threshold=0.5,
            n_seeds=3,
            max_branches=_BRANCH_LIMIT,
            cut_edges=(),
            repair_symmetry=True,
            chirality=False,
            profile=profile,
        )
        cap_hits = sum(item.get("result") == "subtree_branch_cap"
                       for item in profile)
        counts["cap_hits"] += cap_hits
        if result.best is None or result.best.score[0] > _MAX_EVENTS:
            counts["event_filtered"] += 1
            continue
        if _USE_CHIRALITY:
            scored = [
                ((candidate.score[0], _chirality_violations(
                    molecule, candidate.mapping, candidate.broken,
                    candidate.formed)), candidate)
                for candidate in result.candidates
                if candidate.score[0] <= _MAX_EVENTS
            ]
            if _REQUIRE_ZERO_CHIRALITY:
                scored = [item for item in scored if item[0][1] == 0]
                if not scored:
                    counts["chirality_filtered"] += 1
                    continue
            score, best = min(scored, key=lambda item: item[0])
        else:
            best = result.best
            score = best.score
        counts["matched_precursors"] += 1
        records.append({
            "schema": "rxn_core.rearrangement_search/v1",
            "row_index": row_index,
            "precursor_id": precursor_id,
            "smiles": smiles,
            "excess_elements": excess_elements,
            "candidate_count": len(result.candidates),
            "cap_hits": cap_hits,
            "complete": cap_hits == 0,
            "score": list(score),
            "chirality_violations": int(score[1]),
            "mapping": sorted([int(a), int(b)] for a, b in best.mapping.items()),
            "raw_mapping": sorted(
                [int(a), int(b)] for a, b in best.raw_mapping.items()),
            "broken": [_event_record(item) for item in best.broken],
            "formed": [_event_record(item) for item in best.formed],
        })
    return dict(counts), records


def _batches(path, shard_index, shard_count, batch_size, id_column):
    batch = []
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as stream:
        for row_index, row in enumerate(csv.DictReader(stream)):
            if row_index % shard_count != shard_index:
                continue
            smiles, precursor_id = row.get("SMILES"), row.get(id_column)
            if not smiles or not precursor_id:
                continue
            batch.append((row_index, smiles, precursor_id))
            if len(batch) >= batch_size:
                yield batch
                batch = []
    if batch:
        yield batch


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--target-smiles", required=True)
    parser.add_argument("--catalog", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--shard-index", type=int, default=0)
    parser.add_argument("--shard-count", type=int, default=1)
    parser.add_argument("--batch-size", type=int, default=32)
    parser.add_argument("--id-column", default="Mcule ID")
    parser.add_argument("--minimum-excess-atoms", type=int, default=1)
    parser.add_argument("--maximum-excess-atoms", type=int, default=3)
    parser.add_argument("--maximum-events", type=int, default=6)
    parser.add_argument("--branch-limit", type=int, default=100)
    parser.add_argument("--chirality", action=argparse.BooleanOptionalAction,
                        default=False)
    parser.add_argument("--require-zero-chirality",
                        action=argparse.BooleanOptionalAction, default=False)
    args = parser.parse_args()
    if args.require_zero_chirality and not args.chirality:
        parser.error("--require-zero-chirality requires --chirality")

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    totals = Counter()
    started = time.perf_counter()
    context = mp.get_context("fork")
    with gzip.open(output, "wt", encoding="utf-8") as sink:
        with context.Pool(
            processes=max(1, args.workers),
            initializer=_worker_init,
            initargs=(args.target_smiles, args.minimum_excess_atoms,
                      args.maximum_excess_atoms,
                      args.maximum_events, args.branch_limit,
                      args.chirality, args.require_zero_chirality),
        ) as pool:
            work = pool.imap_unordered(
                _search_batch,
                _batches(args.catalog, args.shard_index, args.shard_count,
                         args.batch_size, args.id_column),
                chunksize=1,
            )
            for counts, records in work:
                totals.update(counts)
                for record in records:
                    sink.write(json.dumps(record, separators=(",", ":")) + "\n")
    elapsed = time.perf_counter() - started
    summary = {
        "schema": "rxn_core.rearrangement_summary/v1",
        "target_smiles": args.target_smiles,
        "shard_index": args.shard_index,
        "shard_count": args.shard_count,
        "workers": args.workers,
        "branch_limit": args.branch_limit,
        "maximum_events": args.maximum_events,
        "minimum_excess_atoms": args.minimum_excess_atoms,
        "maximum_excess_atoms": args.maximum_excess_atoms,
        "chirality": args.chirality,
        "require_zero_chirality": args.require_zero_chirality,
        "explicit_hydrogens": True,
        "elapsed_seconds": elapsed,
        "counts": dict(totals),
    }
    output.with_suffix(output.suffix + ".summary.json").write_text(
        json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
