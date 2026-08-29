#!/usr/bin/env python3
"""Check explicit-H AAM sensitivity to the number of seed orders."""
from __future__ import annotations

import argparse
import json
import time

from rdkit import Chem

import search_mcule_rearrangement as search
from rxn_core.alignment.api import match_wbo_graphs


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-smiles", required=True)
    parser.add_argument("--target-smiles", required=True)
    parser.add_argument("--n-seeds", type=int, required=True)
    parser.add_argument("--branch-limit", type=int, default=100)
    args = parser.parse_args()

    search._worker_init(
        args.target_smiles, 0, 0, 10_000, args.branch_limit, True, True)
    source = Chem.AddHs(Chem.MolFromSmiles(args.source_smiles))
    elements, wbo_source = search._mol_graph(source)
    profile = []
    started = time.perf_counter()
    result = match_wbo_graphs(
        elements, wbo_source, search._TARGET_ELEMENTS, search._TARGET_WBO,
        graph_floor=0.2, iso_tol=1.0, dwbo_threshold=0.5,
        n_seeds=args.n_seeds, max_branches=args.branch_limit,
        cut_edges=(), repair_symmetry=True, chirality=False,
        profile=profile,
    )
    scored = []
    for candidate in result.candidates:
        chirality = search._chirality_violations(
            source, candidate.mapping, candidate.broken, candidate.formed)
        if chirality == 0:
            scored.append((len(candidate.broken) + len(candidate.formed),
                           candidate))
    score, best = min(scored, key=lambda item: item[0])
    print(json.dumps({
        "n_seeds": args.n_seeds,
        "branch_limit": args.branch_limit,
        "elapsed_seconds": time.perf_counter() - started,
        "candidate_count": len(result.candidates),
        "cap_hits": sum(item.get("result") == "subtree_branch_cap"
                        for item in profile),
        "best_zero_chirality_edits": score,
        "broken": len(best.broken),
        "formed": len(best.formed),
    }, indent=2))


if __name__ == "__main__":
    main()
