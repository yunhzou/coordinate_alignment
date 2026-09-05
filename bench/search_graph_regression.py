"""Persist reusable semantic/timing snapshots for the search-graph migration."""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
from cases import CASES


def snapshot(result):
    return sorted(({
        "key": repr(item.key),
        "representative": list(item.representative.images),
        "branches": sorted(({
            "mapping": list(branch.representative.images),
            "fragments": [
                {"atoms": list(fragment.r_atoms),
                 "mapping": list(fragment.representative_assignments),
                 "deferred": list(fragment.deferred_edges)}
                for fragment in branch.hierarchy.fragments
            ],
        } for branch in item.branches), key=lambda x: json.dumps(x, sort_keys=True)),
    } for item in result.mechanisms), key=lambda x: x["key"])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", default="tetraphenyl", choices=CASES)
    parser.add_argument("--output", required=True)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--package-root", type=Path, default=ROOT)
    parser.add_argument("--legacy-api", action="store_true")
    parser.add_argument("--repeats", type=int, default=1)
    parser.add_argument("--save-intermediates", action="store_true")
    args = parser.parse_args()
    sys.path.insert(0, str(args.package_root / "src"))
    from rxn_core import AAMProblem, AAMSearchConfig, search_aam
    if not args.legacy_api:
        from rxn_core.mechanisms import group_mechanisms
    reactant, product = CASES[args.case]()
    timings = []
    for _ in range(args.repeats):
        options = {}
        if args.save_intermediates:
            options['intermediate_dir'] = Path(args.output).parent / (Path(args.output).stem + '_chunks')
        start = time.perf_counter()
        result = search_aam(AAMProblem(reactant, product, args.case),
                            AAMSearchConfig(), workers=args.workers, **options)
        grouped = result if args.legacy_api else group_mechanisms(result)
        timings.append(time.perf_counter() - start)
    report = {"case": args.case, "workers": args.workers,
              "seconds": min(timings), "all_seconds": timings,
              "mechanisms": snapshot(grouped), "metrics": vars(result.metrics)}
    if not args.legacy_api:
        report['graph'] = {'states': len(result.graph.states),
                           'transitions': len(result.graph.transitions),
                           'json_bytes': len(json.dumps(result.graph.to_record())),
                           'grouping_seconds': grouped.elapsed_seconds}
    path = Path(args.output)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(report, indent=2) + "\n")
    if not args.legacy_api:
        from rxn_core.artifacts import aam_record
        path.with_suffix('.aam.json').write_text(json.dumps(aam_record(result)) + '\n')
    print(json.dumps({k: v for k, v in report.items() if k != "mechanisms"}))


if __name__ == "__main__":
    main()
