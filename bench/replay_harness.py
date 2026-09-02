"""Layer-0 replay harness: run search_aam serially on an example, record call
counts, a running digest of every extension / dedupe / growth result, and the
full completed-branch pool.  Compare two runs with --compare.

Usage:
  .venv/bin/python bench/replay_harness.py record  <out.json> [--case tempo]
  .venv/bin/python bench/replay_harness.py compare <a.json> <b.json>
"""
from __future__ import annotations

import hashlib
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import pynauty  # noqa: E402

import rxn_core.alignment.branch as branch_mod  # noqa: E402
import rxn_core.growth.island as island_mod  # noqa: E402
import rxn_core.matcher.extend as extend_mod  # noqa: E402
from rxn_core.aam import search_aam  # noqa: E402
from rxn_core.cli import _endpoint_cache  # noqa: E402
from rxn_core.domain import AAMProblem, AAMSearchConfig  # noqa: E402
from rxn_core.matcher.state import _SymCand  # noqa: E402

CASES = {
    "tempo": (
        ROOT / "docs/example_runs/pr1.tempo_ts3/work/endpoints/R",
        ROOT / "docs/example_runs/pr1.tempo_ts3/work/endpoints/P",
    ),
}


def _cand_key(cand):
    if isinstance(cand, _SymCand):
        return (
            tuple(sorted(cand.mapping.items())),
            tuple((b.r_atoms, b.p_atoms, bool(b.extendable)) for b in cand.blocks),
            tuple(sorted(cand.exact_fixed)),
            int(cand.multiplicity),
            tuple((b.r_atoms, b.p_atoms) for b in cand.automorph_blocks),
        )
    return ("dict", tuple(sorted(dict(cand).items())))


def _iso_key(iso):
    return (
        tuple(sorted(iso.items())),
        tuple(sorted(tuple(e) for e in getattr(iso, "deferred_edges", ()))),
        tuple(sorted(getattr(iso, "fragment", ()))),
        json.dumps(getattr(iso, "symmetry", {}), sort_keys=True, default=list),
    )


class Recorder:
    def __init__(self):
        self.counts = {
            "extend_calls": 0, "children_returned": 0, "single_child_steps": 0,
            "zero_child_steps": 0, "dedupe_calls": 0, "dedupe_inputs": 0,
            "dedupe_single_input": 0, "certificate_calls": 0, "autgrp_calls": 0,
            "grow_calls": 0, "grow_isos": 0, "find_islands_calls": 0,
        }
        self.h_extend = hashlib.blake2b(digest_size=16)
        self.h_grow = hashlib.blake2b(digest_size=16)
        self.h_find = hashlib.blake2b(digest_size=16)
        self.t_extend = 0.0
        self.t_dedupe = 0.0

    def install(self):
        rec = self
        orig_extend = extend_mod._extend_sym_cands
        orig_dedupe = extend_mod._dedup_sym_cands
        orig_grow = island_mod.grow_island
        orig_find = branch_mod.find_islands
        orig_cert = pynauty.certificate
        orig_autgrp = pynauty.autgrp

        def extend(*a, **k):
            t0 = time.perf_counter()
            out = orig_extend(*a, **k)
            rec.t_extend += time.perf_counter() - t0
            rec.counts["extend_calls"] += 1
            rec.counts["children_returned"] += len(out)
            if len(out) == 1:
                rec.counts["single_child_steps"] += 1
            if not out:
                rec.counts["zero_child_steps"] += 1
            rec.h_extend.update(repr([_cand_key(c) for c in out]).encode())
            return out

        def dedupe(cands, *a, **k):
            cands = list(cands)
            t0 = time.perf_counter()
            out = orig_dedupe(cands, *a, **k)
            rec.t_dedupe += time.perf_counter() - t0
            rec.counts["dedupe_calls"] += 1
            rec.counts["dedupe_inputs"] += len(cands)
            if len(cands) == 1:
                rec.counts["dedupe_single_input"] += 1
            return out

        def grow(*a, **k):
            out = orig_grow(*a, **k)
            rec.counts["grow_calls"] += 1
            rec.counts["grow_isos"] += len(out)
            rec.h_grow.update(repr([_iso_key(i) for i in out]).encode())
            return out

        def find(*a, **k):
            out = orig_find(*a, **k)
            rec.counts["find_islands_calls"] += 1
            rec.h_find.update(repr([
                (tuple(sorted(b.mapping.items())),
                 tuple(sorted(b.islands_R.items())),
                 tuple(sorted(tuple(e) for e in b.deferred_edges)),
                 json.dumps(b.symmetry_paths, sort_keys=True, default=list))
                for b in out]).encode())
            return out

        def cert(*a, **k):
            rec.counts["certificate_calls"] += 1
            return orig_cert(*a, **k)

        def autgrp(*a, **k):
            rec.counts["autgrp_calls"] += 1
            return orig_autgrp(*a, **k)

        # the growth loop imported these names into its own namespace
        extend_mod._extend_sym_cands = extend
        island_mod._extend_sym_cands = extend
        extend_mod._dedup_sym_cands = dedupe
        island_mod._dedup_sym_cands = dedupe
        island_mod.grow_island = grow
        branch_mod.grow_island = grow
        branch_mod.find_islands = find
        import rxn_core.alignment.sweep as sweep_mod
        sweep_mod.find_islands = find
        pynauty.certificate = cert
        pynauty.autgrp = autgrp


def pool_record(result):
    mechanisms = []
    for mech in result.mechanisms:
        branches = []
        for br in mech.branches:
            branches.append({
                "mapping": list(br.representative.images),
                "hierarchy": br.hierarchy.to_record(),
                "encounter_count": br.encounter_count,
                "cuts": [list(c) for c in br.cuts],
                "covered_path_count": br.covered_path_count,
            })
        mechanisms.append({
            "key": json.dumps(mech.key, default=list),
            "representative": list(mech.representative.images),
            "cuts": [list(c) for c in mech.cuts],
            "includes_uncut_search": mech.includes_uncut_search,
            "encounter_count": mech.encounter_count,
            "branches": branches,
        })
    return mechanisms


def record(out_path, case="tempo", workers=1):
    r_dir, p_dir = CASES[case]
    reactant = _endpoint_cache(str(r_dir), "R")
    product = _endpoint_cache(str(p_dir), "P")
    problem = AAMProblem(reactant, product, name=case)
    rec = Recorder()
    rec.install()
    t0 = time.perf_counter()
    result = search_aam(problem, AAMSearchConfig(), workers=workers)
    elapsed = time.perf_counter() - t0
    payload = {
        "case": case,
        "workers": workers,
        "elapsed_seconds": elapsed,
        "extend_seconds": rec.t_extend,
        "dedupe_seconds_inside_extend": rec.t_dedupe,
        "counts": rec.counts,
        "digest_extend": rec.h_extend.hexdigest(),
        "digest_grow": rec.h_grow.hexdigest(),
        "digest_find_islands": rec.h_find.hexdigest(),
        "metrics": result.metrics.__dict__,
        "pool": pool_record(result),
    }
    Path(out_path).write_text(json.dumps(payload, indent=1, default=list))
    summary = {k: v for k, v in payload.items() if k != "pool"}
    summary["mechanisms"] = len(payload["pool"])
    summary["branches"] = sum(len(m["branches"]) for m in payload["pool"])
    print(json.dumps(summary, indent=1, default=str))


def compare(a_path, b_path):
    a = json.loads(Path(a_path).read_text())
    b = json.loads(Path(b_path).read_text())
    same_pool = a["pool"] == b["pool"]
    print("pool identical           :", same_pool)
    for key in ("digest_extend", "digest_grow", "digest_find_islands"):
        print(f"{key:25s}:", "identical" if a[key] == b[key] else "DIFFERENT")
    print(f"elapsed  {a['elapsed_seconds']:.2f} s -> {b['elapsed_seconds']:.2f} s  "
          f"(x{a['elapsed_seconds']/max(b['elapsed_seconds'],1e-9):.2f})")
    print(f"extend   {a['extend_seconds']:.2f} s -> {b['extend_seconds']:.2f} s")
    print(f"dedupe   {a['dedupe_seconds_inside_extend']:.2f} s -> {b['dedupe_seconds_inside_extend']:.2f} s")
    for k in a["counts"]:
        if a["counts"][k] != b["counts"].get(k):
            print(f"  count {k:22s} {a['counts'][k]} -> {b['counts'].get(k)}")
    if not same_pool:
        # locate first difference
        for i, (ma, mb) in enumerate(zip(a["pool"], b["pool"])):
            if ma != mb:
                print("first differing mechanism index", i, "key", ma["key"][:80])
                break
        sys.exit(1)


if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "record":
        case = "tempo"
        workers = 1
        args = sys.argv[3:]
        if "--case" in args:
            case = args[args.index("--case") + 1]
        if "--workers" in args:
            workers = int(args[args.index("--workers") + 1])
        record(sys.argv[2], case=case, workers=workers)
    elif mode == "compare":
        compare(sys.argv[2], sys.argv[3])
    else:
        raise SystemExit(__doc__)
