"""Command-line composition of the typed rxn_core API."""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from .artifacts import write_rp_bundle, write_ts_record
from .chemistry_computations import load_cached_xtb
from .domain import (
    AAMProblem, AAMSearchConfig, MolecularEndpoint, TransitionStateTarget,
    VibrationalModes,
)
from .rp import align_reaction
from .ts import analyze_transition_state


def _endpoint_npz(path, label):
    with np.load(path, allow_pickle=False) as data:
        return MolecularEndpoint(
            tuple(str(item) for item in data["elements"]),
            data["coordinates"], data["wbo"], label=label)


def _endpoint_cache(path, label):
    elements, coordinates, wbo, _xyz = load_cached_xtb(path)
    return MolecularEndpoint(tuple(elements), coordinates, wbo, label=label)


def _endpoint(args, prefix, label):
    npz = getattr(args, f"{prefix}_npz")
    cache = getattr(args, f"{prefix}_cache")
    if bool(npz) == bool(cache):
        raise SystemExit(
            f"provide exactly one --{prefix}-npz or --{prefix}-cache")
    return _endpoint_npz(npz, label) if npz else _endpoint_cache(cache, label)


def _config(args):
    anchors = tuple(tuple(map(int, value.split(":"))) for value in args.anchor)
    return AAMSearchConfig(
        cut_floor=args.cut_floor,
        graph_floor=args.graph_floor,
        iso_tolerance=args.iso_tolerance,
        event_threshold=args.event_threshold,
        metal_event_threshold=args.metal_event_threshold,
        seed_count=args.seed_count,
        branch_limit=args.branch_limit,
        anchors=anchors,
    )


def _parser():
    parser = argparse.ArgumentParser(prog="rxn-core")
    parser.add_argument("--reactant-npz")
    parser.add_argument("--product-npz")
    parser.add_argument("--reactant-cache")
    parser.add_argument("--product-cache")
    parser.add_argument("--target-npz", action="append", default=[])
    parser.add_argument("--name", default="reaction")
    parser.add_argument("--output", default="rxn_core_output")
    parser.add_argument("--workers", type=int, default=1,
                        help="CPU processes used by the cut sweep")
    parser.add_argument("--post-workers", type=int, default=None,
                        help="CPU processes used for analytical family compilation")
    parser.add_argument("--cut-floor", type=float, default=0.2)
    parser.add_argument("--graph-floor", type=float, default=0.2)
    parser.add_argument("--iso-tolerance", type=float, default=1.0)
    parser.add_argument("--event-threshold", type=float, default=0.5)
    parser.add_argument("--metal-event-threshold", type=float, default=0.3)
    parser.add_argument("--seed-count", type=int, default=3)
    parser.add_argument("--branch-limit", type=int, default=100)
    parser.add_argument("--anchor", action="append", default=[], metavar="R:P")
    return parser


def main(argv=None):
    args = _parser().parse_args(argv)
    reactant = _endpoint(args, "reactant", "R")
    product = _endpoint(args, "product", "P")
    config = _config(args)
    rp = align_reaction(
        AAMProblem(reactant, product, name=args.name),
        search_config=config,
        workers=max(1, args.workers),
        post_workers=args.post_workers)
    output = write_rp_bundle(rp, args.output)
    for index, path in enumerate(args.target_npz, 1):
        with np.load(path, allow_pickle=False) as data:
            target = TransitionStateTarget(
                MolecularEndpoint(
                    tuple(str(item) for item in data["elements"]),
                    data["coordinates"], data["wbo"],
                    label=Path(path).stem),
                VibrationalModes(data["frequencies"], data["modes"]))
        ts = analyze_transition_state(rp, target, search_config=config)
        write_ts_record(ts, output / f"ts_{index:03d}.json")
    print(output.resolve())


if __name__ == "__main__":
    main()
