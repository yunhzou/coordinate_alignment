#!/usr/bin/env python3
"""Build the all-case R.xyz/P_final.xyz package and offline viewer."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

if __package__:
    from .neb_support.neb_orientation_package import build_deliverable
else:
    from neb_support.neb_orientation_package import build_deliverable


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path, required=True)
    parser.add_argument("--orientation-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument(
        "--3dmol-js", dest="renderer_js", type=Path, required=True)
    parser.add_argument(
        "--3dmol-license", dest="renderer_license", type=Path)
    parser.add_argument("--archive", type=Path)
    parser.add_argument(
        "--interpolation-images",
        type=int,
        help=(
            "also generate this many endpoint-inclusive IDPP internal-"
            "coordinate path images and add piecewise path diagnostics/"
            "viewer controls"
        ),
    )
    return parser


def main(argv=None) -> int:
    args = _parser().parse_args(argv)
    result = build_deliverable(
        args.source_root,
        args.orientation_root,
        args.output_root,
        args.renderer_js,
        renderer_license_path=args.renderer_license,
        archive_path=args.archive,
        interpolation_image_count=args.interpolation_images,
    )
    print(json.dumps({
        "output_root": str(args.output_root.resolve()),
        "archive": (
            None if args.archive is None else str(args.archive.resolve())),
        "case_count": result["case_count"],
        "mechanism_count": result["mechanism_count"],
        "viewer": result["viewer"],
    }, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
