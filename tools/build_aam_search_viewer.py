"""Create an offline graph/path viewer from saved AAM, without rerunning it."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from rxn_core import aam_from_record, write_aam_bundle


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('aam_json', type=Path)
    parser.add_argument('output_directory', type=Path)
    args = parser.parse_args()
    result = aam_from_record(json.loads(args.aam_json.read_text()))
    output = write_aam_bundle(result, args.output_directory)
    print((output / 'search.html').resolve())


if __name__ == '__main__':
    main()
