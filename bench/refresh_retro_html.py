"""Refresh viewer layout without changing saved mappings, coordinates, or scores."""
import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'tools'))
from build_retro_db_viewer import _html


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('html', type=Path, nargs='+')
    args = parser.parse_args()
    for path in args.html:
        saved = path.read_text().split('<script>const data=', 1)[1]
        payload, _ = json.JSONDecoder().raw_decode(saved)
        updated = _html(payload)
        decoded, _ = json.JSONDecoder().raw_decode(updated.split('<script>const data=', 1)[1])
        assert decoded == payload
        path.write_text(updated)
        print(path)


if __name__ == '__main__':
    main()
