"""Refresh viewer layout without changing saved mappings, coordinates, or scores."""
import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'tools'))
from build_retro_db_viewer import _html, _fragment_colors, _group_precursors


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('html', type=Path, nargs='+')
    parser.add_argument('--results', type=Path,
                        help='Saved viewer report supplying exact fragment partitions')
    args = parser.parse_args()
    report = json.loads(args.results.read_text()) if args.results else None
    for path in args.html:
        saved = path.read_text().split('<script>const data=', 1)[1]
        payload, _ = json.JSONDecoder().raw_decode(saved)
        if report is not None:
            for assembly in payload['assemblies']:
                rank = assembly['rank']
                if isinstance(rank, int):
                    raw = report['assemblies'][rank - 1]
                elif rank.startswith('validation '):
                    raw = report['validation_assemblies'][int(rank.split()[1]) - 1]
                elif rank == 'diagnostic':
                    raw = report['diagnostic_assembly']
                elif rank == 'ground truth':
                    raw = report['expected_mapping']
                else:
                    raise ValueError(f'Unknown saved assembly rank: {rank}')
                assembly['fragments'] = _fragment_colors(
                    _group_precursors(raw['precursors']), assembly['models'])
            if payload['unassembled_target'] is not None:
                payload['unassembled_target'].update(fragmentStyles=[], fragmentAlternatives=[])
        updated = _html(payload)
        decoded, _ = json.JSONDecoder().raw_decode(updated.split('<script>const data=', 1)[1])
        assert decoded == payload
        path.write_text(updated)
        print(path)


if __name__ == '__main__':
    main()
