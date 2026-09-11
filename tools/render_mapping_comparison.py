"""Render saved comparison JSON through the shared original R/P/TS presentation."""
import argparse
import json
from pathlib import Path

from rxn_core.viewers import collection_html, comparison_document


def render(source, output):
    cases = json.loads(Path(source).read_text())
    documents = []
    for case in cases:
        if 'recovery' in case:
            case['note'] = ('R → P; sweep, cap 1000. Missing-pattern recovery: ' + '; '.join(
                f"{r['seeds']} seeds: {'recovered' if r['status'] == 'represented' else 'absent from saved families'}"
                for r in case['recovery']) + '. One-seed AAM alternatives and actual recovered witnesses are separate choices.')
        documents.append(comparison_document(case))
    Path(output).write_text(collection_html(documents))
    return len(documents)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    print(f'Rendered {render(args.source, args.output)} saved cases; no mapping search.')
