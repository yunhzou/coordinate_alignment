"""Regenerate out/mode_analysis/<step>.csv from the patched per-step
HTMLs. Pure data extraction — no alignment, no xtb."""
from __future__ import annotations
import sys as _sys
from pathlib import Path as _Path
_HERE = _Path(__file__).resolve().parent
_sys.path.insert(0, str(_HERE.parent / "src"))
_sys.path.insert(0, str(_HERE))
PROJECT_ROOT = _HERE.parent  # _RXN_CORE_PATH_SETUP

import csv
import json
import re
from pathlib import Path

SRC = PROJECT_ROOT / "out" / "mode_viewer"
DST = PROJECT_ROOT / "out" / "mode_analysis"
DST.mkdir(exist_ok=True)


def main():
    files = sorted(SRC.glob("*.html"))
    files = [f for f in files if f.name not in ('index.html', 'flat_view.html')]
    print(f"Regenerating {len(files)} CSVs from patched HTMLs...")
    n_ok = 0
    fieldnames = ['step', 'ts_label', 'mode_idx', 'freq',
                  'bond_overlap', 'rxn_overlap', 'core_fraction',
                  'mode_rank', 'n_imag', 'n_modes_total', 'n_core_atoms',
                  'core_atoms']
    for hp in files:
        text = hp.read_text()
        m = re.search(r"const DATA = (\{.*?\});\n", text, re.DOTALL)
        if not m: continue
        data = json.loads(m.group(1))
        rows = []
        for ts in data['ts_list']:
            imag = [(i, mm) for i, mm in enumerate(ts['modes']) if mm['freq'] < 0]
            imag.sort(key=lambda t: -t[1].get('bond_overlap', 0))
            for rank, (_, mm) in enumerate(imag):
                rows.append({
                    'step': data['step'],
                    'ts_label': ts['label'],
                    'mode_idx': mm['idx'],
                    'freq': mm['freq'],
                    'bond_overlap': mm.get('bond_overlap', 0),
                    'rxn_overlap': mm.get('rxn_overlap', 0),
                    'core_fraction': mm.get('core_fraction', 0),
                    'mode_rank': rank,
                    'n_imag': ts['n_imag'],
                    'n_modes_total': ts['n_modes_total'],
                    'n_core_atoms': len(data['core_atoms']),
                    'core_atoms': ','.join(map(str, data['core_atoms'])),
                })
        sanitized = re.sub(r"[^A-Za-z0-9._-]", "_", data['step'])
        with (DST / f"{sanitized}.csv").open('w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=fieldnames)
            w.writeheader(); w.writerows(rows)
        n_ok += 1
    print(f"Done. {n_ok} CSVs regenerated at {DST}")


if __name__ == "__main__":
    main()
