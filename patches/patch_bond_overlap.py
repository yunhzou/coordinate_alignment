"""
Patch existing per-step viewer HTMLs in-place to add `bond_overlap` to
every mode and re-rank the default selection by it.

NO xtb invocation, NO PQ alignment recompute. The per-step HTMLs at
out/mode_viewer/<step>.html already contain (after the rxn_overlap
build):
  - panels[i].xyz_coords       (TS coords reindexed to R-frame)
  - panels[i].modes[k].disp    (mode displacement reindexed to R-frame)
  - broken_bonds   (R-frame indices)
  - formed_bonds_R (R-frame indices)
  - core_atoms     (R-frame indices)

That's everything needed to compute bond_overlap = |d · V̂| / ||d||
where V is the bond-stretch (broken) / bond-compress (formed) vector
at TS coords. We just read each HTML, compute, splice back in.
"""
from __future__ import annotations
import sys as _sys
from pathlib import Path as _Path
_HERE = _Path(__file__).resolve().parent
_sys.path.insert(0, str(_HERE.parent / "src"))
_sys.path.insert(0, str(_HERE))
PROJECT_ROOT = _HERE.parent  # _RXN_CORE_PATH_SETUP

import json
import re
import sys
import time
from pathlib import Path

import numpy as np

from analyze_core_modes import bond_reaction_vector, bond_overlap_per_mode

OUT_DIR = PROJECT_ROOT / "out" / "mode_viewer"


def patch_step(html_path):
    text = html_path.read_text()
    match = re.search(r"const DATA = (\{.*?\});\n", text, re.DOTALL)
    if not match:
        return False
    data = json.loads(match.group(1))

    broken = [tuple(b) for b in data['broken_bonds']]
    formed = [tuple(b) for b in data['formed_bonds_R']]

    for ts in data['ts_list']:
        if not ts['modes']:
            continue
        xyz = np.asarray(ts['xyz_coords'])
        V = bond_reaction_vector(xyz, broken, formed)
        disps = np.asarray([mm['disp'] for mm in ts['modes']])  # (k, n_atoms, 3)
        bond_ov = bond_overlap_per_mode(disps, V)               # (k,)
        for i, mm in enumerate(ts['modes']):
            mm['bond_overlap'] = round(float(bond_ov[i]), 4)
        # Default = imag mode with highest bond_overlap (else any).
        imag_positions = [i for i, mm in enumerate(ts['modes']) if mm['freq'] < 0]
        if imag_positions:
            best_pos = max(imag_positions, key=lambda i: ts['modes'][i]['bond_overlap'])
        else:
            best_pos = max(range(len(ts['modes'])),
                           key=lambda i: ts['modes'][i]['bond_overlap'])
        ts['default_mode_idx'] = best_pos

    new_data = json.dumps(data)
    new_text = text[:match.start()] + f"const DATA = {new_data};\n" + text[match.end():]

    # Also patch the JS to display bond_overlap.
    new_text = new_text.replace(
        '<th>rxn ovlp</th><th>core frac</th>',
        '<th>bond ovlp</th><th>rxn ovlp</th><th>core frac</th>',
    )
    new_text = new_text.replace(
        '<td>${m.freq.toFixed(2)}</td><td>${m.rxn_overlap.toFixed(3)}</td><td>${m.core_fraction.toFixed(3)}</td>',
        '<td>${m.freq.toFixed(2)}</td><td>${m.bond_overlap.toFixed(3)}</td><td>${m.rxn_overlap.toFixed(3)}</td><td>${m.core_fraction.toFixed(3)}</td>',
    )
    new_text = new_text.replace(
        '`mode #${m.idx} · ${tag} freq=${f} cm⁻¹ · rxn_overlap=${m.rxn_overlap.toFixed(3)} · core_frac=${m.core_fraction.toFixed(3)}`',
        '`mode #${m.idx} · ${tag} freq=${f} cm⁻¹ · bond_ovlp=${m.bond_overlap.toFixed(3)} · rxn_ovlp=${m.rxn_overlap.toFixed(3)} · core_frac=${m.core_fraction.toFixed(3)}`',
    )

    html_path.write_text(new_text)
    return True


def main():
    files = sorted(OUT_DIR.glob("*.html"))
    files = [f for f in files if f.name not in ('index.html', 'flat_view.html')]
    print(f"Patching {len(files)} per-step HTMLs (no recompute)...")
    t0 = time.time()
    n_ok = n_fail = 0
    for i, hp in enumerate(files, 1):
        try:
            if patch_step(hp):
                n_ok += 1
            else:
                n_fail += 1
                print(f"  [{i}] {hp.name}: no DATA found")
        except Exception as e:
            n_fail += 1
            print(f"  [{i}] {hp.name}: {e}")
    print(f"\nDone in {time.time()-t0:.1f}s. {n_ok} ok, {n_fail} failed.")


if __name__ == "__main__":
    main()
