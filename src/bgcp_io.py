"""
Shared BGCP IO helpers: paths, charge/mult lookup, xyz concatenation.

Used by run_pq_bgcp.py, build_pq_viewer.py, build_pq_regression_viewers.py,
build_pq_regression_traces.py.
"""
from __future__ import annotations
import json
import re
from pathlib import Path

import numpy as np


BGCP_ROOT = Path("/Users/yunhengz/empty_for_claude/rxn_core/Benchmark_Guesses_Collective_Package")
TSDISCO_INDEX = Path(
    "/Users/yunhengz/empty_for_claude/tsdisco_benchmark_visualization_plain_portable/index.html"
)
ROOT = Path(__file__).parent
WORK = ROOT / "work_bgcp"
OUT = ROOT / "out"
OUT.mkdir(parents=True, exist_ok=True)
WORK.mkdir(parents=True, exist_ok=True)


def load_tsdisco_chg_uhf_lookup():
    """step_id -> (charge, uhf) by parsing tsdisco index.html.
    Falls back to empty dict if the index isn't accessible."""
    if not TSDISCO_INDEX.exists():
        return {}
    html = TSDISCO_INDEX.read_text()
    m = re.search(r"const TSDISCO_DATA = (\{.*?\});\s*\n", html, re.DOTALL)
    if m is None:
        return {}
    data = json.loads(m.group(1))
    out = {}
    for s in data["steps"]:
        chg = s.get("charge", 0) or 0
        mult = s.get("multiplicity", 1) or 1
        out[s["step_id"]] = (int(chg), max(0, int(mult) - 1))
    return out


LOOKUP = load_tsdisco_chg_uhf_lookup()


def list_step_dirs():
    return sorted(d for d in BGCP_ROOT.iterdir()
                  if d.is_dir() and not d.name.startswith("."))


def list_initial_guesses(step_dir):
    return sorted((step_dir / "initial_guess").glob("*.xyz"))


def iter_num(path):
    """Extract iter number from filename like '...iter12_xxxx.xyz'."""
    m = re.search(r"_iter(\d+)_", path.name)
    return int(m.group(1)) if m else 9999


def _parse_xyz_text(txt):
    lines = txt.strip().splitlines()
    n = int(lines[0].split()[0])
    elements, coords = [], []
    for ln in lines[2:2 + n]:
        parts = ln.split()
        elements.append(parts[0])
        coords.append([float(x) for x in parts[1:4]])
    return elements, np.array(coords)


def concat_xyz(xyz_texts, gap=10.0):
    """Concatenate multiple xyz texts. Multi-fragment inputs are translated
    apart along the x-axis with `gap` Å between bounding boxes, so xtb sees
    them as well-separated species (no spurious inter-fragment WBOs)."""
    if len(xyz_texts) == 1:
        return xyz_texts[0]
    parsed = [_parse_xyz_text(t) for t in xyz_texts]
    elements, coords_list, cursor = [], [], 0.0
    for els, coords in parsed:
        if coords.size == 0:
            continue
        bb_min = coords.min(axis=0)
        bb_max = coords.max(axis=0)
        shift = np.array([
            cursor - bb_min[0],
            -(bb_max[1] + bb_min[1]) / 2.0,
            -(bb_max[2] + bb_min[2]) / 2.0,
        ])
        coords_list.append(coords + shift)
        elements.extend(els)
        cursor += (bb_max[0] - bb_min[0]) + gap
    all_coords = np.vstack(coords_list)
    body = "\n".join(f"{el}  {c[0]:.6f}  {c[1]:.6f}  {c[2]:.6f}"
                     for el, c in zip(elements, all_coords))
    return f"{len(elements)}\nconcatenated_separated\n{body}\n"


def read_xyzs(dir_path):
    """Concatenated xyz text for all *.xyz under dir_path (or None)."""
    files = sorted(dir_path.glob("*.xyz"))
    if not files:
        return None
    return concat_xyz([f.read_text() for f in files])
