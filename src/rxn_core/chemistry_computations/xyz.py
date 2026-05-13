"""XYZ parsing and formatting utilities."""
from __future__ import annotations

from pathlib import Path

import numpy as np


def parse_xyz(path):
    """Parse an XYZ file into ``(elements, coords)``."""
    lines = Path(path).read_text().strip().splitlines()
    n = int(lines[0])
    elements, coords = [], []
    for ln in lines[2:2 + n]:
        parts = ln.split()
        elements.append(parts[0])
        coords.append([float(x) for x in parts[1:4]])
    return elements, np.array(coords)


def parse_xyz_file(path: Path):
    """Compatibility wrapper around :func:`parse_xyz`."""
    return parse_xyz(path)


def read_first_xyz(directory: Path):
    """Read the first ``*.xyz`` file in a directory."""
    files = sorted(Path(directory).glob("*.xyz"))
    if not files:
        raise FileNotFoundError(f"no xyz in {directory}")
    return parse_xyz(files[0])


def write_xyz_str(elements, coords, comment=""):
    """Format coordinates as a standard XYZ string."""
    out = [str(len(elements)), comment]
    for el, (x, y, z) in zip(elements, coords):
        out.append(f"{el}  {x:.6f}  {y:.6f}  {z:.6f}")
    return "\n".join(out) + "\n"


def xyz_block(elements, xyz, comment=""):
    """Format coordinates as an XYZ block with wider numeric fields."""
    lines = [str(len(elements)), comment]
    for el, c in zip(elements, xyz):
        lines.append(f"{el:<3s}  {c[0]:14.8f}  {c[1]:14.8f}  {c[2]:14.8f}")
    return "\n".join(lines) + "\n"


def xyz_with_disp(elements, xyz, disp, comment=""):
    """Format an extended XYZ block: ``element x y z dx dy dz``."""
    lines = [str(len(elements)), comment]
    for el, c, d in zip(elements, xyz, disp):
        lines.append(
            f"{el:<3s}  {c[0]:14.8f}  {c[1]:14.8f}  {c[2]:14.8f}"
            f"   {d[0]:10.6f}  {d[1]:10.6f}  {d[2]:10.6f}"
        )
    return "\n".join(lines) + "\n"
