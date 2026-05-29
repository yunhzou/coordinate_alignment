"""Cached xtb execution and WBO loading."""
from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

import numpy as np

from .xyz import parse_xyz


def read_wbo_file(path, n_atoms):
    """Read an xtb ``wbo`` file into a symmetric WBO matrix."""
    wbo = np.zeros((n_atoms, n_atoms))
    for ln in Path(path).read_text().splitlines():
        parts = ln.split()
        if len(parts) < 3:
            continue
        i, j = int(parts[0]) - 1, int(parts[1]) - 1
        value = float(parts[2])
        wbo[i, j] = value
        wbo[j, i] = value
    return wbo


def load_cached_xtb(workdir):
    """Load ``(elements, coords, wbo, xyz_path)`` from an xtb cache dir."""
    workdir = Path(workdir)
    xyz_files = [f for f in workdir.iterdir() if f.suffix == ".xyz"]
    if len(xyz_files) != 1:
        raise RuntimeError(f"expected 1 xyz in {workdir}, found {len(xyz_files)}")
    elements, coords = parse_xyz(xyz_files[0])
    wf = workdir / "wbo"
    if not wf.exists():
        raise RuntimeError(f"no wbo file in {workdir}")
    return elements, coords, read_wbo_file(wf, len(elements)), xyz_files[0]


def _xtb_command(local_name, mode, charge=0, uhf=0):
    cmd = ["xtb", local_name, "--gfn", "2", mode]
    if charge:
        cmd += ["--chrg", str(charge)]
    if uhf:
        cmd += ["--uhf", str(uhf)]
    return cmd


def run_xtb(xyz_path, workdir, charge=0, uhf=0, omp_threads=1):
    """Run cached xtb GFN2 single-point and return ``(elements, coords, wbo)``."""
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    local = workdir / Path(xyz_path).name
    src_text = Path(xyz_path).read_text()
    cached_text = local.read_text() if local.exists() else None
    wf = workdir / "wbo"
    cached = (cached_text == src_text) and wf.exists()
    if not cached:
        if Path(xyz_path).resolve() != local.resolve():
            shutil.copy(xyz_path, local)
        env = os.environ.copy()
        env["OMP_NUM_THREADS"] = str(max(1, int(omp_threads)))
        res = subprocess.run(
            _xtb_command(local.name, "--sp", charge=charge, uhf=uhf),
            cwd=workdir,
            capture_output=True,
            text=True,
            env=env,
        )
        (workdir / "xtb.stdout").write_text(res.stdout or "")
        (workdir / "xtb.stderr").write_text(res.stderr or "")
        if res.returncode != 0:
            raise RuntimeError(f"xtb failed: {res.stderr[-500:]}")
        if not wf.exists():
            raise RuntimeError("no wbo file")
    elements, coords = parse_xyz(local)
    return elements, coords, read_wbo_file(wf, len(elements))


def run_xtb_hess(xyz_path: Path, workdir: Path, charge: int = 0, uhf: int = 0,
                 omp_threads: int = 1):
    """Run cached ``xtb --hess`` and return geometry, WBO, frequencies, modes."""
    from ..modes import parse_g98_modes

    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    local = workdir / Path(xyz_path).name
    src_text = Path(xyz_path).read_text()
    cached_text = local.read_text() if local.exists() else None
    g98 = workdir / "g98.out"
    wbo = workdir / "wbo"
    cached = (cached_text == src_text) and g98.exists() and wbo.exists()
    if not cached:
        if Path(xyz_path).resolve() != local.resolve():
            shutil.copy(xyz_path, local)
        env = os.environ.copy()
        env["OMP_NUM_THREADS"] = str(max(1, int(omp_threads)))
        res = subprocess.run(
            _xtb_command(local.name, "--hess", charge=charge, uhf=uhf),
            cwd=workdir,
            capture_output=True,
            text=True,
            env=env,
        )
        (workdir / "xtbhess.stdout").write_text(res.stdout or "")
        (workdir / "xtbhess.stderr").write_text(res.stderr or "")
        if res.returncode != 0:
            raise RuntimeError(
                f"xtb hess failed in {workdir}: {res.stderr[-500:]}"
            )
        if not g98.exists() or not wbo.exists():
            raise RuntimeError(f"missing g98.out or wbo in {workdir}")
    elements, coords = parse_xyz(local)
    freqs, modes_ts = parse_g98_modes(g98)
    return elements, coords, read_wbo_file(wbo, len(elements)), freqs, modes_ts
