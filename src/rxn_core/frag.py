"""
Low-level utilities used by rxn_core_pq:

  parse_xyz / write_xyz_str — xyz IO
  run_xtb                   — GFN2 single-point with WBO matrix (cached)
  build_graph               — WBO graph (edge iff WBO >= bond_cut)
  expand_mapping            — element-multiset pairing of unmapped neighbors
  classify_bonds            — broken/formed bond classification by ΔWBO

The previous fragment-anchor / find_islands algorithm has been removed;
rxn_core_pq is the live alignment implementation.
"""
from __future__ import annotations

import shutil
import subprocess
from collections import defaultdict
from pathlib import Path

import numpy as np
import networkx as nx


# -------------------- xyz IO --------------------

def parse_xyz(path):
    lines = Path(path).read_text().strip().splitlines()
    n = int(lines[0])
    elements, coords = [], []
    for ln in lines[2:2 + n]:
        parts = ln.split()
        elements.append(parts[0])
        coords.append([float(x) for x in parts[1:4]])
    return elements, np.array(coords)


def write_xyz_str(elements, coords, comment=""):
    out = [str(len(elements)), comment]
    for el, (x, y, z) in zip(elements, coords):
        out.append(f"{el}  {x:.6f}  {y:.6f}  {z:.6f}")
    return "\n".join(out) + "\n"


# -------------------- xtb GFN2 single-point --------------------

def run_xtb(xyz_path, workdir, charge=0, uhf=0):
    """Run an xtb GFN2 single-point and return (elements, coords, wbo).
    Cached: re-uses an existing `wbo` file in workdir when the input xyz
    matches what's already there. xtb is the dominant cost so cache hits
    make repeated runs near-instant."""
    workdir = Path(workdir); workdir.mkdir(parents=True, exist_ok=True)
    local = workdir / Path(xyz_path).name
    src_text = Path(xyz_path).read_text()
    cached_text = local.read_text() if local.exists() else None
    wf = workdir / "wbo"
    cached = (cached_text == src_text) and wf.exists()
    if not cached:
        shutil.copy(xyz_path, local)
        cmd = ["xtb", local.name, "--gfn", "2", "--sp"]
        if charge: cmd += ["--chrg", str(charge)]
        if uhf:    cmd += ["--uhf",  str(uhf)]
        res = subprocess.run(cmd, cwd=workdir, capture_output=True, text=True)
        if res.returncode != 0:
            raise RuntimeError(f"xtb failed: {res.stderr[-500:]}")
        if not wf.exists():
            raise RuntimeError("no wbo file")
    elements, coords = parse_xyz(local)
    n = len(elements)
    wbo = np.zeros((n, n))
    for ln in wf.read_text().splitlines():
        parts = ln.split()
        if len(parts) < 3: continue
        i, j = int(parts[0]) - 1, int(parts[1]) - 1
        v = float(parts[2])
        wbo[i, j] = v; wbo[j, i] = v
    return elements, coords, wbo


# -------------------- WBO graph --------------------

def build_graph(elements, wbo, bond_cut=0.5):
    """Connectivity graph with element on each node and WBO weight on each
    edge. Bond exists iff WBO >= bond_cut. The PQ algorithm uses
    bond_cut = graph_floor = 0.2 to admit weak/partial bonds."""
    g = nx.Graph()
    for i, e in enumerate(elements):
        g.add_node(i, element=e)
    n = len(elements)
    for i in range(n):
        for j in range(i + 1, n):
            if wbo[i, j] >= bond_cut:
                g.add_edge(i, j, wbo=float(wbo[i, j]))
    return g


# -------------------- mapping post-processing --------------------

def expand_mapping(mapping, g_R, g_P):
    """Pure-connectivity expansion from already-mapped atoms. For each
    mapped pair (u, v), pair u's unmapped R-neighbors with v's unmapped
    P-neighbors element-by-element. If counts match exactly, commit the
    pairing in arbitrary order (symmetric atoms like methyl Hs are
    interchangeable). If counts differ, leave them unmapped — that's a
    real connectivity change at this atom (i.e. reaction-core).
    Loops until no further progress."""
    inv = {v: k for k, v in mapping.items()}
    while True:
        progressed = False
        for u in list(mapping.keys()):
            v = mapping[u]
            r_groups = defaultdict(list)
            for w in g_R.neighbors(u):
                if w in mapping: continue
                r_groups[g_R.nodes[w]['element']].append(w)
            p_groups = defaultdict(list)
            for x in g_P.neighbors(v):
                if x in inv: continue
                p_groups[g_P.nodes[x]['element']].append(x)
            for el, rs in r_groups.items():
                ps = p_groups.get(el, [])
                if len(ps) != len(rs):
                    continue
                for w, x in zip(rs, ps):
                    mapping[w] = x
                    inv[x] = w
                    progressed = True
        if not progressed:
            break
    return mapping


def classify_bonds(mapping, wbo_R, wbo_P, dwbo_threshold=0.5):
    """Bond classification by WBO change.

        broken iff (WBO_R - WBO_P) >= dwbo_threshold
        formed iff (WBO_P - WBO_R) >= dwbo_threshold

    Since WBO is non-negative, requiring (wR - wP) >= dwbo_threshold
    automatically implies wR >= dwbo_threshold (wP ≥ 0), so a single
    threshold both gates "is this a bond worth considering" and "did
    its order change enough." For pairs with one or both endpoints
    unmapped (no image bond defined) we treat the missing wP as 0 and
    apply the same threshold to wR.

    Returns (broken_list, formed_list, core_R, core_P) where each bond
    record is (i, j, wbo_R_or_None, wbo_P_or_None)."""
    inv = {v: k for k, v in mapping.items()}
    nR, nP = wbo_R.shape[0], wbo_P.shape[0]
    broken, formed = [], []
    for i in range(nR):
        for j in range(i + 1, nR):
            wR = wbo_R[i, j]
            if wR < dwbo_threshold: continue
            if i not in mapping or j not in mapping:
                broken.append((i, j, float(wR), None)); continue
            wP = wbo_P[mapping[i], mapping[j]]
            if wR - wP >= dwbo_threshold:
                broken.append((i, j, float(wR), float(wP)))
    for ip in range(nP):
        for jp in range(ip + 1, nP):
            wP = wbo_P[ip, jp]
            if wP < dwbo_threshold: continue
            if ip not in inv or jp not in inv:
                formed.append((ip, jp, None, float(wP))); continue
            wR = wbo_R[inv[ip], inv[jp]]
            if wP - wR >= dwbo_threshold:
                formed.append((ip, jp, float(wR), float(wP)))
    core_R = sorted({i for (i, j, _, _) in broken}
                    | {j for (i, j, _, _) in broken})
    core_P = sorted({i for (i, j, _, _) in formed}
                    | {j for (i, j, _, _) in formed})
    return broken, formed, core_R, core_P
