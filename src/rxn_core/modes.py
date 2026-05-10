"""
Per-mode reactive features.

Given a normal-mode displacement matrix in R-frame indexing, computes:
  - bond_overlap (beta)  : projection on the broken/formed bond-axis vector V
  - rxn_overlap  (rho)   : projection on the reaction-coord vector Delta restricted to core
  - core_fraction (kappa): fraction of mode energy on core atoms (computed inline by callers)

All features are sign-blind, in [0, 1], and zero on spectator atoms by construction
(V_hat is sparse on non-core; Delta_hat_core is sparse on non-core).

Helpers also include the g98.out parser, Kabsch alignment, the reaction-coord
delta builder, the bond-reaction-vector builder, and mode reindexing into
R-frame.
"""
from __future__ import annotations
import re
from pathlib import Path

import numpy as np


def parse_g98_modes(path):
    """Parse vibration modes from xtb's g98.out (Gaussian 98 format).

    Returns:
      freqs:   (n_modes,) float array, cm^-1, negative = imaginary
      modes:   (n_modes, n_atoms, 3) float array of Cartesian displacements
    """
    text = Path(path).read_text()
    lines = text.splitlines()
    freqs_all = []
    blocks = []
    i = 0
    while i < len(lines):
        line = lines[i]
        if line.lstrip().startswith('Frequencies --'):
            tokens = line.split('--', 1)[1].split()
            block_freqs = [float(t) for t in tokens]
            j = i + 1
            while j < len(lines) and not lines[j].lstrip().startswith('Atom'):
                j += 1
            j += 1  # skip header
            atoms_disp = []  # rows of [d1x,d1y,d1z,d2x,d2y,d2z,d3x,d3y,d3z]
            while j < len(lines):
                ln = lines[j].strip()
                if not ln or ln.startswith('Frequencies'):
                    break
                parts = ln.split()
                if len(parts) < 2 + 3 * len(block_freqs):
                    break
                vals = [float(x) for x in parts[2:2 + 3 * len(block_freqs)]]
                atoms_disp.append(vals)
                j += 1
            arr = np.asarray(atoms_disp, dtype=float)
            n_atoms = arr.shape[0]
            for k in range(len(block_freqs)):
                disp = arr[:, 3 * k:3 * (k + 1)]  # (n_atoms, 3)
                freqs_all.append(block_freqs[k])
                blocks.append(disp)
            i = j
        else:
            i += 1
    if not blocks:
        return np.zeros(0), np.zeros((0, 0, 3))
    freqs = np.asarray(freqs_all, dtype=float)
    modes = np.stack(blocks, axis=0)  # (n_modes, n_atoms, 3)
    return freqs, modes


def kabsch(P, Q):
    """Rigid alignment of Q onto P. Returns (R, t) so that aligned = (Q - <Q>) @ R + <P>.
    Both P, Q are (n, 3) arrays."""
    P = np.asarray(P, dtype=float); Q = np.asarray(Q, dtype=float)
    Pc = P.mean(0); Qc = Q.mean(0)
    H = (Q - Qc).T @ (P - Pc)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    t = Pc - Qc @ R.T
    return R, t


def reaction_coord_delta(xyzR, xyzP, mapping_RP):
    """Per-atom reaction-coordinate displacement Delta_i = P[m(i)] - R[i]
    after Kabsch-aligning P (in R-atom order) to R. Mapped atoms only;
    unmapped get 0. Returns (n_R, 3) array."""
    n_R = len(xyzR)
    p_in_R_order = np.zeros((n_R, 3))
    mapped = []
    for r, p in mapping_RP.items():
        p_in_R_order[r] = xyzP[p]
        mapped.append(r)
    if not mapped:
        return np.zeros((n_R, 3))
    idx = np.array(mapped)
    R, t = kabsch(np.asarray(xyzR)[idx], p_in_R_order[idx])
    p_aligned = (p_in_R_order @ R.T) + t
    delta = p_aligned - np.asarray(xyzR)
    in_map = np.zeros(n_R, dtype=bool); in_map[idx] = True
    delta[~in_map] = 0
    return delta


def bond_reaction_vector(xyz_TS_in_R, broken_bonds_R, formed_bonds_R):
    """Per-atom reaction direction at TS coordinates.

    For each broken bond (i, j): atoms should be moving APART --
    V[i] -= unit, V[j] += unit  where unit = (xyz[j]-xyz[i])/||.||.
    For each formed bond (i, j): atoms should be moving TOGETHER --
    V[i] += unit, V[j] -= unit.

    Sign convention is chosen so that a true concerted reaction mode
    accumulates contributions coherently. Returns (n_atoms, 3) ndarray,
    NOT normalized."""
    xyz = np.asarray(xyz_TS_in_R, dtype=float)
    V = np.zeros_like(xyz)
    for i, j in broken_bonds_R:
        v = xyz[j] - xyz[i]
        n = float(np.linalg.norm(v))
        if n < 1e-9: continue
        u = v / n
        V[i] -= u
        V[j] += u
    for i, j in formed_bonds_R:
        v = xyz[j] - xyz[i]
        n = float(np.linalg.norm(v))
        if n < 1e-9: continue
        u = v / n
        V[i] += u
        V[j] -= u
    return V


def bond_overlap_per_mode(modes_R, V):
    """|mode . V_hat| / ||mode||  in [0, 1].

    V is the bond-reaction-direction vector from bond_reaction_vector().
    Numerator captures alignment with the bond-stretch/contract motion;
    denominator (full mode norm) penalizes amplitude on spectator atoms."""
    n_modes = modes_R.shape[0]
    if n_modes == 0:
        return np.zeros(0)
    v_flat = np.asarray(V).reshape(-1)
    v_norm = float(np.linalg.norm(v_flat))
    if v_norm < 1e-9:
        return np.zeros(n_modes)
    v_unit = v_flat / v_norm
    m_flat = modes_R.reshape(n_modes, -1)
    m_norm = np.linalg.norm(m_flat, axis=1)
    dots = np.abs(m_flat @ v_unit)
    out = np.zeros(n_modes)
    valid = m_norm > 1e-9
    out[valid] = dots[valid] / m_norm[valid]
    return out


def rxn_overlap_per_mode(modes_R, delta, core_atoms):
    """Reaction-mode overlap: projection of each full-molecule mode onto
    the unit reaction-coord vector restricted to core atoms.

    Builds Delta_hat_core (zero outside core, normalized over core),
    then q_m = |d_m . Delta_hat_core| / ||d_m|| in [0, 1]. Penalizes
    spectator-atom motion (in denominator) and core motion not aligned
    with Delta (in numerator)."""
    n_modes = modes_R.shape[0]
    if n_modes == 0 or not core_atoms:
        return np.zeros(n_modes)
    core = np.asarray(core_atoms, dtype=int)
    delta_pad = np.zeros_like(modes_R[0])
    delta_pad[core] = delta[core]
    d_flat = delta_pad.reshape(-1)
    d_norm = float(np.linalg.norm(d_flat))
    if d_norm < 1e-9:
        return np.zeros(n_modes)
    d_unit = d_flat / d_norm
    m_flat = modes_R.reshape(n_modes, -1)
    m_norm = np.linalg.norm(m_flat, axis=1)
    dots = np.abs(m_flat @ d_unit)
    out = np.zeros(n_modes)
    valid = m_norm > 1e-9
    out[valid] = dots[valid] / m_norm[valid]
    return out


def core_atoms_in_R_frame(mapping_R_to_P, broken, formed):
    """Atoms touching broken (R indices) or formed (P, mapped back via inv).
    Returns sorted list of R-frame indices."""
    inv = {v: k for k, v in mapping_R_to_P.items()}
    core = set()
    for (i, j, _, _) in broken:
        core.add(int(i)); core.add(int(j))
    for (i, j, _, _) in formed:
        if i in inv: core.add(int(inv[i]))
        if j in inv: core.add(int(inv[j]))
    return sorted(core)


def reindex_modes_to_R(modes_TS, mapping_R_to_TS, n_R):
    """Reindex mode displacements from TS-frame to R-frame indexing,
    so modes_R[k][i_R] == modes_TS[k][mapping_R_to_TS[i_R]]. Returns
    (n_modes, n_R, 3)."""
    n_modes = modes_TS.shape[0]
    out = np.zeros((n_modes, n_R, 3))
    for r, t in mapping_R_to_TS.items():
        out[:, r] = modes_TS[:, t]
    return out
