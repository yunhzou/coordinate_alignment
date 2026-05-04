"""
WBO-anchor atom mapping & reaction-core detection.

From-first-principles algorithm:
  1. Run xtb single-point on reactant and product to get Wiberg Bond Orders.
  2. Compute Morgan-like atom signatures, but the recursive shell uses
     (rounded WBO, neighbor_signature) pairs instead of bond-type labels.
     The signature is therefore continuous and naturally symmetry-breaking.
  3. For each atom, find the largest radius at which its signature is
     unique within its own complex AND uniquely matches a partner atom in
     the other complex. That match becomes an anchor.
  4. Propagate the mapping outward from anchors along bonds. Greedy match
     at each step uses (signature-radius-match, |delta-WBO|).
  5. Atoms whose signature never matches at any radius -> reaction core.
  6. Bonds present in R but not in P (via mapping) -> broken; vice versa -> formed.
"""

from __future__ import annotations

import shutil
import subprocess
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from scipy.optimize import linear_sum_assignment


# -------------------- IO --------------------

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


def run_xtb(xyz_path, workdir, charge=0, uhf=0):
    """Run xtb GFN2 single-point. Returns (elements, coords, wbo NxN)."""
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    local = workdir / Path(xyz_path).name
    shutil.copy(xyz_path, local)
    cmd = ["xtb", local.name, "--gfn", "2", "--sp"]
    if charge:
        cmd += ["--chrg", str(charge)]
    if uhf:
        cmd += ["--uhf", str(uhf)]
    res = subprocess.run(cmd, cwd=workdir, capture_output=True, text=True)
    if res.returncode != 0:
        raise RuntimeError(f"xtb failed: {res.stderr[-500:]}")
    elements, coords = parse_xyz(local)
    n = len(elements)
    wbo = np.zeros((n, n))
    wbo_file = workdir / "wbo"
    if not wbo_file.exists():
        raise RuntimeError("xtb did not produce a wbo file")
    for ln in wbo_file.read_text().splitlines():
        parts = ln.split()
        if len(parts) < 3:
            continue
        i, j = int(parts[0]) - 1, int(parts[1]) - 1
        v = float(parts[2])
        wbo[i, j] = v
        wbo[j, i] = v
    return elements, coords, wbo


# -------------------- signatures --------------------

def neighbor_list(wbo, eps=0.05):
    """eps is a numerical-noise floor only; not a chemistry threshold."""
    n = wbo.shape[0]
    return [
        [(j, float(wbo[i, j])) for j in range(n) if i != j and wbo[i, j] > eps]
        for i in range(n)
    ]


def signatures(elements, wbo, max_radius=4, wbo_round=1, eps=0.05):
    """
    Per-atom Morgan-like hashes at radii 0..max_radius.

    radius 0 = element only.
    radius r = (sig_at_r-1[i],  sorted_tuple( (round(WBO[i,j],wbo_round), sig_at_r-1[j])  for j in nb[i] ))
    """
    n = len(elements)
    nb = neighbor_list(wbo, eps=eps)
    sigs = {0: {i: elements[i] for i in range(n)}}
    for r in range(1, max_radius + 1):
        prev = sigs[r - 1]
        cur = {}
        for i in range(n):
            shell = tuple(sorted(
                (round(w, wbo_round), prev[j]) for j, w in nb[i]
            ))
            cur[i] = (prev[i], shell)
        sigs[r] = cur
    return sigs, nb


# -------------------- anchor finding --------------------

def find_anchors(sigs_R, sigs_P, elements_R, elements_P, max_radius,
                 min_anchor_radius=2):
    """
    Per-element global Hungarian on signature-match-radius. Commit a pair
    (i_R, j_P) as an anchor only if their signatures agree at radius
    >= min_anchor_radius -- this gates against atoms in the reaction core,
    whose local environment has changed.

    For symmetry-equivalent atoms (e.g. ortho carbons of a phenyl ring),
    Hungarian arbitrarily picks one of the equivalent pairings, which is
    fine because the atoms are chemically interchangeable. Propagation
    then extends consistently from there.
    """
    by_el_R = defaultdict(list)
    by_el_P = defaultdict(list)
    for i, e in enumerate(elements_R):
        by_el_R[e].append(i)
    for j, e in enumerate(elements_P):
        by_el_P[e].append(j)
    anchors = []
    BIG = 1e6
    for el, rs in by_el_R.items():
        ps = by_el_P.get(el, [])
        if not ps:
            continue
        n = max(len(rs), len(ps))
        match_r = np.full((n, n), -1)
        cost = np.full((n, n), BIG)
        for a, i in enumerate(rs):
            for b, j in enumerate(ps):
                m = -1
                for r in range(max_radius, -1, -1):
                    if sigs_R[r][i] == sigs_P[r][j]:
                        m = r
                        break
                match_r[a, b] = m
                cost[a, b] = -m
        row, col = linear_sum_assignment(cost)
        for a, b in zip(row, col):
            if a < len(rs) and b < len(ps):
                mr = int(match_r[a, b])
                if mr >= min_anchor_radius:
                    anchors.append((rs[a], ps[b], mr))
    return anchors


# -------------------- propagation --------------------

def propagate(anchors, sigs_R, sigs_P, nb_R, nb_P,
              elements_R, elements_P, max_radius,
              min_match_radius=1):
    """
    BFS outward from anchors. For each unmapped R-neighbor of a mapped atom,
    pick the best partner among unmapped P-neighbors of the corresponding
    mapped P atom. Score = (largest-radius signature match, -|delta WBO|).

    Commit if signatures agree at radius >= min_match_radius. Default 0 means
    element-only is enough to commit -- the signature radius only ranks among
    candidates. This is correct: when an atom's local environment changes
    because of the reaction, we still want to map the atom (otherwise we
    cannot compute the WBO delta), but the ranking still pulls toward the
    chemically-closer choice when several are available.
    """
    mapping = {i: j for (i, j, _) in anchors}
    inv = {j: i for i, j in mapping.items()}
    frontier = [i for (i, _, _) in anchors]

    def pair_score(k, wR, l, wP):
        """Larger is better. (match_r, -|delta wbo|) lifted to a scalar."""
        match_r = -1
        for r in range(max_radius, -1, -1):
            if sigs_R[r][k] == sigs_P[r][l]:
                match_r = r
                break
        return match_r * 100.0 - abs(wR - wP)

    while frontier:
        next_frontier = []
        for i in frontier:
            j = mapping[i]
            unmapped_r = [(k, w) for (k, w) in nb_R[i] if k not in mapping]
            unmapped_p = [(l, w) for (l, w) in nb_P[j] if l not in inv]
            # Group by element on both sides so we only match like to like.
            by_el_R = defaultdict(list)
            by_el_P = defaultdict(list)
            for (k, w) in unmapped_r:
                by_el_R[elements_R[k]].append((k, w))
            for (l, w) in unmapped_p:
                by_el_P[elements_P[l]].append((l, w))
            # For each element, do Hungarian assignment between R-side and
            # P-side candidates.  This avoids the cascading errors of
            # greedy local picks (e.g. mapping a carbonyl-O to an
            # iodide-coordinating-O just because it was processed first).
            for el, rs in by_el_R.items():
                ps = by_el_P.get(el, [])
                if not ps:
                    continue
                M, N = len(rs), len(ps)
                cost = np.zeros((M, N))
                for a, (k, wR) in enumerate(rs):
                    for b, (l, wP) in enumerate(ps):
                        cost[a, b] = -pair_score(k, wR, l, wP)
                row_ind, col_ind = linear_sum_assignment(cost)
                for a, b in zip(row_ind, col_ind):
                    k, wR = rs[a]
                    l, wP = ps[b]
                    # recover the (match_r, -|deltaW|) for the gating check
                    s = -cost[a, b]
                    match_r = int(round((s + abs(wR - wP)) / 100.0))
                    if match_r >= min_match_radius:
                        mapping[k] = l
                        inv[l] = k
                        next_frontier.append(k)
        frontier = next_frontier
    return mapping


# -------------------- core phase: 3D-based residual mapping --------------------

def kabsch(P, Q):
    """Return rotation R and translation t that minimize ||R*P + t - Q||."""
    Pc, Qc = P.mean(0), Q.mean(0)
    H = (P - Pc).T @ (Q - Qc)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1, 1, d])
    R = Vt.T @ D @ U.T
    t = Qc - R @ Pc
    return R, t


def connected_components(nb, n):
    seen = [False] * n
    comps = []
    for s in range(n):
        if seen[s]:
            continue
        comp = []
        stack = [s]
        while stack:
            x = stack.pop()
            if seen[x]:
                continue
            seen[x] = True
            comp.append(x)
            for (j, _) in nb[x]:
                if not seen[j]:
                    stack.append(j)
        comps.append(sorted(comp))
    return comps


def map_core_by_geometry(mapping, coords_R, coords_P, elements_R, elements_P, nb_R,
                         k_local=8):
    """
    For atoms unmapped after the spectator phase: project each one into
    P-space using a LOCAL Kabsch built from its k nearest already-mapped
    neighbors in R. Then per-element Hungarian on 3D distances.

    Local-instead-of-global Kabsch matters when the spectator atoms cluster
    in one region of the molecule (e.g. a substrate tail far from the
    reactive site) -- a single global rotation overfits that cluster and
    mis-projects atoms in the geometrically-distorted reactive region. A
    local frame built from nearby spectators is robust to such distortion.

    Iterates: each pass commits some atoms, which then become local-frame
    seeds for further passes. This sweeps inward from the spectator region
    toward the core.
    """
    if not mapping:
        return mapping
    inv = {j: i for i, j in mapping.items()}

    for _outer_pass in range(20):
        unmapped_R = [i for i in range(len(elements_R)) if i not in mapping]
        unmapped_P = [j for j in range(len(elements_P)) if j not in inv]
        if not unmapped_R or not unmapped_P:
            break

        mapped_R_idx = np.array(sorted(mapping.keys()))
        if len(mapped_R_idx) < 3:
            break

        by_el_R = defaultdict(list)
        by_el_P = defaultdict(list)
        for i in unmapped_R:
            by_el_R[elements_R[i]].append(i)
        for j in unmapped_P:
            by_el_P[elements_P[j]].append(j)

        progressed = False
        for el, rs in by_el_R.items():
            ps = by_el_P.get(el, [])
            if not ps:
                continue
            # Local Kabsch projection per R atom
            projections = []
            confidence = []  # smaller = better
            for i in rs:
                d = np.linalg.norm(coords_R[mapped_R_idx] - coords_R[i], axis=1)
                k = min(k_local, len(mapped_R_idx))
                order = np.argsort(d)[:k]
                nn_R = mapped_R_idx[order]
                nn_P = np.array([mapping[int(x)] for x in nn_R])
                Rmat, tvec = kabsch(coords_R[nn_R], coords_P[nn_P])
                projections.append(Rmat @ coords_R[i] + tvec)
                # confidence = max distance to local spectator (smaller = denser local frame)
                confidence.append(d[order[-1]])
            cost = np.zeros((len(rs), len(ps)))
            for a, proj in enumerate(projections):
                for b, j in enumerate(ps):
                    cost[a, b] = np.linalg.norm(proj - coords_P[j])
            row_ind, col_ind = linear_sum_assignment(cost)
            # Sort commits by confidence: most confident (smallest local frame) first
            commit = sorted(zip(row_ind, col_ind, confidence), key=lambda x: x[2])
            # Commit half per pass (so confident matches inform later passes)
            n_commit = max(1, len(commit) // 2) if _outer_pass < 5 else len(commit)
            for a, b, _ in commit[:n_commit]:
                mapping[rs[a]] = ps[b]
                inv[ps[b]] = rs[a]
                progressed = True
        if not progressed:
            break
    return mapping


# -------------------- cleanup pass --------------------

def cleanup_pass(mapping, nb_R, nb_P, elements_R, elements_P, max_iter=10):
    """
    For atoms still unmapped after BFS propagation: if at least one of their
    bonded neighbors *is* mapped, look at the mapped partner's unmapped
    neighbors of the same element. If a unique candidate exists, commit it.
    Repeats until no further progress.
    """
    inv = {j: i for i, j in mapping.items()}
    nR, nP = len(elements_R), len(elements_P)
    for _ in range(max_iter):
        progressed = False
        for k in range(nR):
            if k in mapping:
                continue
            mapped_neighbors = [(j, w) for (j, w) in nb_R[k] if j in mapping]
            if not mapped_neighbors:
                continue
            cand_set = None
            for (j, _) in mapped_neighbors:
                jP = mapping[j]
                opts = {l for (l, _) in nb_P[jP]
                        if l not in inv and elements_P[l] == elements_R[k]}
                cand_set = opts if cand_set is None else (cand_set & opts)
                if not cand_set:
                    break
            if cand_set and len(cand_set) == 1:
                l = cand_set.pop()
                mapping[k] = l
                inv[l] = k
                progressed = True
        if not progressed:
            break
    return mapping


# -------------------- bond classification --------------------

def classify_bonds(mapping, wbo_R, wbo_P, bond_cut=0.5):
    """
    Returns:
      broken : (i_R, j_R, wbo_R_value, wbo_P_value_or_None)
      formed : (i_P, j_P, wbo_R_value_or_None, wbo_P_value)
      core_R, core_P : sets of atoms involved in change or unmapped
    """
    inv = {j: i for i, j in mapping.items()}
    nR, nP = wbo_R.shape[0], wbo_P.shape[0]
    broken, formed = [], []

    # R-side bonds
    for i in range(nR):
        for j in range(i + 1, nR):
            wR = wbo_R[i, j]
            if wR < bond_cut:
                continue
            if i not in mapping or j not in mapping:
                broken.append((i, j, float(wR), None))
                continue
            wP = wbo_P[mapping[i], mapping[j]]
            if wP < bond_cut:
                broken.append((i, j, float(wR), float(wP)))

    # P-side bonds
    for iP in range(nP):
        for jP in range(iP + 1, nP):
            wP = wbo_P[iP, jP]
            if wP < bond_cut:
                continue
            if iP not in inv or jP not in inv:
                formed.append((iP, jP, None, float(wP)))
                continue
            wR = wbo_R[inv[iP], inv[jP]]
            if wR < bond_cut:
                formed.append((iP, jP, float(wR), float(wP)))

    core_R = set(i for i in range(nR) if i not in mapping)
    core_P = set(j for j in range(nP) if j not in inv)
    for (i, j, _, _) in broken:
        core_R.add(i); core_R.add(j)
    for (iP, jP, _, _) in formed:
        core_P.add(iP); core_P.add(jP)

    return broken, formed, sorted(core_R), sorted(core_P)


# -------------------- swap refinement --------------------

def refine_mapping_by_swaps(mapping, wbo_R, wbo_P, elements_R, elements_P,
                            bond_cut=0.5, max_iters=20):
    """
    Hill-climb: try pairwise swaps of same-element atoms that reduce the total
    number of broken+formed bonds. Local Kabsch mapping can permute spectator
    pairs whose 3D positions are close (e.g. an ester carbonyl-C and its
    adjacent O-CH3 carbon); those swaps show up as paired broken/formed
    bonds whose fingerprints don't match (so they're not caught by the
    symmetry-reconciliation step). A swap of the two atoms restores
    consistency in O(local-bond-count) operations per attempt.
    """
    nR, nP = wbo_R.shape[0], wbo_P.shape[0]
    inv = {v: k for k, v in mapping.items()}

    def total_changes(m):
        invm = {v: k for k, v in m.items()}
        b = f = 0
        for i in range(nR):
            for j in range(i + 1, nR):
                wR = wbo_R[i, j]
                if wR < bond_cut:
                    continue
                if i not in m or j not in m:
                    b += 1
                    continue
                if wbo_P[m[i], m[j]] < bond_cut:
                    b += 1
        for ip in range(nP):
            for jp in range(ip + 1, nP):
                wP = wbo_P[ip, jp]
                if wP < bond_cut:
                    continue
                if ip not in invm or jp not in invm:
                    f += 1
                    continue
                if wbo_R[invm[ip], invm[jp]] < bond_cut:
                    f += 1
        return b + f

    def problematic(m):
        invm = {v: k for k, v in m.items()}
        atoms = set()
        for i in range(nR):
            for j in range(i + 1, nR):
                if wbo_R[i, j] >= bond_cut:
                    if i not in m or j not in m or wbo_P[m[i], m[j]] < bond_cut:
                        atoms.add(i); atoms.add(j)
        for ip in range(nP):
            for jp in range(ip + 1, nP):
                if wbo_P[ip, jp] >= bond_cut:
                    if ip in invm and jp in invm:
                        if wbo_R[invm[ip], invm[jp]] < bond_cut:
                            atoms.add(invm[ip]); atoms.add(invm[jp])
        return atoms

    cur = total_changes(mapping)
    for _ in range(max_iters):
        bad = problematic(mapping)
        if not bad:
            break
        moved = False
        for i in list(bad):
            if i not in mapping:
                continue
            best_j, best_score = None, cur
            for j in list(mapping.keys()):
                if i == j or elements_R[i] != elements_R[j]:
                    continue
                # tentative swap
                mapping[i], mapping[j] = mapping[j], mapping[i]
                s = total_changes(mapping)
                # revert
                mapping[i], mapping[j] = mapping[j], mapping[i]
                if s < best_score:
                    best_score = s
                    best_j = j
            if best_j is not None:
                mapping[i], mapping[best_j] = mapping[best_j], mapping[i]
                cur = best_score
                moved = True
                break
        if not moved:
            break
    return mapping


# -------------------- symmetry reconciliation --------------------

def reconcile_symmetric(broken, formed, sigs_R, sigs_P, fp_radius=1):
    """
    Cancel broken+formed bond pairs whose chemistry fingerprints match.

    A "fingerprint" is (frozen multiset of endpoint signatures at fp_radius,
    rounded WBO). When a bond appears as broken in R and an identical
    fingerprint appears as formed in P, the most likely explanation is that
    the atom mapping permuted equivalent atoms inside a symmetric region
    (e.g. the two ortho carbons of a phenyl ring) -- the bond is conserved.
    """
    def fp(i, j, w, sigs):
        return (frozenset((sigs[fp_radius][i], sigs[fp_radius][j])), round(w, 1))

    formed_by_fp = defaultdict(list)
    for idx, (iP, jP, _, wP) in enumerate(formed):
        formed_by_fp[fp(iP, jP, wP, sigs_P)].append(idx)

    cancel_broken, cancel_formed = set(), set()
    for idx, (i, j, wR, _) in enumerate(broken):
        bucket = formed_by_fp.get(fp(i, j, wR, sigs_R), [])
        if bucket:
            cancel_broken.add(idx)
            cancel_formed.add(bucket.pop())

    new_broken = [b for i, b in enumerate(broken) if i not in cancel_broken]
    new_formed = [f for i, f in enumerate(formed) if i not in cancel_formed]
    return new_broken, new_formed


# -------------------- top-level driver --------------------

def analyze(reactant_xyz, product_xyz, workdir,
            charge=0, uhf=0, max_radius=4, bond_cut=0.5):
    workdir = Path(workdir)
    elR, xyzR, wboR = run_xtb(reactant_xyz, workdir / "R", charge=charge, uhf=uhf)
    elP, xyzP, wboP = run_xtb(product_xyz, workdir / "P", charge=charge, uhf=uhf)
    if len(elR) != len(elP):
        raise ValueError(f"Atom-count mismatch: R={len(elR)} P={len(elP)}")
    if Counter(elR) != Counter(elP):
        raise ValueError(f"Element composition mismatch: {Counter(elR)} vs {Counter(elP)}")

    sigs_R, nb_R = signatures(elR, wboR, max_radius=max_radius)
    sigs_P, nb_P = signatures(elP, wboP, max_radius=max_radius)

    anchors = find_anchors(sigs_R, sigs_P, elR, elP, max_radius,
                           min_anchor_radius=2)
    # phase 1: BFS propagate from anchors. Hungarian per shell uses
    # (signature-match-radius, 3D distance after running Kabsch) -- the 3D
    # distance term breaks symmetry between same-element neighbors whose
    # local environment is similar.
    mapping = propagate(
        anchors, sigs_R, sigs_P, nb_R, nb_P,
        elR, elP, max_radius=max_radius, min_match_radius=0,
        coords_R=xyzR, coords_P=xyzP,
    )
    mapping = cleanup_pass(mapping, nb_R, nb_P, elR, elP)
    n_spectator = len(mapping)
    # phase 2: 3D-distance Hungarian via local Kabsch for residual atoms
    mapping = map_core_by_geometry(mapping, xyzR, xyzP, elR, elP, nb_R)
    # phase 3: pairwise-swap hill climb (catches simple permutation errors)
    mapping = refine_mapping_by_swaps(mapping, wboR, wboP, elR, elP, bond_cut=bond_cut)
    broken, formed, core_R, core_P = classify_bonds(
        mapping, wboR, wboP, bond_cut=bond_cut
    )
    broken, formed = reconcile_symmetric(broken, formed, sigs_R, sigs_P)
    # rebuild core atoms after reconciliation
    inv = {j: i for i, j in mapping.items()}
    core_R = set(i for i in range(len(elR)) if i not in mapping)
    core_P = set(j for j in range(len(elP)) if j not in inv)
    for (i, j, _, _) in broken:
        core_R.add(i); core_R.add(j)
    for (iP, jP, _, _) in formed:
        core_P.add(iP); core_P.add(jP)
    core_R, core_P = sorted(core_R), sorted(core_P)

    return dict(
        elements_R=elR, coords_R=xyzR, wbo_R=wboR,
        elements_P=elP, coords_P=xyzP, wbo_P=wboP,
        anchors=anchors, mapping=mapping,
        broken=broken, formed=formed,
        core_R=core_R, core_P=core_P,
        n_spectator=n_spectator,
    )
