"""
Fragment-based reaction-core finder.

Algorithm (from first principles, WBO-graph fragment matching):

  1. Build a weighted graph from each WBO matrix: nodes = atoms (with element
     label), edges = atom pairs with WBO above a small noise floor (edge
     weight = continuous WBO).

  2. For each seed atom u in R and v in P with the same element, extract the
     radius-r connected subgraph rooted at u (and v). Search for the largest
     r at which the two subgraphs are graph-isomorphic, requiring:
        * node element equality
        * edge WBO equality within tolerance
        * the root atoms are mapped to each other (u -> v)
     The isomorphism provides an atom-by-atom correspondence inside the
     fragment -- ALL atoms inside the matched fragment become mounted, not
     just the root.

  3. Merge fragment matches: when two matches overlap and agree on the
     atoms in their intersection, take their union as a confident anchor.
     If they disagree, the boundary marks the reaction core.

  4. Expand: from the merged anchor region, walk outward through bonds. An
     unmapped neighbor is added if its local environment (the small
     subgraph linking it to the already-anchored region) is consistent
     between R and P.

  5. Atoms never anchored = reaction core. Bonds in R between mapped atoms
     whose mapped image isn't bonded in P = broken. Symmetrically for
     formed.
"""

from __future__ import annotations

import shutil
import subprocess
from collections import defaultdict, Counter
from pathlib import Path

import numpy as np
import networkx as nx
from networkx.algorithms.isomorphism import GraphMatcher


# -------------------- IO + xtb (same as before) --------------------

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
    workdir = Path(workdir); workdir.mkdir(parents=True, exist_ok=True)
    local = workdir / Path(xyz_path).name
    shutil.copy(xyz_path, local)
    cmd = ["xtb", local.name, "--gfn", "2", "--sp"]
    if charge: cmd += ["--chrg", str(charge)]
    if uhf: cmd += ["--uhf", str(uhf)]
    res = subprocess.run(cmd, cwd=workdir, capture_output=True, text=True)
    if res.returncode != 0:
        raise RuntimeError(f"xtb failed: {res.stderr[-500:]}")
    elements, coords = parse_xyz(local)
    n = len(elements)
    wbo = np.zeros((n, n))
    wf = workdir / "wbo"
    if not wf.exists():
        raise RuntimeError("no wbo file")
    for ln in wf.read_text().splitlines():
        parts = ln.split()
        if len(parts) < 3: continue
        i, j = int(parts[0]) - 1, int(parts[1]) - 1
        v = float(parts[2])
        wbo[i, j] = v; wbo[j, i] = v
    return elements, coords, wbo


# -------------------- WBO graph --------------------

def build_graph(elements, wbo, bond_cut=0.5):
    """
    Connectivity graph with WBO weights on edges. Bond exists iff
    WBO >= bond_cut (binary cutoff just for graph existence). The WBO
    value is kept on each edge and used during fragment matching with
    a tolerance of 0.5 (|delta WBO| > 0.5 -> not equivalent).
    """
    g = nx.Graph()
    for i, e in enumerate(elements):
        g.add_node(i, element=e)
    n = len(elements)
    for i in range(n):
        for j in range(i + 1, n):
            if wbo[i, j] >= bond_cut:
                g.add_edge(i, j, wbo=float(wbo[i, j]))
    return g


def k_ball(g, root, radius):
    """Return the subgraph induced by atoms within `radius` bonds of root."""
    levels = {root: 0}
    frontier = [root]
    for _ in range(radius):
        nxt = []
        for u in frontier:
            for v in g.neighbors(u):
                if v not in levels:
                    levels[v] = levels[u] + 1
                    nxt.append(v)
        frontier = nxt
    return g.subgraph(levels.keys())


# -------------------- fragment matching --------------------

def make_node_match(g_R, g_P):
    return lambda nR, nP: nR['element'] == nP['element']


def make_edge_match(wbo_tol):
    def em(eR, eP):
        return abs(eR['wbo'] - eP['wbo']) <= wbo_tol
    return em


def _frag_invariant(sub, root, wbo_round=0.5):
    """Cheap invariant for fast fragment-pair pre-filtering."""
    el_count = tuple(sorted(sub.nodes[n]['element'] for n in sub.nodes))
    deg_seq = tuple(sorted(d for _, d in sub.degree()))
    edge_hist = tuple(sorted(
        round(sub[u][v]['wbo'] / wbo_round) * wbo_round
        for u, v in sub.edges()
    ))
    root_neighbors = tuple(sorted(
        (round(sub[root][nbr]['wbo'] / wbo_round) * wbo_round,
         sub.nodes[nbr]['element'])
        for nbr in sub.neighbors(root)
    ))
    return (el_count, deg_seq, edge_hist, root_neighbors)


def match_fragment(g_R, g_P, root_R, root_P, radius, wbo_tol, max_isos=4):
    """
    All graph isomorphisms between the radius-r ball around root_R in g_R
    and the radius-r ball around root_P in g_P that map root_R to root_P.

    Pre-filters with cheap invariants (element multiset, degree sequence,
    edge WBO histogram, root-neighbor signature) to skip GraphMatcher
    invocations on obviously-incompatible fragment pairs -- the heavy
    isomorphism search only runs when the invariants already agree.
    """
    sub_R = k_ball(g_R, root_R, radius)
    sub_P = k_ball(g_P, root_P, radius)
    if sub_R.number_of_nodes() != sub_P.number_of_nodes():
        return []
    if sub_R.number_of_edges() != sub_P.number_of_edges():
        return []
    if _frag_invariant(sub_R, root_R) != _frag_invariant(sub_P, root_P):
        return []
    matcher = GraphMatcher(sub_R, sub_P,
                           node_match=make_node_match(g_R, g_P),
                           edge_match=make_edge_match(wbo_tol))
    found = []
    for m in matcher.isomorphisms_iter():
        if m.get(root_R) == root_P:
            found.append(dict(m))
            if len(found) >= max_isos:
                break
    return found


# -------------------- anchor finding --------------------

def find_fragment_anchors(g_R, g_P, max_radius=4, min_radius=1, wbo_tol=0.2,
                          min_fragment_size=2):
    """
    For each (u in R, v in P) with same element, find the largest radius r
    at which their r-balls are isomorphic with u<->v root. Stores ALL
    isomorphisms (not just one) so downstream merging can pick the
    correspondence that's consistent with already-anchored atoms.

    Filters out trivially-small fragments (size < min_fragment_size) since
    those match too many places to be useful anchors.
    """
    anchors = []
    n_R = g_R.number_of_nodes()
    n_P = g_P.number_of_nodes()
    for u in range(n_R):
        eR = g_R.nodes[u]['element']
        for v in range(n_P):
            if g_P.nodes[v]['element'] != eR:
                continue
            for r in range(max_radius, min_radius - 1, -1):
                isos = match_fragment(g_R, g_P, u, v, r, wbo_tol)
                if isos:
                    size = len(isos[0])
                    if size < min_fragment_size:
                        break
                    anchors.append({
                        'root_R': u, 'root_P': v, 'radius': r,
                        'isomorphisms': isos, 'size': size,
                        'internal_unique': len(isos) == 1,
                    })
                    break
    return anchors


def globally_unique_anchors(anchors):
    """
    Mark anchors with TWO uniqueness flags:
      root_unique:   root_R has only one P-partner across all anchors AND
                     root_P has only one R-partner -- so the root pair
                     (u <-> v) is definite even if the fragment has
                     internal symmetry.
      internal_unique: copy of `internal_unique` (only one iso of root pair).
      globally_unique: both root_unique AND internal_unique.

    Phase-1 of merging will use root_unique to lock just the root atoms,
    which is enough to start breaking symmetry of overlapping fragments.
    """
    by_root_R = defaultdict(set)
    by_root_P = defaultdict(set)
    for a in anchors:
        by_root_R[a['root_R']].add(a['root_P'])
        by_root_P[a['root_P']].add(a['root_R'])
    out = []
    for a in anchors:
        ru = (len(by_root_R[a['root_R']]) == 1
              and len(by_root_P[a['root_P']]) == 1)
        a = dict(a)
        a['root_unique'] = ru
        a['globally_unique'] = ru and a['internal_unique']
        out.append(a)
    return out


# -------------------- merging anchors into a global mapping --------------------

def merge_anchors(anchors, n_R, n_P, g_R=None, g_P=None):
    """
    Constraint propagation from globally-unique anchor seeds.

    Each anchor stores all isomorphisms of its fragment; we track which
    ones are still "alive" (consistent with the current mapping). Each
    commitment to mapping[r] = p eliminates inconsistent isos, which
    shrinks the candidate set of other atoms. Atoms whose candidate set
    collapses to exactly one P atom commit next. The cascade stops when
    no atom has a uniquely-determined candidate -- those atoms are the
    reaction core (or they're in a symmetric spectator region the
    cascade cannot disambiguate, which is fine: their bonds will appear
    as preserved at the multiset level even if individual atoms are
    interchangeable).
    """
    mapping = {}
    inv = {}
    alive_isos = [set(range(len(a['isomorphisms']))) for a in anchors]

    # PRE-SEED: any element appearing exactly once in both R and P. This
    # handles cases like a single Pd or single I where the atom's local
    # environment may change (so no fragment match), but the atom is still
    # uniquely identified by element alone.
    if g_R is not None and g_P is not None:
        el_count_R = Counter(g_R.nodes[i]['element'] for i in g_R.nodes)
        el_count_P = Counter(g_P.nodes[i]['element'] for i in g_P.nodes)
        for el, c in el_count_R.items():
            if c == 1 and el_count_P.get(el) == 1:
                r_atom = next(i for i in g_R.nodes if g_R.nodes[i]['element'] == el)
                p_atom = next(i for i in g_P.nodes if g_P.nodes[i]['element'] == el)
                if r_atom not in mapping and p_atom not in inv:
                    mapping[r_atom] = p_atom
                    inv[p_atom] = r_atom

    def candidate_set(r_atom):
        if r_atom in mapping:
            return {mapping[r_atom]}
        cs = set()
        for ai, a in enumerate(anchors):
            for idx in alive_isos[ai]:
                iso = a['isomorphisms'][idx]
                if r_atom in iso:
                    cs.add(iso[r_atom])
        return cs - set(inv.keys())

    def commit(r_atom, p_atom):
        mapping[r_atom] = p_atom
        inv[p_atom] = r_atom
        # Kill isos inconsistent with this commit
        for ai, a in enumerate(anchors):
            kill = set()
            for idx in alive_isos[ai]:
                iso = a['isomorphisms'][idx]
                if iso.get(r_atom, p_atom) != p_atom:
                    kill.add(idx)
                    continue
                # also check against full mapping
                for r, p in iso.items():
                    if r in mapping and mapping[r] != p:
                        kill.add(idx); break
                    if p in inv and inv[p] != r:
                        kill.add(idx); break
            alive_isos[ai] -= kill

    # SEED PHASE A: globally-unique fragment anchors -- the rock-solid
    # mounting points where both root pair AND iso are unique.
    for ai, a in enumerate(anchors):
        if not a['globally_unique']:
            continue
        iso = a['isomorphisms'][0]
        for r, p in iso.items():
            if r not in mapping and p not in inv:
                commit(r, p)

    # SEED PHASE B: largest-internally-unique-fragment first. Walk anchors
    # by descending size. Commit a fragment if its iso is consistent with
    # the current mapping (no conflicts) AND the iso has at least one
    # atom not yet mapped (so we make progress). This implements the
    # "find the biggest mounting island first, lock it all at once"
    # strategy: a single-atom commit would be at the mercy of local
    # Hungarian permutations, but committing a whole 8-atom rigid
    # fragment locks its internal correspondence completely.
    progressed = True
    while progressed:
        progressed = False
        for ai, a in sorted(enumerate(anchors),
                            key=lambda x: -x[1]['size']):
            if not a['internal_unique']:
                continue
            if not alive_isos[ai]:
                continue
            iso = a['isomorphisms'][next(iter(alive_isos[ai]))]
            # Skip if iso conflicts with current mapping
            ok = True
            for r, p in iso.items():
                if r in mapping and mapping[r] != p:
                    ok = False; break
                if p in inv and inv[p] != r:
                    ok = False; break
            if not ok:
                continue
            # Skip if iso adds nothing
            new_atoms = [(r, p) for r, p in iso.items()
                         if r not in mapping and p not in inv]
            if not new_atoms:
                continue
            # Require that the (root_R, root_P) pair is locally distinguished:
            # either it's already in mapping, or this anchor's root has the
            # FEWEST P-candidates among remaining anchors that mention root_R.
            # Cheap heuristic: only commit if root pair is currently unique
            # OR fragment size >= some big threshold (8+ atoms).
            if a['size'] >= 6 or a['root_R'] in mapping or a['root_P'] in inv:
                for r, p in new_atoms:
                    commit(r, p)
                progressed = True

    # PROPAGATION: any atom whose alive isos all agree on its image
    while True:
        progressed = False
        # Find atom with candidate set of size 1 that isn't mapped yet
        for r in range(n_R):
            if r in mapping:
                continue
            cs = candidate_set(r)
            if len(cs) == 1:
                p = next(iter(cs))
                if p not in inv:
                    commit(r, p)
                    progressed = True
                    break
        if not progressed:
            break

    # SECOND PROPAGATION: also use isos that locally agree even if root pair
    # wasn't initially globally unique. After the seed cascade, many fragment
    # isos have only one alive option among those overlapping the mapping --
    # commit those too.
    while True:
        progressed = False
        for ai, a in enumerate(anchors):
            if len(alive_isos[ai]) != 1:
                continue
            iso = a['isomorphisms'][next(iter(alive_isos[ai]))]
            # Require some overlap with current mapping (so it's not floating)
            overlap = sum(1 for r, p in iso.items()
                          if r in mapping and mapping[r] == p)
            if overlap == 0:
                continue
            for r, p in iso.items():
                if r in mapping or p in inv:
                    continue
                commit(r, p)
                progressed = True
        if not progressed:
            break

    return mapping


# -------------------- expansion + core detection --------------------

def expand_mapping(mapping, g_R, g_P, max_passes=20):
    """
    Pure-connectivity expansion from already-mapped atoms. For each mapped
    pair (u, v), pair u's unmapped R-neighbors with v's unmapped
    P-neighbors element-by-element. If counts match exactly, commit the
    pairing in arbitrary order (symmetric atoms like methyl Hs are
    interchangeable). If counts differ, leave them unmapped -- that's a
    real connectivity change at this atom, i.e. reaction-core.
    """
    inv = {v: k for k, v in mapping.items()}
    for _ in range(max_passes):
        progressed = False
        for u in list(mapping.keys()):
            v = mapping[u]
            r_groups = defaultdict(list)
            for w in g_R.neighbors(u):
                if w in mapping:
                    continue
                r_groups[g_R.nodes[w]['element']].append(w)
            p_groups = defaultdict(list)
            for x in g_P.neighbors(v):
                if x in inv:
                    continue
                p_groups[g_P.nodes[x]['element']].append(x)
            for el, rs in r_groups.items():
                ps = p_groups.get(el, [])
                if len(ps) != len(rs):
                    continue  # connectivity change here, leave unmapped
                for w, x in zip(rs, ps):
                    mapping[w] = x
                    inv[x] = w
                    progressed = True
        if not progressed:
            break
    return mapping


# -------------------- grow-until-unique island finder --------------------

def _set_unique(cands):
    """True iff every candidate iso covers the same SET of P atoms.
    Symmetry-equivalent permutations (different sequences, same set)
    count as unique. Returns False if cands is empty."""
    if not cands:
        return False
    if len(cands) == 1:
        return True
    s0 = frozenset(cands[0].values())
    return all(frozenset(c.values()) == s0 for c in cands[1:])


def _filter_cands_by_cross_bonds(cands, fragment, mapping, g_R, g_P, wbo_tol):
    """Drop any cand whose specific atom assignments would CONFLICT with
    the cross-bonds connecting fragment atoms to atoms in existing
    locked islands.

    For each cand and each (r in fragment, r' not in fragment but in
    mapping, R-bond between r-r'), require that cand[r] is bonded to
    mapping[r'] in g_P with |dWBO| <= wbo_tol. Cands that fail are
    chemically inconsistent with the existing islands and must be
    discarded.

    Returns the filtered list. If empty, the fragment is genuinely
    isolated from the existing islands at the chemistry boundary.
    """
    out = []
    for cand in cands:
        ok = True
        for r in fragment:
            p = cand[r]
            for r2 in g_R.neighbors(r):
                if r2 in fragment:
                    continue
                if r2 not in mapping:
                    continue
                p2 = mapping[r2]
                wR = g_R[r][r2]['wbo']
                if not g_P.has_edge(p, p2):
                    ok = False; break
                wP = g_P[p][p2]['wbo']
                if abs(wR - wP) > wbo_tol:
                    ok = False; break
            if not ok:
                break
        if ok:
            out.append(cand)
    return out


def _absorb_mapped_frontier(g_R, g_P, fragment, candidates, mapping, inv,
                             wbo_tol, distance=None):
    """Repeatedly absorb any frontier atom that is already in `mapping`
    (a previously-locked atom). Each absorption pins that atom's image
    to mapping[n] and narrows cands accordingly. Returns updated
    (fragment, candidates). Stops when no more mapped frontier atoms
    can be absorbed."""
    while True:
        absorbed = False
        # Frontier: atoms outside fragment, bonded to something inside.
        for u in list(fragment):
            for n in g_R.neighbors(u):
                if n in fragment:
                    continue
                if n not in mapping:
                    continue
                # Try to absorb n with pinned image mapping[n]
                bonded = [b for b in g_R.neighbors(n) if b in fragment]
                r_wbos = [(b, g_R[b][n]['wbo']) for b in bonded]
                n_pinned = mapping[n]
                n_el = g_R.nodes[n]['element']
                if g_P.nodes[n_pinned]['element'] != n_el:
                    continue  # element mismatch — shouldn't happen
                new_cands = []
                for cand in candidates:
                    used_p = set(cand.values())
                    if n_pinned in used_p:
                        continue
                    # Check that n_pinned is a P-neighbor of cand[b] for all b
                    ok = True
                    for (b, wR) in r_wbos:
                        if not g_P.has_edge(cand[b], n_pinned):
                            ok = False; break
                        if abs(wR - g_P[cand[b]][n_pinned]['wbo']) > wbo_tol:
                            ok = False; break
                    if ok:
                        nc = dict(cand); nc[n] = n_pinned
                        new_cands.append(nc)
                if new_cands:
                    fragment.add(n)
                    candidates = new_cands
                    if distance is not None:
                        distance[n] = 1 + min(distance[b] for b in bonded
                                               if b in distance)
                    absorbed = True
                    break
            if absorbed:
                break
        if not absorbed:
            break
    return fragment, candidates


def grow_island(g_R, g_P, seed, mapping, inv,
                wbo_tol=0.5,
                growth_min_wbo=0.6,
                top_degen=0.1,
                max_lock_cands=100,
                min_lock_size=2):
    """
    Grow a fragment from `seed` in g_R and lock it as an island when it
    matches uniquely (or with a small symmetry-equivalent set) in g_P.

    --- Inputs ---
    g_R, g_P     : NetworkX undirected graphs of R and P (nodes have
                   `element`, edges have `wbo`).
    seed         : the R atom to start from (int).
    mapping, inv : current global R->P / P->R maps. Atoms already in
                   `mapping` participate in soft-merge: the iso must
                   respect their existing image when growth absorbs them.

    --- Algorithm parameters ---
    wbo_tol         (0.5) : edge-WBO match tolerance. An iso is valid
                            iff every fragment bond's |WBO_R - WBO_P|
                            is within this. Anything bigger -> chemistry
                            (broken/formed bond).
    growth_min_wbo  (0.6) : weak-bond cutoff. The fragment refuses to
                            traverse any bond below this WBO.
    top_degen       (0.1) : at each step, only frontier atoms whose
                            connecting WBO is within `top_degen` of the
                            shell's max-WBO are considered (top-1 with
                            small-tie tolerance for chemically-equivalent
                            bonds).
    max_lock_cands  (8)   : how many candidate isos we still accept as
                            "symmetry-degenerate" at growth halt. >8 is
                            treated as too ambiguous -> reject.
    min_lock_size   (2)   : require the fragment to have grown beyond
                            just the seed before it can lock. Singletons
                            (size 1) are always rejected.

    --- Termination ---
    Returns the iso (dict R->P) when:
      * cands == 1, or
      * growth halts (frontier exhausted / no strong frontier / no
        extension keeps any cand alive) AND cands <= max_lock_cands AND
        fragment_size >= min_lock_size.
    Returns None otherwise (chemistry-boundary or too-ambiguous seed).
    """
    # Hardcoded safety caps (rarely tweaked).
    max_cands_hard = 2000  # cap during growth to prevent memory blowup
    max_iters = g_R.number_of_nodes()  # outer-loop safety, rarely tight

    # ... (after seed init below, before main loop, we'll do an eager
    # mapped-frontier absorption at every step)
    if seed in mapping:
        return None
    seed_el = g_R.nodes[seed]['element']
    candidates = [{seed: v} for v in g_P.nodes()
                  if v not in inv and g_P.nodes[v]['element'] == seed_el]
    if not candidates:
        return None
    fragment = {seed}
    distance = {seed: 0}  # graph distance from seed within fragment
    for _ in range(max_iters):
        # PHASE A: eagerly absorb any frontier atoms that are already in
        # `mapping` (previously-locked island atoms). Each absorption
        # pins its image and narrows cands.
        fragment, candidates = _absorb_mapped_frontier(
            g_R, g_P, fragment, candidates, mapping, inv, wbo_tol,
            distance=distance)
        if not candidates:
            return None
        # Set-uniqueness: lock when all isos describe the same P-atom set
        # (symmetry-induced sequence permutations are accepted).
        if _set_unique(candidates):
            return candidates[0]
        frontier = set()
        for u in fragment:
            for n in g_R.neighbors(u):
                if n not in fragment:
                    frontier.add(n)
        if not frontier:
            # Whole connected component absorbed. Lock if set-unique
            # (all cands cover same P-atom set) and fragment grew.
            if _set_unique(candidates) and len(fragment) >= min_lock_size:
                return candidates[0]
            return None

        # Only consider strong-bonded frontier atoms.
        frontier_info = {}
        for n in frontier:
            bonded = [u for u in g_R.neighbors(n) if u in fragment]
            max_w = max(g_R[u][n]['wbo'] for u in bonded)
            if max_w < growth_min_wbo:
                continue
            # Distance from seed = 1 + min distance of any bonded fragment atom
            dist = 1 + min(distance[u] for u in bonded if u in distance)
            frontier_info[n] = (dist, max_w)
        if not frontier_info:
            # Only weak (<0.6) frontier remaining. Lock if set-unique.
            if _set_unique(candidates) and len(fragment) >= min_lock_size:
                return candidates[0]
            return None

        # BFS-by-distance-from-seed FIRST: grow seed's direct neighbors
        # before going to second-shell atoms. Within the same distance
        # shell, apply top-1 WBO with `top_degen` tolerance.
        min_dist = min(d for (d, _) in frontier_info.values())
        same_shell = {n: w for n, (d, w) in frontier_info.items() if d == min_dist}
        top_w = max(same_shell.values())
        top_band = [n for n, w in same_shell.items() if w >= top_w - top_degen]

        best_n = None
        best_cands = None
        for n in top_band:
            n_el = g_R.nodes[n]['element']
            bonded_in_frag = [u for u in g_R.neighbors(n) if u in fragment]
            r_wbos = [(u, g_R[u][n]['wbo']) for u in bonded_in_frag]
            n_pinned = mapping.get(n, None)  # soft merge constraint
            new_cands = []
            over = False
            for cand in candidates:
                used_p = set(cand.values())
                v_set = set(g_P.neighbors(cand[bonded_in_frag[0]]))
                for u in bonded_in_frag[1:]:
                    v_set &= set(g_P.neighbors(cand[u]))
                v_set -= used_p
                if n_pinned is not None:
                    v_set = v_set & {n_pinned}
                for v in v_set:
                    if g_P.nodes[v]['element'] != n_el:
                        continue
                    if all(abs(w - g_P[cand[u]][v]['wbo']) <= wbo_tol
                           for u, w in r_wbos):
                        nc = dict(cand); nc[n] = v
                        new_cands.append(nc)
                        if len(new_cands) > max_cands_hard:
                            over = True; break
                if over: break
            if over or not new_cands:
                continue
            if best_cands is None or len(new_cands) < len(best_cands):
                best_n = n
                best_cands = new_cands
        if best_n is None:
            # No extension keeps any candidate alive. Lock if set-unique.
            if _set_unique(candidates) and len(fragment) >= min_lock_size:
                return candidates[0]
            return None
        fragment.add(best_n)
        # Update distance for the newly committed atom
        bonded_in_frag = [u for u in g_R.neighbors(best_n) if u in fragment - {best_n}]
        distance[best_n] = 1 + min(distance[u] for u in bonded_in_frag)
        candidates = best_cands
    # End of max_iters loop
    if _set_unique(candidates) and len(fragment) >= min_lock_size:
        return candidates[0]
    return None


def merge_touching_islands(g_R, g_P, mapping, atom_island_R, atom_island_P, wbo_tol=0.5):
    """
    After all seeds are tried, find every pair of islands (A, B) that
    share at least one direct edge in g_R (i.e., an atom in A is bonded
    to an atom in B). For each such pair, verify the merge is valid:
    every cross-island bond R[a]-R[b] (a in A, b in B) must have a
    corresponding bond in g_P between mapping[a] and mapping[b] with
    |WBO_R - WBO_P| <= wbo_tol. If valid -> merge by relabeling B's
    atoms into A's island id. If invalid -> keep separate (the boundary
    is reaction chemistry).
    Iterates until no more merges happen (A merges into B may enable
    a transitive (A+B)-C merge).
    """
    progressed = True
    while progressed:
        progressed = False
        # Find all (A, B) pairs that touch (and aren't the same island)
        seen = set()
        for u, v in g_R.edges():
            iA = atom_island_R.get(u)
            iB = atom_island_R.get(v)
            if iA is None or iB is None or iA == iB:
                continue
            pair = (min(iA, iB), max(iA, iB))
            if pair in seen:
                continue
            seen.add(pair)
            # Validate merge: every cross-island bond must agree on WBO
            ok = True
            for x, y in g_R.edges():
                ixA = atom_island_R.get(x)
                ixB = atom_island_R.get(y)
                if {ixA, ixB} != {iA, iB}:
                    continue
                wR = g_R[x][y]['wbo']
                px, py = mapping[x], mapping[y]
                if not g_P.has_edge(px, py):
                    ok = False; break
                wP = g_P[px][py]['wbo']
                if abs(wR - wP) > wbo_tol:
                    ok = False; break
            if not ok:
                continue
            # Merge: pick the smaller id and relabel
            keep, drop = pair
            for r, idx in list(atom_island_R.items()):
                if idx == drop:
                    atom_island_R[r] = keep
                    atom_island_P[mapping[r]] = keep
            progressed = True
            break
    return mapping, atom_island_R, atom_island_P


def find_islands(g_R, g_P, wbo_tol=0.5):
    """
    Iteratively grow islands with soft merge. Each pass walks every
    unmapped seed; if grow_island returns a unique iso, commit it
    (which may merge with already-locked islands automatically because
    grow_island treats mapped neighbors as constraints). Stop when a
    pass commits nothing. Atoms never mapped are the reaction core.
    """
    mapping = {}
    inv = {}
    n_islands = 0
    while True:
        progressed = False
        for seed in sorted(g_R.nodes()):
            if seed in mapping:
                continue
            iso = grow_island(g_R, g_P, seed, mapping, inv, wbo_tol=wbo_tol)
            if iso is None:
                continue
            for r, p in iso.items():
                if r not in mapping and p not in inv:
                    mapping[r] = p
                    inv[p] = r
            n_islands += 1
            progressed = True
        if not progressed:
            break
    return mapping, n_islands


# -------------------- fallback: 3D-distance Hungarian for residuals --------------------

def kabsch(P, Q):
    Pc, Qc = P.mean(0), Q.mean(0)
    H = (P - Pc).T @ (Q - Qc)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1, 1, d])
    R = Vt.T @ D @ U.T
    return R, Qc - R @ Pc


def map_residual_by_3d(mapping, coords_R, coords_P, elements_R, elements_P, k_local=8):
    """
    For atoms unmapped after fragment matching + expansion: use a local
    Kabsch built from the k nearest already-mapped neighbors to project
    each unmapped R atom into P-space, then Hungarian on 3D distance per
    element. Strictly a visualization fallback -- atoms mapped here are
    still considered part of the reaction core.
    """
    from scipy.optimize import linear_sum_assignment
    if not mapping:
        return mapping, set()
    inv = {v: k for k, v in mapping.items()}
    unmapped_R = [i for i in range(len(elements_R)) if i not in mapping]
    unmapped_P = [j for j in range(len(elements_P)) if j not in inv]
    if not unmapped_R or not unmapped_P:
        return mapping, set()
    mapped_R_idx = np.array(sorted(mapping.keys()))
    if len(mapped_R_idx) < 3:
        return mapping, set()
    by_el_R = defaultdict(list); by_el_P = defaultdict(list)
    for i in unmapped_R: by_el_R[elements_R[i]].append(i)
    for j in unmapped_P: by_el_P[elements_P[j]].append(j)
    fallback_atoms = set()
    for el, rs in by_el_R.items():
        ps = by_el_P.get(el, [])
        if not ps: continue
        projections = []
        for i in rs:
            d = np.linalg.norm(coords_R[mapped_R_idx] - coords_R[i], axis=1)
            k = min(k_local, len(mapped_R_idx))
            order = np.argsort(d)[:k]
            nn_R = mapped_R_idx[order]
            nn_P = np.array([mapping[int(x)] for x in nn_R])
            Rmat, tvec = kabsch(coords_R[nn_R], coords_P[nn_P])
            projections.append(Rmat @ coords_R[i] + tvec)
        cost = np.zeros((len(rs), len(ps)))
        for a, proj in enumerate(projections):
            for b, j in enumerate(ps):
                cost[a, b] = np.linalg.norm(proj - coords_P[j])
        row, col = linear_sum_assignment(cost)
        for a, b in zip(row, col):
            mapping[rs[a]] = ps[b]
            inv[ps[b]] = rs[a]
            fallback_atoms.add(rs[a])
    return mapping, fallback_atoms


# -------------------- bond classification --------------------

def classify_bonds(mapping, wbo_R, wbo_P, bond_high=0.6, bond_low=0.3):
    """
    Hysteresis-style classification to avoid false positives from
    WBO drift across a single threshold:
      broken iff WBO_R >= bond_high and WBO_P < bond_low
      formed iff WBO_R < bond_low and WBO_P >= bond_high
    Bonds whose WBO sits in the [bond_low, bond_high) gray zone on either
    side are treated as preserved partial bonds, not flagged.
    """
    inv = {v: k for k, v in mapping.items()}
    nR, nP = wbo_R.shape[0], wbo_P.shape[0]
    broken, formed = [], []
    for i in range(nR):
        for j in range(i + 1, nR):
            wR = wbo_R[i, j]
            if wR < bond_high: continue
            if i not in mapping or j not in mapping:
                broken.append((i, j, float(wR), None)); continue
            wP = wbo_P[mapping[i], mapping[j]]
            if wP < bond_low:
                broken.append((i, j, float(wR), float(wP)))
    for ip in range(nP):
        for jp in range(ip + 1, nP):
            wP = wbo_P[ip, jp]
            if wP < bond_high: continue
            if ip not in inv or jp not in inv:
                formed.append((ip, jp, None, float(wP))); continue
            wR = wbo_R[inv[ip], inv[jp]]
            if wR < bond_low:
                formed.append((ip, jp, float(wR), float(wP)))
    core_R = set()
    core_P = set()
    for (i, j, _, _) in broken: core_R.add(i); core_R.add(j)
    for (ip, jp, _, _) in formed: core_P.add(ip); core_P.add(jp)
    return broken, formed, sorted(core_R), sorted(core_P)


# -------------------- driver --------------------

def analyze(reactant_xyz, product_xyz, workdir,
            charge=0, uhf=0,
            graph_bond_cut=0.5, bond_high=0.6, bond_low=0.3):
    workdir = Path(workdir)
    elR, xyzR, wboR = run_xtb(reactant_xyz, workdir / "R", charge=charge, uhf=uhf)
    elP, xyzP, wboP = run_xtb(product_xyz, workdir / "P", charge=charge, uhf=uhf)
    if Counter(elR) != Counter(elP):
        raise ValueError(f"composition mismatch: {Counter(elR)} vs {Counter(elP)}")

    g_R = build_graph(elR, wboR, bond_cut=graph_bond_cut)
    g_P = build_graph(elP, wboP, bond_cut=graph_bond_cut)

    # Strict grow-until-unique with WBO tol=0.5 and soft merge.
    mapping, n_islands = find_islands(g_R, g_P, wbo_tol=0.5)
    n_after_merge = len(mapping)
    # Element-counted expansion for symmetric neighbors of mapped atoms.
    mapping = expand_mapping(mapping, g_R, g_P)
    n_after_expand = len(mapping)
    fallback_atoms = set()
    anchors = []

    broken, formed, core_R, core_P = classify_bonds(mapping, wboR, wboP,
                                                    bond_high=bond_high, bond_low=bond_low)
    # core_R / core_P now contains ONLY atoms touching a broken/formed bond.
    # Fallback-mapped atoms whose bonds are preserved aren't core. (We keep
    # fallback_atoms in the result dict for inspection if needed.)

    return dict(
        elements_R=elR, coords_R=xyzR, wbo_R=wboR,
        elements_P=elP, coords_P=xyzP, wbo_P=wboP,
        anchors=anchors, mapping=mapping,
        broken=broken, formed=formed,
        core_R=core_R, core_P=core_P,
        n_after_merge=n_after_merge,
        n_after_expand=n_after_expand,
        n_anchors=len(anchors),
        fallback_atoms=sorted(fallback_atoms),
    )
