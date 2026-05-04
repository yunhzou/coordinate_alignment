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

def expand_mapping(mapping, g_R, g_P,
                   events=None, atom_island_R=None, atom_island_P=None):
    """
    Pure-connectivity expansion from already-mapped atoms. For each mapped
    pair (u, v), pair u's unmapped R-neighbors with v's unmapped
    P-neighbors element-by-element. If counts match exactly, commit the
    pairing in arbitrary order (symmetric atoms like methyl Hs are
    interchangeable). If counts differ, leave them unmapped -- that's a
    real connectivity change at this atom, i.e. reaction-core.
    Loops until no further progress is made.

    Optional trace recording:
      events         : list to receive island_locked events (one per
                       (parent, element) group).
      atom_island_R  : dict to receive R-atom -> island_id map. Newly
                       paired atoms inherit the parent's island.
      atom_island_P  : dict to receive P-atom -> island_id map.
    """
    record = events is not None
    track = atom_island_R is not None and atom_island_P is not None
    inv = {v: k for k, v in mapping.items()}
    pass_no = 0
    while True:
        progressed = False
        pass_no += 1
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
                if track:
                    parent_island = atom_island_R.get(u)
                    if parent_island is None:
                        parent_island = _next_island_id(atom_island_R)
                        atom_island_R[u] = parent_island
                        atom_island_P[v] = parent_island
                paired = []
                for w, x in zip(rs, ps):
                    mapping[w] = x
                    inv[x] = w
                    progressed = True
                    if track:
                        atom_island_R[w] = parent_island
                        atom_island_P[x] = parent_island
                    if record:
                        paired.append((int(w), int(x)))
                if record and paired:
                    events.append({
                        'type': 'island_locked',
                        'island_idx': parent_island if track else None,
                        'pairs': paired,
                        'merged_with': [],
                        'relabeled': [],
                        'mapped_total': len(mapping),
                        'expand_pass': pass_no,
                        'parent_atom': int(u),
                    })
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


def _p_atoms_from_cands(cands):
    """Union of P atoms across all candidate isos. Used for trace events
    so the visualizer can color every P-atom currently still-possible as
    an image of the growing fragment."""
    out = set()
    for c in cands:
        for v in c.values():
            out.add(int(v))
    return sorted(out)


def _next_island_id(atom_island_R):
    """Smallest unused island id. Used by find_islands and expand_mapping
    when allocating ids during trace recording."""
    return max(atom_island_R.values(), default=0) + 1


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
                min_lock_size=2,
                events=None):
    """
    Grow a fragment from `seed` in g_R and lock it as an island when
    growth halts. Any cand surviving the chemistry-boundary stop has
    the same fragment and the same WBO-consistent assignment up to
    symmetric atoms (e.g. methyl Hs); picking candidates[0] is safe
    because all surviving cands are chemically equivalent. No cand-
    count cap.

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
    min_lock_size   (2)   : minimum fragment size for a lock. Default 2
                            prevents singleton seeds from polluting the
                            mapping during main growth (a low-idx singleton
                            could steal symmetric atoms that a higher-idx
                            big seed needs for its iso). Truly leftover
                            atoms are paired by `residual_pair` after
                            expand_mapping.
    events          (None): if a list is passed, every seed_start /
                            commit / seed_end decision is appended to it
                            for visualization. Behavior is identical
                            with or without recording.

    --- Termination ---
    Returns the iso (dict R->P) when:
      * cands are set-unique (early lock), or
      * growth halts (frontier exhausted / no strong frontier / no
        extension keeps any cand alive) AND fragment_size >= min_lock_size.
    Returns None otherwise (chemistry boundary kicked all cands).
    """
    max_cands_hard = 2000
    max_iters = g_R.number_of_nodes()
    record = events is not None

    if seed in mapping:
        return None
    seed_el = g_R.nodes[seed]['element']
    candidates = [{seed: v} for v in g_P.nodes()
                  if v not in inv and g_P.nodes[v]['element'] == seed_el]
    if not candidates:
        if record:
            events.append({'type': 'seed_start', 'seed': seed,
                           'init_cands': 0, 'p_atoms': []})
            events.append({'type': 'seed_end', 'result': 'no_initial_cands'})
        return None
    fragment = {seed}
    distance = {seed: 0}
    if record:
        events.append({
            'type': 'seed_start',
            'seed': seed,
            'init_cands': len(candidates),
            'fragment': sorted(fragment),
            'p_atoms': _p_atoms_from_cands(candidates),
        })

    def _try_lock(reason):
        # POLICY: when growth has halted (frontier exhausted, no strong
        # frontier, or all extensions cut) and cands remain non-empty,
        # commit candidates[0]. The cands all share the same fragment
        # and the same WBO-consistent assignment; any remaining
        # ambiguity is symmetry over interchangeable atoms (methyl Hs,
        # equivalent ligands), where every cand gives the same bond-
        # classification result. We do not score / tiebreak — pick any.
        ok = bool(candidates) and len(fragment) >= min_lock_size
        if record:
            events.append({
                'type': 'seed_end',
                'result': 'success' if ok else reason,
                'final_cands': len(candidates),
                'fragment': sorted(fragment),
                'iso': ({int(k): int(v) for k, v in candidates[0].items()}
                        if ok else None),
            })
        return candidates[0] if ok else None

    for _ in range(max_iters):
        # PHASE A: eagerly absorb any frontier atoms already in `mapping`
        # (previously-locked island atoms). Each absorption pins its
        # image and narrows cands.
        fragment, candidates = _absorb_mapped_frontier(
            g_R, g_P, fragment, candidates, mapping, inv, wbo_tol,
            distance=distance)
        if not candidates:
            if record:
                events.append({
                    'type': 'seed_end', 'result': 'absorb_failed',
                    'final_cands': 0,
                    'fragment': sorted(fragment),
                    'iso': None,
                })
            return None
        # Set-uniqueness: lock when all isos cover the same P-atom set
        # (symmetry-induced sequence permutations accepted).
        if _set_unique(candidates):
            if record:
                events.append({
                    'type': 'seed_end', 'result': 'success',
                    'final_cands': len(candidates),
                    'fragment': sorted(fragment),
                    'iso': {int(k): int(v) for k, v in candidates[0].items()},
                })
            return candidates[0]
        frontier = set()
        for u in fragment:
            for n in g_R.neighbors(u):
                if n not in fragment:
                    frontier.add(n)
        if not frontier:
            return _try_lock('no_frontier')

        # Only consider strong-bonded frontier atoms.
        frontier_info = {}
        filter_reason = {}
        for n in frontier:
            bonded = [u for u in g_R.neighbors(n) if u in fragment]
            max_w = max(g_R[u][n]['wbo'] for u in bonded)
            dist = 1 + min(distance[u] for u in bonded if u in distance)
            if max_w < growth_min_wbo:
                if record:
                    filter_reason[n] = (f'wbo<{growth_min_wbo}', max_w, dist)
                continue
            frontier_info[n] = (dist, max_w)
        if not frontier_info:
            return _try_lock('no_strong_frontier')

        # BFS-by-distance-from-seed PREFERRED: try closest shell first,
        # but if no extension works there, fall back to further shells
        # (any atom in fragment can propagate). "If 16 ways, keep
        # propagating from any atom in island."
        shells_sorted = sorted(set(d for (d, _) in frontier_info.values()))
        best_n = None
        best_cands = None
        tries = []
        chosen_shell = None
        chosen_top_w = None
        for cur_dist in shells_sorted:
            same_shell = {n: w for n, (d, w) in frontier_info.items() if d == cur_dist}
            top_w = max(same_shell.values())
            tier_top = [n for n, w in same_shell.items() if w >= top_w - top_degen]
            tier_rest = [n for n, w in same_shell.items() if w < top_w - top_degen]
            # Two tiers per shell: top-WBO band first, then lower-WBO
            # fallback. If top band all cut (e.g. chemistry-boundary
            # high-WBO bond like a C=C → C-C), try other propagation
            # directions from the fragment before giving up the shell.
            for tier_atoms in (tier_top, tier_rest):
                for n in tier_atoms:
                    n_el = g_R.nodes[n]['element']
                    bonded_in_frag = [u for u in g_R.neighbors(n) if u in fragment]
                    r_wbos = [(u, g_R[u][n]['wbo']) for u in bonded_in_frag]
                    n_pinned = mapping.get(n, None)
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
                    if record:
                        decision = 'CUT' if (over or not new_cands) else 'ok'
                        wbo_str = ', '.join(f'{u}:{w:.2f}' for u, w in r_wbos)
                        max_w_to_frag = max(w for _, w in r_wbos) if r_wbos else 0.0
                        min_d = 1 + min(distance[u] for u, _ in r_wbos
                                        if u in distance)
                        tries.append({
                            'atom': int(n),
                            'element': g_R.nodes[n]['element'],
                            'new_cands': len(new_cands),
                            'over': over,
                            'decision': decision,
                            'max_wbo_to_frag': round(max_w_to_frag, 3),
                            'wbo_bonds': wbo_str,
                            'distance_from_seed': min_d,
                            'shell_attempted': cur_dist,
                        })
                    if over or not new_cands:
                        continue
                    if best_cands is None or len(new_cands) < len(best_cands):
                        best_n = n
                        best_cands = new_cands
                        chosen_shell = cur_dist
                        chosen_top_w = top_w
                if best_n is not None:
                    break  # extension found at this tier -- commit it
            if best_n is not None:
                break  # extension found at this shell -- commit it

        if best_n is None:
            return _try_lock('all_cut')

        # Build step_info for this commit BEFORE we extend.
        step_info = None
        if record:
            if chosen_shell is None:
                chosen_shell = shells_sorted[0] if shells_sorted else 0
                chosen_top_w = max((w for (_, w) in frontier_info.values()),
                                   default=0.0)
            filtered = [{
                'atom': int(n),
                'element': g_R.nodes[n]['element'],
                'max_wbo_to_frag': round(w, 3),
                'distance_from_seed': d,
                'filtered_reason': reason,
            } for n, (reason, w, d) in filter_reason.items()]
            step_info = {
                'cands_before': len(candidates),
                'shell': chosen_shell,
                'top_wbo': round(chosen_top_w, 3),
                'tried': tries,
                'filtered': filtered,
            }

        fragment.add(best_n)
        bonded_in_frag = [u for u in g_R.neighbors(best_n)
                          if u in fragment - {best_n}]
        distance[best_n] = 1 + min(distance[u] for u in bonded_in_frag)
        candidates = best_cands
        if record:
            commit_bonds = [(u, round(g_R[u][best_n]['wbo'], 3))
                            for u in bonded_in_frag]
            events.append({
                'type': 'commit',
                'added': int(best_n),
                'element': g_R.nodes[best_n]['element'],
                'cands': len(candidates),
                'fragment': sorted(fragment),
                'p_atoms': _p_atoms_from_cands(candidates),
                'distance_from_seed': distance[best_n],
                'bonds_to_fragment': commit_bonds,
                'step_info': step_info,
            })

    return _try_lock('max_iters')


def merge_touching_islands(g_R, g_P, mapping, atom_island_R, atom_island_P,
                           wbo_tol=0.5, events=None):
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

    Optional `events` list receives island_locked events with
    island_island_merge=True and cross_bonds metadata.
    """
    record = events is not None
    progressed = True
    while progressed:
        progressed = False
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
            ok = True
            cross_bonds = [] if record else None
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
                if record:
                    cross_bonds.append((int(x), int(y),
                                        round(wR, 3), round(wP, 3)))
            if not ok:
                continue
            keep, drop = pair
            relabeled = []
            for r, idx in list(atom_island_R.items()):
                if idx == drop:
                    if record:
                        relabeled.append((int(r), int(drop)))
                    atom_island_R[r] = keep
                    atom_island_P[mapping[r]] = keep
            if record:
                events.append({
                    'type': 'island_locked',
                    'island_idx': keep,
                    'pairs': [],
                    'merged_with': [drop],
                    'relabeled': relabeled,
                    'mapped_total': len(mapping),
                    'island_island_merge': True,
                    'cross_bonds': cross_bonds,
                })
            progressed = True
            break
    return mapping, atom_island_R, atom_island_P


def find_islands(g_R, g_P, wbo_tol=0.5,
                 events=None, atom_island_R=None, atom_island_P=None):
    """
    Iteratively grow islands with soft merge. Each pass walks every
    unmapped seed; if grow_island returns a unique iso, commit it
    (which may merge with already-locked islands automatically because
    grow_island treats mapped neighbors as constraints). Stop when a
    pass commits nothing. Atoms never mapped are the reaction core.

    Optional trace recording:
      events         : list to receive pass_start / seed_* / commit /
                       island_locked events.
      atom_island_R  : dict to receive R-atom -> island_id map.
      atom_island_P  : dict to receive P-atom -> island_id map.
    All three are independent — pass any subset.
    """
    record = events is not None
    track = atom_island_R is not None and atom_island_P is not None
    mapping = {}
    inv = {}
    n_islands = 0
    pass_no = 0
    while True:
        progressed = False
        pass_no += 1
        if record:
            events.append({'type': 'pass_start', 'pass': pass_no,
                           'mapped': len(mapping)})
        for seed in sorted(g_R.nodes()):
            if seed in mapping:
                continue
            iso = grow_island(g_R, g_P, seed, mapping, inv,
                              wbo_tol=wbo_tol, events=events)
            if iso is None:
                continue

            # Determine which existing islands this iso touches and pick
            # the surviving id. When not tracking, just allocate fresh.
            if track:
                touched = set()
                for r in iso.keys():
                    if r in atom_island_R:
                        touched.add(atom_island_R[r])
                if touched:
                    merged_id = min(touched)
                else:
                    merged_id = _next_island_id(atom_island_R)
            else:
                touched = set()
                merged_id = None

            committed_new = []
            relabeled = []
            for r, p in iso.items():
                if r not in mapping:
                    mapping[r] = p
                    inv[p] = r
                    if track:
                        atom_island_R[r] = merged_id
                        atom_island_P[p] = merged_id
                    if record:
                        committed_new.append((int(r), int(p)))
                elif track and atom_island_R.get(r) != merged_id:
                    if record:
                        relabeled.append((int(r), int(atom_island_R[r])))
                    atom_island_R[r] = merged_id
                    atom_island_P[mapping[r]] = merged_id
            if track:
                for r in list(atom_island_R.keys()):
                    if (atom_island_R[r] in touched
                            and atom_island_R[r] != merged_id):
                        if record:
                            relabeled.append((int(r), int(atom_island_R[r])))
                        atom_island_R[r] = merged_id
                        atom_island_P[mapping[r]] = merged_id

            n_islands += 1
            progressed = True
            if record:
                events.append({
                    'type': 'island_locked',
                    'island_idx': merged_id if track else n_islands,
                    'pairs': committed_new,
                    'merged_with': (sorted(touched - {merged_id})
                                    if track else []),
                    'relabeled': relabeled,
                    'mapped_total': len(mapping),
                })
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
