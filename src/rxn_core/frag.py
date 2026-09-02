"""
Low-level utilities used by rxn_core:

  build_graph               — WBO graph (edge iff WBO >= bond_cut)
  expand_mapping            — element-multiset pairing of unmapped neighbors
  classify_bonds            — broken/formed bond classification by ΔWBO

Alignment logic lives in the `alignment`, `growth`, and `matcher` packages.
XYZ and xtb helpers live in `chemistry_computations` and are re-exported
here for compatibility with existing imports.
"""
from __future__ import annotations

import functools
from collections import defaultdict
from collections.abc import Mapping
from dataclasses import dataclass, field
from typing import Any

import numpy as np
import networkx as nx

from .chemistry_computations import parse_xyz, run_xtb, write_xyz_str


METAL_ELEMENTS = frozenset({
    "Li", "Be", "Na", "Mg", "Al", "K", "Ca", "Sc", "Ti", "V", "Cr",
    "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Rb", "Sr", "Y", "Zr",
    "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Cs",
    "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy",
    "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir",
    "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "Fr", "Ra", "Ac", "Th",
    "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md",
    "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn",
    "Nh", "Fl", "Mc", "Lv",
})


def is_metal_element(element):
    """Return True for elements treated with the metal WBO-event cutoff."""
    if element is None:
        return False
    s = str(element).strip()
    if not s:
        return False
    return s[0].upper() + s[1:].lower() in METAL_ELEMENTS


def bond_event_threshold(elements, i, j, *,
                         default_threshold=0.5,
                         metal_threshold=None):
    """Pair-specific delta-WBO threshold for mechanism events.

    The default threshold is used for ordinary covalent/organic pairs.  If a
    metal threshold is supplied and either endpoint is a metal, the lower
    metal-aware cutoff is used.
    """
    if elements is None or metal_threshold is None:
        return float(default_threshold)
    if is_metal_element(elements[i]) or is_metal_element(elements[j]):
        return float(metal_threshold)
    return float(default_threshold)


# -------------------- WBO graph --------------------


@dataclass
class WeightedNode:
    """Node record for generalized weighted-graph matching."""

    element: str | None = None
    label: str | None = None
    features: dict[str, Any] = field(default_factory=dict)
    attrs: dict[str, Any] = field(default_factory=dict)

    def as_attrs(self):
        out = dict(self.attrs)
        if self.element is not None:
            out["element"] = self.element
        if self.label is not None:
            out["label"] = self.label
        out["features"] = dict(self.features)
        for key, value in self.features.items():
            out.setdefault(key, value)
        return out


@dataclass
class WeightedGraph:
    """Typed weighted graph used by generalized subgraph matching."""

    nodes: list[WeightedNode | Mapping[str, Any] | str]
    weights: np.ndarray
    weight_name: str = "wbo"
    coords: np.ndarray | None = None
    metadata: dict[str, Any] = field(default_factory=dict)

    def to_networkx(self, bond_cut=0.5):
        return build_weighted_graph(
            self.nodes, self.weights, bond_cut=bond_cut,
            weight_name=self.weight_name, coords=self.coords,
            metadata=self.metadata)


def _node_attrs(node):
    if isinstance(node, WeightedNode):
        return node.as_attrs()
    if isinstance(node, str):
        return {"element": node, "features": {}}
    if isinstance(node, Mapping):
        out = dict(node)
        features = dict(out.get("features") or {})
        out["features"] = features
        for key, value in features.items():
            out.setdefault(key, value)
        return out
    raise TypeError(f"unsupported weighted graph node: {node!r}")


def build_weighted_graph(nodes, weights, bond_cut=0.5, *,
                         weight_name="wbo", coords=None, metadata=None):
    """Build a NetworkX graph from arbitrary node descriptors and weights."""
    weights = np.asarray(weights, dtype=float)
    if weights.ndim != 2 or weights.shape[0] != weights.shape[1]:
        raise ValueError("weights must be a square matrix")
    if weights.shape[0] != len(nodes):
        raise ValueError("number of nodes must match weights shape")
    g = nx.Graph()
    g.graph["wbo_matrix"] = weights
    g.graph["bond_cut"] = float(bond_cut)
    g.graph["weight_name"] = str(weight_name)
    if coords is not None:
        g.graph["coords"] = np.asarray(coords, dtype=float)
    g.graph.update(dict(metadata or {}))
    for i, node in enumerate(nodes):
        g.add_node(i, **_node_attrs(node))
    n = len(nodes)
    for i in range(n):
        for j in range(i + 1, n):
            if weights[i, j] >= bond_cut:
                w = float(weights[i, j])
                g.add_edge(i, j, **{weight_name: w, "wbo": w})
    return g


def build_graph(elements, wbo, bond_cut=0.5):
    """Connectivity graph with element on each node and WBO weight on each
    edge. Bond exists iff WBO >= bond_cut; the full WBO matrix is retained on
    the graph so scoring and trace diagnostics can inspect exact WBO values."""
    g = nx.Graph()
    g.graph["wbo_matrix"] = np.asarray(wbo, dtype=float)
    g.graph["bond_cut"] = float(bond_cut)
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


def _classify_bonds_reference(mapping, wbo_R, wbo_P, dwbo_threshold=0.5,
                              elements_R=None, elements_P=None,
                              metal_dwbo_threshold=None):
    """Scalar-loop reference implementation of :func:`classify_bonds`.

    Kept verbatim as the exactness oracle for the vectorised version and as
    the fallback for inputs outside its fast path.

    Bond classification by WBO change.

        broken iff (WBO_R - WBO_P) >= pair_threshold
        formed iff (WBO_P - WBO_R) >= pair_threshold

    `pair_threshold` is normally `dwbo_threshold`.  If element lists and
    `metal_dwbo_threshold` are supplied, any pair containing a metal uses the
    metal threshold instead.  Since WBO is non-negative, requiring
    (wR - wP) >= pair_threshold automatically implies wR >= pair_threshold
    (wP >= 0), so a single
    threshold both gates "is this a bond worth considering" and "did
    its order change enough." For pairs with one or both endpoints
    unmapped (no image bond defined) we treat the missing wP as 0 and
    apply the same threshold to wR.

    Returns (broken_list, formed_list, core_R, core_P) where each bond
    record is (i, j, wbo_R_or_None, wbo_P_or_None)."""
    inv = {v: k for k, v in mapping.items()}
    nR, nP = wbo_R.shape[0], wbo_P.shape[0]
    default_threshold = float(dwbo_threshold)
    metal_threshold = (
        None if metal_dwbo_threshold is None
        else float(metal_dwbo_threshold))
    metal_R = (
        None if elements_R is None or metal_threshold is None
        else tuple(is_metal_element(element) for element in elements_R))
    metal_P = (
        None if elements_P is None or metal_threshold is None
        else tuple(is_metal_element(element) for element in elements_P))
    broken, formed = [], []
    for i in range(nR):
        for j in range(i + 1, nR):
            threshold = (
                metal_threshold
                if metal_R is not None and (metal_R[i] or metal_R[j])
                else default_threshold)
            wR = wbo_R[i, j]
            if wR < threshold: continue
            if i not in mapping or j not in mapping:
                broken.append((i, j, float(wR), None)); continue
            wP = wbo_P[mapping[i], mapping[j]]
            if wR - wP >= threshold:
                broken.append((i, j, float(wR), float(wP)))
    for ip in range(nP):
        for jp in range(ip + 1, nP):
            threshold = (
                metal_threshold
                if metal_P is not None and (metal_P[ip] or metal_P[jp])
                else default_threshold)
            wP = wbo_P[ip, jp]
            if wP < threshold: continue
            if ip not in inv or jp not in inv:
                formed.append((ip, jp, None, float(wP))); continue
            wR = wbo_R[inv[ip], inv[jp]]
            if wP - wR >= threshold:
                formed.append((ip, jp, float(wR), float(wP)))
    core_R = sorted({i for (i, j, _, _) in broken}
                    | {j for (i, j, _, _) in broken})
    core_P = sorted({i for (i, j, _, _) in formed}
                    | {j for (i, j, _, _) in formed})
    return broken, formed, core_R, core_P


@functools.lru_cache(maxsize=64)
def _upper_pair_indices(n):
    """Row-major ``i < j`` index pairs, i.e. the nested-loop order."""
    return np.triu_indices(int(n), 1)


def _integer_mapping_items(mapping):
    """Return ``[(int key, int value), ...]`` or None if any is not integral."""
    items = []
    for key, value in mapping.items():
        if (not isinstance(key, (int, np.integer))
                or not isinstance(value, (int, np.integer))):
            return None
        items.append((int(key), int(value)))
    return items


def _classify_bond_side(w_self, w_other, mapped, image, metal,
                        default_threshold, metal_threshold):
    """One side of :func:`classify_bonds`, vectorised.

    Emits ``(i, j, w_self, w_other_or_None)`` for pairs ``i < j`` of ``w_self``
    in row-major order.  Every arithmetic step mirrors the reference loop on
    float64 data: the same per-pair threshold (metal if either endpoint is a
    metal), the same ``not (w < threshold)`` gate (so NaN passes through as
    before), ``w_other`` looked up through the same integer images (negative
    images wrap exactly as scalar indexing did), and the same
    ``w_self - w_other >= threshold`` test.  Values reach Python via
    ``ndarray.tolist()``, which yields the same floats as ``float(scalar)``.
    """
    n = w_self.shape[0]
    if n < 2:
        return []
    iu, ju = _upper_pair_indices(n)
    w = w_self[iu, ju]
    if metal is None:
        threshold = np.full(w.shape, default_threshold, dtype=float)
    else:
        threshold = np.where(
            metal[iu] | metal[ju], metal_threshold, default_threshold)
    keep = ~(w < threshold)
    if not keep.any():
        return []
    i_s = iu[keep]
    j_s = ju[keep]
    w_s = w[keep]
    t_s = threshold[keep]
    both = mapped[i_s] & mapped[j_s]
    emit = ~both
    other = np.full(w_s.shape, np.nan, dtype=float)
    idx = np.flatnonzero(both)
    if idx.size:
        looked_up = w_other[image[i_s[idx]], image[j_s[idx]]]
        other[idx] = looked_up
        emit[idx[(w_s[idx] - looked_up) >= t_s[idx]]] = True
    out = []
    for i, j, ws, wo, is_both, do_emit in zip(
            i_s.tolist(), j_s.tolist(), w_s.tolist(), other.tolist(),
            both.tolist(), emit.tolist()):
        if do_emit:
            out.append((i, j, ws, wo if is_both else None))
    return out


def classify_bonds(mapping, wbo_R, wbo_P, dwbo_threshold=0.5,
                   elements_R=None, elements_P=None,
                   metal_dwbo_threshold=None):
    """Bond classification by WBO change.

        broken iff (WBO_R - WBO_P) >= pair_threshold
        formed iff (WBO_P - WBO_R) >= pair_threshold

    `pair_threshold` is normally `dwbo_threshold`.  If element lists and
    `metal_dwbo_threshold` are supplied, any pair containing a metal uses the
    metal threshold instead.  Since WBO is non-negative, requiring
    (wR - wP) >= pair_threshold automatically implies wR >= pair_threshold
    (wP >= 0), so a single
    threshold both gates "is this a bond worth considering" and "did
    its order change enough." For pairs with one or both endpoints
    unmapped (no image bond defined) we treat the missing wP as 0 and
    apply the same threshold to wR.

    Returns (broken_list, formed_list, core_R, core_P) where each bond
    record is (i, j, wbo_R_or_None, wbo_P_or_None).

    Vectorised over all ``i < j`` pairs; identical output to
    :func:`_classify_bonds_reference` (same record order, floats, ``None``
    handling and threshold selection).  Inputs outside the fast path (non
    float64 base-class ndarrays, non-integer mapping entries) are delegated to
    the reference implementation unchanged.
    """
    if (type(wbo_R) is not np.ndarray or type(wbo_P) is not np.ndarray
            or wbo_R.dtype != np.float64 or wbo_P.dtype != np.float64
            or wbo_R.ndim != 2 or wbo_P.ndim != 2):
        return _classify_bonds_reference(
            mapping, wbo_R, wbo_P, dwbo_threshold=dwbo_threshold,
            elements_R=elements_R, elements_P=elements_P,
            metal_dwbo_threshold=metal_dwbo_threshold)
    items = _integer_mapping_items(mapping)
    if items is None:
        return _classify_bonds_reference(
            mapping, wbo_R, wbo_P, dwbo_threshold=dwbo_threshold,
            elements_R=elements_R, elements_P=elements_P,
            metal_dwbo_threshold=metal_dwbo_threshold)
    nR, nP = wbo_R.shape[0], wbo_P.shape[0]
    default_threshold = float(dwbo_threshold)
    metal_threshold = (
        None if metal_dwbo_threshold is None
        else float(metal_dwbo_threshold))
    metal_R = (
        None if elements_R is None or metal_threshold is None
        else np.array([is_metal_element(e) for e in elements_R], dtype=bool))
    metal_P = (
        None if elements_P is None or metal_threshold is None
        else np.array([is_metal_element(e) for e in elements_P], dtype=bool))
    # ``i in mapping`` for i in range(nR): only in-range integer keys matter.
    # ``inv`` is the last-writer-wins inverse, reproduced by assignment order.
    mapped_R = np.zeros(nR, dtype=bool)
    image_R = np.zeros(nR, dtype=np.intp)
    mapped_P = np.zeros(nP, dtype=bool)
    image_P = np.zeros(nP, dtype=np.intp)
    for r, p in items:
        if 0 <= r < nR:
            mapped_R[r] = True
            image_R[r] = p
        if 0 <= p < nP:
            mapped_P[p] = True
            image_P[p] = r
    broken = _classify_bond_side(
        wbo_R, wbo_P, mapped_R, image_R, metal_R,
        default_threshold, metal_threshold)
    formed = [
        (ip, jp, wR, wP)
        for ip, jp, wP, wR in _classify_bond_side(
            wbo_P, wbo_R, mapped_P, image_P, metal_P,
            default_threshold, metal_threshold)
    ]
    core_R = sorted({i for (i, j, _, _) in broken}
                    | {j for (i, j, _, _) in broken})
    core_P = sorted({i for (i, j, _, _) in formed}
                    | {j for (i, j, _, _) in formed})
    return broken, formed, core_R, core_P
