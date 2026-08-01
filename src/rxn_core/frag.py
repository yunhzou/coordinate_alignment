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
