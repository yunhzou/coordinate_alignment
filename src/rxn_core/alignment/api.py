"""Public API for WBO graph matching."""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from pathlib import Path

from ..chemistry_computations import run_xtb
from ..frag import build_graph, expand_mapping, classify_bonds
from .branch import (
    _chirality_violations,
    _generate_seed_orders,
    find_islands,
    symmetry_repair_mapping,
)
from ..matcher import _nauty_orbits


@dataclass
class MatchCandidate:
    """One scored molecule-level alignment candidate.

    `mapping` is the symmetry-aware selected representative used for mechanism
    scoring. `raw_mapping` is the deterministic witness emitted by fragment
    growth before optional symmetry repair. `symmetry_fragments` records the
    compressed fragment states that produced the witness.
    """
    seed_index: int
    seed_atom: int
    branch_index: int
    cut_edges: tuple
    raw_mapping: dict
    mapping: dict
    broken: list
    formed: list
    chirality_violations: int
    score: tuple
    deferred_edges: tuple
    symmetry_fragments: list
    events: list | None = None
    search_path: object = None


@dataclass
class MatchResult:
    """Result returned by `match_wbo_graphs`."""
    candidates: list
    best: MatchCandidate | None
    graph_floor: float
    iso_tol: float
    n_seeds: int
    cut_edges: tuple


def _apply_cut_edges(g_R, cut_edges):
    if not cut_edges:
        return g_R
    g = g_R.copy()
    g.graph.update(g_R.graph)
    for i, j in cut_edges:
        if g.has_edge(i, j):
            g.remove_edge(i, j)
    return g


def cut_edges_above_floor(wboR, floor=0.2):
    """All R atom pairs whose WBO is at or above the growth floor."""
    n = int(wboR.shape[0])
    return tuple((i, j) for i in range(n) for j in range(i + 1, n)
                 if float(wboR[i, j]) >= floor)


def match_wbo_graphs(elR, wboR, elP, wboP, *,
                     xyzR=None, xyzP=None,
                     graph_floor=0.2, iso_tol=1.0,
                     dwbo_threshold=0.5, metal_dwbo_threshold=0.3,
                     symmetry_wbo_tol=0.2,
                     n_seeds=3, max_branches=1_000_000,
                     cut_edges=(), repair_symmetry=True,
                     chirality=False, capture_events=False, profile=None):
    """Symmetry-centric molecule match.

    This is the public low-level match function.  One match means:

    1. Build WBO graphs.
    2. Run exactly `n_seeds` seed orderings by default.
    3. For each seed, grow fragments with compressed symmetry candidates.
    4. Score candidates by least broken+formed bonds.

    `cut_edges` changes only the R growth graph.  Bond classification and
    symmetry repair still use the full WBO matrices.
    """
    if Counter(elR) != Counter(elP):
        raise ValueError(f"composition mismatch: {Counter(elR)} vs {Counter(elP)}")
    g_R_full = build_graph(elR, wboR, bond_cut=graph_floor)
    g_R = _apply_cut_edges(g_R_full, tuple(tuple(e) for e in cut_edges or ()))
    g_P = build_graph(elP, wboP, bond_cut=graph_floor)
    symmetry_wbo_tol = float(iso_tol)
    p_orbits = _nauty_orbits(g_P, wbo_tol=iso_tol)

    candidates = []
    orders = _generate_seed_orders(g_R, n_seeds)
    for seed_index, order in enumerate(orders):
        events = [] if capture_events else None
        search_graph = find_islands(
            g_R, g_P, order,
            graph_floor=graph_floor, iso_tol=iso_tol,
            dwbo_threshold=dwbo_threshold,
            metal_dwbo_threshold=metal_dwbo_threshold,
            symmetry_wbo_tol=symmetry_wbo_tol,
            max_branches=max_branches, events=events, profile=profile)
        for branch_index, branch in enumerate(search_graph.paths()):
            raw_mapping = expand_mapping(branch.mapping, g_R, g_P)
            mapping = dict(raw_mapping)
            if repair_symmetry:
                mapping = symmetry_repair_mapping(
                    mapping, wboR, wboP, g_R_full, g_P, p_orbits,
                    dwbo_threshold=dwbo_threshold,
                    metal_dwbo_threshold=metal_dwbo_threshold)
            broken, formed, _, _ = classify_bonds(
                mapping, wboR, wboP, dwbo_threshold=dwbo_threshold,
                elements_R=elR, elements_P=elP,
                metal_dwbo_threshold=metal_dwbo_threshold)
            chir = 0
            if chirality and xyzR is not None and xyzP is not None:
                chir = _chirality_violations(
                    mapping, xyzR, xyzP, broken, formed, elR)
            score = (len(broken) + len(formed), chir)
            candidates.append(MatchCandidate(
                seed_index=seed_index,
                seed_atom=int(order[0]) if order else -1,
                branch_index=branch_index,
                cut_edges=tuple(tuple(map(int, e)) for e in cut_edges or ()),
                raw_mapping=raw_mapping,
                mapping=mapping,
                broken=broken,
                formed=formed,
                chirality_violations=chir,
                score=score,
                deferred_edges=tuple(sorted(tuple(map(int, e))
                                            for e in branch.deferred_edges)),
                symmetry_fragments=list(branch.fragments),
                events=events if capture_events and branch_index == 0 else None,
                search_path=branch,
            ))
    candidates.sort(key=lambda c: c.score)
    best = candidates[0] if candidates else None
    return MatchResult(
        candidates=candidates, best=best,
        graph_floor=graph_floor, iso_tol=iso_tol, n_seeds=n_seeds,
        cut_edges=tuple(tuple(map(int, e)) for e in cut_edges or ()))


def align_from_arrays(elR, xyzR, wboR, elP, xyzP, wboP,
                      graph_floor=0.2, iso_tol=1,
                      dwbo_threshold=0.5, metal_dwbo_threshold=0.3,
                      symmetry_wbo_tol=0.2,
                      n_seeds=3, max_branches=1_000_000,
                      chirality=True, return_all=False):
    """Pure-graph entry point: assumes (el, xyz, wbo) for R and P are
    already in hand (e.g. loaded from a pre-computed xtb cache). Runs the
    multi-seed symmetry-aware search + scoring."""
    result = match_wbo_graphs(
        elR, wboR, elP, wboP,
        xyzR=xyzR, xyzP=xyzP,
        graph_floor=graph_floor, iso_tol=iso_tol,
        dwbo_threshold=dwbo_threshold,
        metal_dwbo_threshold=metal_dwbo_threshold,
        symmetry_wbo_tol=symmetry_wbo_tol,
        n_seeds=n_seeds, max_branches=max_branches,
        repair_symmetry=True, chirality=chirality)
    if result.best is None:
        raise RuntimeError("no alignment candidates produced")
    best = result.best
    mapping = best.mapping
    broken = best.broken
    formed = best.formed
    chir = best.chirality_violations
    out = dict(
        elements_R=elR, coords_R=xyzR, wbo_R=wboR,
        elements_P=elP, coords_P=xyzP, wbo_P=wboP,
        mapping=mapping, broken=broken, formed=formed,
        n_mapped=len(mapping), n_broken=len(broken), n_formed=len(formed),
        chirality_violations=chir,
        score=best.score,
    )
    if return_all:
        out['all_scored'] = [
            (c.score, c.mapping, c.broken, c.formed, c.chirality_violations)
            for c in result.candidates
        ]
    return out


def analyze_alignment(reactant_xyz, product_xyz, workdir,
                      charge=0, multiplicity=1,
                      graph_floor=0.2, iso_tol=1.0,
                      dwbo_threshold=0.5, metal_dwbo_threshold=0.3,
                      symmetry_wbo_tol=0.2,
                      n_seeds=3, max_branches=1_000_000,
                      chirality=True,
                      return_all=False):
    """Run xtb on R and P, build graphs, align, and score by mechanism."""
    if int(multiplicity) < 1:
        raise ValueError("multiplicity must be >= 1")
    xtb_uhf = int(multiplicity) - 1
    workdir = Path(workdir)
    elR, xyzR, wboR = run_xtb(
        reactant_xyz, workdir / "R", charge=charge, uhf=xtb_uhf)
    elP, xyzP, wboP = run_xtb(
        product_xyz, workdir / "P", charge=charge, uhf=xtb_uhf)
    return align_from_arrays(
        elR, xyzR, wboR, elP, xyzP, wboP,
        graph_floor=graph_floor, iso_tol=iso_tol,
        dwbo_threshold=dwbo_threshold,
        metal_dwbo_threshold=metal_dwbo_threshold,
        symmetry_wbo_tol=symmetry_wbo_tol,
        n_seeds=n_seeds, max_branches=max_branches,
        chirality=chirality, return_all=return_all)
