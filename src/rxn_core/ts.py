"""Typed transition-state processing composed from partial AAM and scoring."""
from __future__ import annotations

import time
from collections import defaultdict

import numpy as np

from .core_aam import search_core_assignments
from .domain import (
    AAMSearchConfig,
    AtomAssignment,
    ReactionContext,
    ResolvedMechanism,
    RPResult,
    TSMechanismResult,
    TSResult,
    TSScore,
    TSScoringConfig,
    TransitionStateTarget,
)
from .modes import bond_overlap_per_mode, reindex_modes_to_R


def _clip01(value):
    return max(0.0, min(1.0, float(value)))


def _score_assignment(rp, mechanism, target, assignment, sources, config):
    reactant = rp.problem.reactant
    product = rp.problem.product
    target_molecule = target.molecule
    mapping_rp = mechanism.mapping.as_dict()
    mapping_rt = assignment.as_dict()
    terms = []
    for kind, bonds in (("broken", mechanism.broken_bonds),
                        ("formed", mechanism.formed_bonds)):
        for left, right in bonds:
            product_left, product_right = mapping_rp[left], mapping_rp[right]
            target_left, target_right = mapping_rt[left], mapping_rt[right]
            w_r = float(reactant.wbo[left, right])
            w_p = float(product.wbo[product_left, product_right])
            w_t = float(target_molecule.wbo[target_left, target_right])
            delta = abs(w_p - w_r)
            if delta < 1e-12:
                progress = 1.0
            elif kind == "formed":
                progress = _clip01((w_t - w_r) / delta)
            else:
                progress = _clip01((w_r - w_t) / delta)
            terms.append({
                "kind": kind,
                "reactant_bond": (left, right),
                "product_bond": (product_left, product_right),
                "target_bond": (target_left, target_right),
                "wbo_reactant": w_r,
                "wbo_product": w_p,
                "wbo_target": w_t,
                "event_weight": delta ** config.event_weight_power,
                "target_progress": progress,
            })

    total_weight = sum(term["event_weight"] for term in terms)
    progress = (1.0 if total_weight < 1e-12 else sum(
        term["event_weight"] * term["target_progress"] for term in terms
    ) / total_weight)
    target_in_r = np.array(reactant.coordinates, copy=True)
    for source, image in assignment.pairs:
        target_in_r[source] = target_molecule.coordinates[image]
    reaction_vector = np.zeros_like(target_in_r)
    for term in terms:
        left, right = term["reactant_bond"]
        vector = target_in_r[right] - target_in_r[left]
        norm = float(np.linalg.norm(vector))
        if norm < 1e-9:
            continue
        unit = term["event_weight"] * vector / norm
        if term["kind"] == "broken":
            reaction_vector[left] -= unit
            reaction_vector[right] += unit
        else:
            reaction_vector[left] += unit
            reaction_vector[right] -= unit

    modes = target.vibrations.displacements
    modes_r = reindex_modes_to_R(modes, mapping_rt, reactant.atom_count)
    mode_norms = np.linalg.norm(modes.reshape(modes.shape[0], -1), axis=1)
    overlaps = bond_overlap_per_mode(
        modes_r, reaction_vector, mode_norms=mode_norms)
    imaginary = np.flatnonzero(target.vibrations.frequencies < 0)
    if not len(imaginary):
        return None
    mode_index = int(max(imaginary, key=lambda index: overlaps[index]))
    overlap = float(overlaps[mode_index])
    return TSScore(
        assignment=assignment,
        sources=frozenset(sources),
        score=overlap * float(progress) ** config.wbo_progress_power,
        overlap=overlap,
        wbo_progress=float(progress),
        mode_index=mode_index,
        frequency=float(target.vibrations.frequencies[mode_index]),
        event_terms=tuple(terms),
    )


def reaction_context_from_rp(rp: RPResult) -> ReactionContext:
    """Project a full analytical R/P result onto the exact TS dependency."""
    if not isinstance(rp, RPResult):
        raise TypeError("reaction_context_from_rp requires an RPResult")
    return ReactionContext(
        problem=rp.analytical.aam.problem,
        config=rp.analytical.aam.config,
        mechanisms=tuple(ResolvedMechanism(
            mapping=item.mapping,
            broken_bonds=item.broken_bonds,
            formed_bonds=item.formed_bonds,
            core_atoms=item.core_atoms,
        ) for item in rp.mechanisms),
    )


def analyze_transition_state(
        reaction: RPResult | ReactionContext, target: TransitionStateTarget, *,
        search_config: AAMSearchConfig | None = None,
        scoring_config: TSScoringConfig | None = None) -> TSResult:
    """Evaluate one TS target under every selected R/P mechanism.

    R->TS and P->TS are independent partial-AAM searches.  P assignments are
    pulled into the R frame through the already selected R/P bijection, then
    exact tuples are unioned with source provenance before mode scoring.
    """
    if isinstance(reaction, RPResult):
        reaction = reaction_context_from_rp(reaction)
    if (not isinstance(reaction, ReactionContext)
            or not isinstance(target, TransitionStateTarget)):
        raise TypeError(
            "typed RPResult/ReactionContext and TransitionStateTarget are required")
    search_config = search_config or reaction.config
    scoring_config = scoring_config or TSScoringConfig()
    started = time.perf_counter()
    results = []
    problem = reaction.problem
    for mechanism in reaction.mechanisms:
        core_r = mechanism.core_atoms
        if not core_r:
            results.append(TSMechanismResult(
                mechanism=mechanism,
                target=target,
                reactant_core_aam=None,
                product_core_aam=None,
                candidates=(),
                selected=None,
                status="no_reaction_core",
                reason=(
                    "the selected R/P mechanism has no broken or formed "
                    "bonds, so reactive-mode scoring is undefined"),
            ))
            continue
        core_p = tuple(mechanism.mapping.images[atom] for atom in core_r)
        from_r = search_core_assignments(
            problem.reactant, target.molecule, core_r,
            config=search_config,
            assignment_limit=scoring_config.core_assignment_limit)
        from_p = search_core_assignments(
            problem.product, target.molecule, core_p,
            config=search_config,
            assignment_limit=scoring_config.core_assignment_limit)
        inverse = mechanism.mapping.inverse().as_dict()
        candidates = defaultdict(set)
        objects = {}
        for assignment in from_r.assignments:
            candidates[assignment.pairs].add("reactant")
            objects[assignment.pairs] = assignment
        for assignment in from_p.assignments:
            pulled = AtomAssignment(tuple(
                (inverse[source_p], target_atom)
                for source_p, target_atom in assignment.pairs))
            candidates[pulled.pairs].add("product")
            objects[pulled.pairs] = pulled
        scores = tuple(filter(None, (
            _score_assignment(
                reaction, mechanism, target, objects[key], candidates[key],
                scoring_config)
            for key in sorted(candidates)
        )))
        eligible = scores
        if scoring_config.prefer_endpoint_consensus:
            consensus = tuple(score for score in scores
                              if len(score.sources) == 2)
            if consensus:
                eligible = consensus
        selected = max(
            eligible,
            key=lambda score: (score.score, -score.mode_index,
                               tuple(-image for _, image in score.assignment.pairs)),
            default=None)
        results.append(TSMechanismResult(
            mechanism=mechanism,
            target=target,
            reactant_core_aam=from_r,
            product_core_aam=from_p,
            candidates=scores,
            selected=selected,
            status="scored",
        ))
    return TSResult(
        reaction=reaction,
        mechanisms=tuple(results),
        elapsed_seconds=time.perf_counter() - started,
    )
