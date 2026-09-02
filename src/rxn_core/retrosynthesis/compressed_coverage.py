"""Recommendation over compressed target-domain claims."""
from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass

TargetDomain = tuple[int, ...]
CoverageSignature = tuple[TargetDomain, ...]
CoveragePattern = tuple[CoverageSignature, ...]

def candidate_target_domains(candidate):
    """Return one allowed target pool for every retained source atom.

    Ordinary AAM symmetry blocks are independent set-to-set assignment
    domains.  Exact automorphism blocks describe a coupled group action and
    therefore remain represented by their deterministic witness here.
    """
    retained = set(map(int, candidate.retained_atoms))
    witness = {int(source): int(target) for source, target in candidate.mapping}
    pools = {source: (target,) for source, target in witness.items()
             if source in retained}
    for fragment in candidate.aam_hierarchy.fragments:
        for domain in fragment.symmetry_domains:
            if domain.source == "exact_automorph_group":
                continue
            target_pool = tuple(map(int, domain.p_atoms))
            for source in domain.r_atoms:
                if source in retained:
                    pools[int(source)] = target_pool
    return tuple((source, pools[source]) for source in sorted(pools))


def coverage_signature(target_domains):
    """Source-index-free multiset signature used to group coverage modules."""
    return tuple(sorted(tuple(pool) for _source, pool in target_domains))


def assign_domain_signatures(signatures, atom_count):
    """Find one injective target witness without enumerating bijections."""
    domains = tuple(
        tuple(map(int, pool))
        for signature in signatures for pool in signature)
    assignment = _extend_domain_assignment((), (), domains)
    if assignment is None or len(domains) > atom_count:
        return None
    by_component = []
    cursor = 0
    for signature in signatures:
        width = len(signature)
        by_component.append(assignment[cursor:cursor + width])
        cursor += width
    return tuple(by_component)


def _extend_domain_assignment(domains, assignment, appended_domains):
    """Extend one maximum matching through augmenting paths."""
    domains = tuple(domains) + tuple(appended_domains)
    prior_count = len(assignment)
    target_owner = {}
    assigned = list(assignment) + [None] * len(appended_domains)
    for variable, target in enumerate(assignment):
        target_owner[target] = variable

    def augment(variable, visited):
        for target in domains[variable]:
            if target in visited:
                continue
            visited.add(target)
            prior = target_owner.get(target)
            if prior is None or augment(prior, visited):
                target_owner[target] = variable
                assigned[variable] = target
                return True
        return False

    new_variables = range(prior_count, len(domains))
    new_variables = sorted(
        new_variables,
        key=lambda variable: (len(domains[variable]), domains[variable]),
    )
    for variable in new_variables:
        if not augment(variable, set()):
            return None
    return tuple(assigned)


def _extend_with_required_target(
        domains, appended_domains, required_target):
    """Return one witness where the appended component owns the pivot."""
    combined = tuple(domains) + tuple(appended_domains)
    appended_start = len(domains)
    forced_variables = [
        appended_start + offset
        for offset, pool in enumerate(appended_domains)
        if required_target in pool
    ]
    for forced in forced_variables:
        assigned = [None] * len(combined)
        assigned[forced] = required_target
        target_owner = {required_target: forced}

        def augment(variable, visited):
            for target in combined[variable]:
                if target in visited:
                    continue
                visited.add(target)
                prior = target_owner.get(target)
                if prior == forced:
                    continue
                if prior is None or augment(prior, visited):
                    target_owner[target] = variable
                    assigned[variable] = target
                    return True
            return False

        variables = sorted(
            (variable for variable in range(len(combined))
             if variable != forced),
            key=lambda variable: (
                len(combined[variable]), combined[variable]),
        )
        if all(augment(variable, set()) for variable in variables):
            return tuple(assigned)
    return None


def assign_candidate_items(items, atom_count):
    """Materialize one final mapping witness for compressed candidate items."""
    signatures = [
        tuple(tuple(pool) for _source, pool in item["target_domains"])
        for item in items
    ]
    assignments = assign_domain_signatures(signatures, atom_count)
    if assignments is None:
        return None
    placed = []
    for item, targets in zip(items, assignments):
        transformed = dict(item)
        mapping = [
            [int(source), int(target)]
            for (source, _pool), target in zip(
                item["target_domains"], targets)
        ]
        target_by_source = dict(mapping)
        transformed.update({
            "mapping": mapping,
            "covered_target_atoms": sorted(target_by_source.values()),
            "attachment_atoms_target": sorted(
                target_by_source[source]
                for source in item["attachment_atoms_source"]
                if source in target_by_source),
        })
        placed.append(transformed)
    return tuple(placed)


@dataclass(frozen=True)
class CoverageRecommendationResult:
    patterns: tuple[CoveragePattern, ...]
    truncated: bool


@dataclass(frozen=True)
class CoverageRecommendationConfig:
    maximum_precursors: int = 3
    frontier_limit: int = 200

    def __post_init__(self):
        if min(self.maximum_precursors, self.frontier_limit) < 1:
            raise ValueError("recommendation limits must be positive")


def recommend_compressed_coverage_patterns(
        signatures, atom_count, rank_pattern, *, result_limit,
        config: CoverageRecommendationConfig | None = None):
    """Recommend full covers without materializing automorphic bijections.

    Each state contains compressed target domains plus one matching witness.
    A constrained uncovered target selects the next family.  Equivalent
    family multisets are merged immediately, and search stops after a bounded
    top-k frontier rather than enumerating every complete cover.
    """
    config = config or CoverageRecommendationConfig()
    signatures = tuple(sorted(set(signatures), key=lambda signature: (
        rank_pattern((signature,), len(signature)), signature)))
    by_target = defaultdict(list)
    for signature in signatures:
        if len(signature) > atom_count:
            continue
        for target in set(target for pool in signature for target in pool):
            by_target[target].append(signature)

    states = [((), (), ())]
    completed = {}
    truncated = False
    for _depth in range(config.maximum_precursors):
        next_by_multiset = {}
        for selected, domains, assignment in states:
            missing = set(range(atom_count)) - set(assignment)
            if not missing:
                completed.setdefault(tuple(sorted(selected)), selected)
                continue
            remaining_atoms = atom_count - len(domains)
            compatible_by_target = {
                target: [signature for signature in by_target[target]
                         if len(signature) <= remaining_atoms]
                for target in missing
            }
            pivot = min(
                missing,
                key=lambda target: (len(compatible_by_target[target]), target),
            )
            for signature in compatible_by_target[pivot]:
                appended = tuple(signature)
                new_assignment = _extend_with_required_target(
                    domains, appended, pivot)
                if new_assignment is None:
                    continue
                new_selected = selected + (signature,)
                new_domains = domains + appended
                if len(new_domains) == atom_count:
                    completed.setdefault(
                        tuple(sorted(new_selected)), new_selected)
                    continue
                key = tuple(sorted(new_selected))
                state = (new_selected, new_domains, new_assignment)
                prior = next_by_multiset.get(key)
                if (prior is None or rank_pattern(
                        new_selected, len(new_domains)) < rank_pattern(
                            prior[0], len(prior[1]))):
                    next_by_multiset[key] = state
        next_states = sorted(
            next_by_multiset.values(),
            key=lambda state: rank_pattern(state[0], len(state[1])),
        )
        if len(next_states) > config.frontier_limit:
            truncated = True
        states = next_states[:config.frontier_limit]
        if not states:
            break
    ranked = sorted(
        completed.values(),
        key=lambda pattern: rank_pattern(pattern, atom_count),
    )
    if len(ranked) > result_limit:
        truncated = True
    return CoverageRecommendationResult(
        patterns=tuple(ranked[:result_limit]),
        truncated=truncated,
    )
