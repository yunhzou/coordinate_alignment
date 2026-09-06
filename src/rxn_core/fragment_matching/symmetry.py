"""Exact correlated actions projected to fragment occupations, not bijections."""
from dataclasses import replace

from ..alignment.post_aam import AAMHierarchy
from ..matcher import _nauty_atom_generators
from ..subgraph import _coerce_graph
from .._group_ops import occupation_orbit, OccupationLimitExceeded


class FragmentOrbitLimitExceeded(RuntimeError):
    def __init__(self, count, limit):
        self.count, self.limit = int(count), int(limit)
        super().__init__(f"fragment occupation limit: {count}>{limit}")


def materialize_target_coverage_orbit(candidate, target, *, iso_tolerance=0.5,
                                      limit=None, generators=None, observed_atoms=None):
    """One witness per correlated fragment-region relation.

    The orbit walk transports only integer images. Hierarchy/generator objects
    are transported once for each final region, never on every group edge.
    Proven equivalent source fragment units form an unordered multiset.
    With observed_atoms, choices permanently outside their generator-reachable
    closure stay compressed in the original hierarchy. Every returned mapping
    remains a complete, correlated witness, including competitor assignments.
    """
    graph = _coerce_graph(target, 0.2)
    generators = tuple(_nauty_atom_generators(graph, wbo_tol=iso_tolerance)
                       if generators is None else generators)
    observed_atoms = None if observed_atoms is None else tuple(observed_atoms)
    families = []
    for derivation in candidate.derivations:
        paths = derivation.initial_paths[:1] + derivation.residual_paths
        mapping = {a: p for path in paths for a, p in path.mapping.items()
                   if a in candidate.retained_atoms}
        action = dict(derivation.target_action)
        mapping = {a: action.get(p, p) for a, p in mapping.items()}
        if set(mapping) != set(candidate.retained_atoms):
            raise ValueError("derivation does not explain retained atoms")
        hierarchy = AAMHierarchy(tuple(f for path in paths for f in path.hierarchy.fragments))
        families.append(replace(candidate, mapping=tuple(sorted(mapping.items())),
            aam_hierarchy=hierarchy.relabel_target(action), derivations=(derivation,),
            covered_target_atoms=tuple(sorted(mapping.values())),
            attachment_atoms_target=tuple(sorted(mapping[a] for a in candidate.attachment_atoms_source))))
    if not candidate.derivations:
        families.append(candidate)  # Explicitly supplied relation, with no search provenance.

    final = {}
    for family in families:
        source_atoms = tuple(a for a, _p in family.mapping)
        positions = {a: i for i, a in enumerate(source_atoms)}
        classes = family.fragment_classes or tuple(range(len(family.retained_fragments)))
        fragment_positions = tuple((label, tuple(positions[a] for a in fragment))
            for label, fragment in zip(classes, family.retained_fragments, strict=True))
        attachments = tuple(positions[a] for a in family.attachment_atoms_source)
        bonds = tuple((positions[a], positions[b]) for a, b in family.preserved_source_bonds)
        local = (() if family.derivations and family.derivations[0].occupation_projected
                 else family.aam_hierarchy.fragments)
        stages = [tuple(dict(enumerate(g.images)) for g in (f.target_generators or ()))
                  for f in reversed(local)] + [generators]
        degree = max([max(graph.nodes(), default=-1) + 1,
                      max(family.covered_target_atoms, default=-1) + 1]
                     + [max(g, default=-1) + 1 for stage in stages for g in stage])
        witness = tuple(p for _a, p in family.mapping)
        def key(images):
            return (tuple(sorted(images)), tuple(sorted(images[i] for i in attachments)),
                    tuple(sorted((label, tuple(sorted(images[i] for i in part)))
                                 for label, part in fragment_positions)),
                    tuple(sorted(tuple(sorted((images[a], images[b]))) for a, b in bonds)))
        actions = tuple(tuple(tuple(g.get(a, a) for a in range(degree)) for g in stage)
                        for stage in stages)
        try:
            states = occupation_orbit(witness, degree, actions, attachments,
                                      fragment_positions, bonds, -1 if limit is None else limit,
                                      observed_atoms)
        except OccupationLimitExceeded:
            raise FragmentOrbitLimitExceeded(limit + 1, limit) from None
        for raw_images, raw_action in states:
            images, action = tuple(raw_images), tuple(raw_action)
            relation = key(images)
            final.setdefault(relation, (family, source_atoms, images, action))
    output = []
    for relation in sorted(final):
        family, atoms, images, action = final[relation]
        action_map = dict(enumerate(action))
        derivations = []
        for derivation in family.derivations:
            prior = dict(derivation.target_action)
            domain = set(prior) | set(action_map)
            combined = tuple(sorted((a, action_map.get(prior.get(a, a), prior.get(a, a)))
                                    for a in domain))
            derivations.append(replace(derivation, target_action=combined))
        mapping = dict(zip(atoms, images))
        output.append(replace(family, mapping=tuple(sorted(mapping.items())),
            covered_target_atoms=relation[0], attachment_atoms_target=relation[1],
            aam_hierarchy=family.aam_hierarchy.relabel_target(action_map),
            derivations=tuple(derivations)))
    return tuple(output)
