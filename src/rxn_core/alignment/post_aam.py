"""Typed symmetry objects shared by AAM and optional mechanism post-processing."""
from __future__ import annotations

from dataclasses import dataclass, field
from functools import cached_property
from typing import Mapping, Protocol, Sequence

import numpy as np


@dataclass(frozen=True)
class AtomPermutation:
    """One exact permutation written as ``atom -> image``."""

    images: tuple[int, ...]

    def __post_init__(self):
        images = tuple(int(value) for value in self.images)
        if set(images) != set(range(len(images))):
            raise ValueError("atom permutation is not bijective")
        object.__setattr__(self, "images", images)

    @classmethod
    def identity(cls, degree: int):
        return cls(tuple(range(int(degree))))

    @classmethod
    def from_mapping(cls, mapping: Mapping[int, int], degree: int | None = None):
        raw = {int(atom): int(image) for atom, image in mapping.items()}
        degree = len(raw) if degree is None else int(degree)
        if set(raw) != set(range(degree)):
            raise ValueError("permutation mapping has an incomplete domain")
        return cls(tuple(raw[atom] for atom in range(degree)))

    @property
    def degree(self):
        return len(self.images)

    @property
    def is_identity(self):
        return all(atom == image for atom, image in enumerate(self.images))

    def inverse(self):
        inverse = [0] * self.degree
        for atom, image in enumerate(self.images):
            inverse[image] = atom
        return AtomPermutation(tuple(inverse))

    def then(self, other: "AtomPermutation"):
        """Return ``other(self(atom))``."""
        if self.degree != other.degree:
            raise ValueError("cannot compose permutations of different degree")
        return AtomPermutation(tuple(other.images[image] for image in self.images))

    def cycles(self):
        cycles, unseen = [], set(range(self.degree))
        while unseen:
            start = min(unseen)
            cycle, atom = [], start
            while atom in unseen:
                unseen.remove(atom)
                cycle.append(atom)
                atom = self.images[atom]
            if len(cycle) > 1:
                cycles.append(tuple(cycle))
        return tuple(cycles)


@dataclass(frozen=True)
class AtomBijection:
    """One complete R-to-P assignment."""

    images: tuple[int, ...]

    def __post_init__(self):
        images = tuple(int(value) for value in self.images)
        if set(images) != set(range(len(images))):
            raise ValueError("atom mapping is not a complete bijection")
        object.__setattr__(self, "images", images)

    @classmethod
    def from_mapping(cls, mapping: Mapping[int, int], degree: int | None = None):
        raw = {int(atom): int(image) for atom, image in mapping.items()}
        degree = len(raw) if degree is None else int(degree)
        if set(raw) != set(range(degree)):
            raise ValueError("atom mapping has an incomplete R domain")
        return cls(tuple(raw[atom] for atom in range(degree)))

    @property
    def degree(self):
        return len(self.images)

    def as_dict(self):
        return {atom: image for atom, image in enumerate(self.images)}

    def inverse(self):
        return AtomBijection(AtomPermutation(self.images).inverse().images)

    def act(self, source: AtomPermutation | None = None,
            target: AtomPermutation | None = None):
        """Return ``target o self o source^-1``.

        Source and target actions are explicit group elements.  No branch
        witness or correspondence search participates in this operation.
        """
        source = source or AtomPermutation.identity(self.degree)
        target = target or AtomPermutation.identity(self.degree)
        if source.degree != self.degree or target.degree != self.degree:
            raise ValueError("group action degree differs from atom mapping")
        source_inverse = source.inverse()
        return AtomBijection(tuple(
            target.images[self.images[source_inverse.images[atom]]]
            for atom in range(self.degree)))

    def product_in_reactant_order(self, product_coords):
        product = np.asarray(product_coords, dtype=float)
        if product.shape != (self.degree, 3):
            raise ValueError("product coordinates disagree with mapping degree")
        return product[np.asarray(self.images, dtype=int)]


@dataclass(frozen=True)
class PermutationGroup:
    """Generator representation of a finite atom-permutation group."""

    degree: int
    generators: tuple[AtomPermutation, ...] = ()

    def __post_init__(self):
        degree = int(self.degree)
        unique = {}
        for generator in self.generators:
            if generator.degree != degree:
                raise ValueError("group generator has the wrong degree")
            if not generator.is_identity:
                unique.setdefault(generator.images, generator)
        object.__setattr__(self, "degree", degree)
        object.__setattr__(self, "generators", tuple(unique.values()))

    @classmethod
    def trivial(cls, degree: int):
        return cls(int(degree), ())

    @classmethod
    def from_generator_mappings(cls, degree: int, generators):
        permutations = []
        for generator in generators:
            if isinstance(generator, Mapping):
                permutations.append(
                    AtomPermutation.from_mapping(generator, degree))
            else:
                permutation = AtomPermutation(tuple(map(int, generator)))
                if permutation.degree != int(degree):
                    raise ValueError("group generator has the wrong degree")
                permutations.append(permutation)
        return cls(int(degree), tuple(permutations))

    def orbits(self):
        parent = list(range(self.degree))

        def find(atom):
            while parent[atom] != atom:
                parent[atom] = parent[parent[atom]]
                atom = parent[atom]
            return atom

        def union(left, right):
            left, right = find(left), find(right)
            if left != right:
                parent[right] = left

        for generator in self.generators:
            for atom, image in enumerate(generator.images):
                union(atom, image)
        groups = {}
        for atom in range(self.degree):
            groups.setdefault(find(atom), []).append(atom)
        return tuple(tuple(members) for members in groups.values())


@dataclass(frozen=True)
class SymmetryDomain:
    """One compressed assignment domain carried by an AAM fragment."""

    r_atoms: tuple[int, ...]
    p_atoms: tuple[int, ...]
    source: str
    extendable: bool = False

    def __post_init__(self):
        object.__setattr__(self, "r_atoms", tuple(sorted(map(int, self.r_atoms))))
        object.__setattr__(self, "p_atoms", tuple(sorted(map(int, self.p_atoms))))


@dataclass(frozen=True)
class FragmentMatch:
    fragment_index: int
    island_index: int
    r_atoms: tuple[int, ...]
    deferred_edges: tuple[tuple[int, int], ...] = ()
    symmetry_domains: tuple[SymmetryDomain, ...] = ()
    representative_assignments: tuple[tuple[int, int], ...] = ()
    exact_fixed: tuple[int, ...] = ()
    multiplicity: int = 1
    automorph_domains: tuple[SymmetryDomain, ...] = ()
    target_generators: tuple[AtomPermutation, ...] | None = None

    def __post_init__(self):
        object.__setattr__(self, "r_atoms", tuple(sorted(map(int, self.r_atoms))))
        object.__setattr__(self, "deferred_edges", tuple(sorted(
            tuple(sorted(map(int, edge))) for edge in self.deferred_edges)))
        object.__setattr__(self, "symmetry_domains", tuple(self.symmetry_domains))
        object.__setattr__(self, "representative_assignments", tuple(sorted(
            (int(source), int(target))
            for source, target in self.representative_assignments)))
        object.__setattr__(self, "exact_fixed", tuple(sorted(
            map(int, self.exact_fixed))))
        object.__setattr__(self, "multiplicity", int(self.multiplicity))
        object.__setattr__(self, "automorph_domains", tuple(
            self.automorph_domains))
        if self.target_generators is not None:
            object.__setattr__(self, "target_generators", tuple(
                self.target_generators))

    @property
    def has_exact_target_group(self):
        """Whether AAM supplied the exact group (including a trivial one)."""
        return self.target_generators is not None


@dataclass(frozen=True)
class AAMHierarchy:
    fragments: tuple[FragmentMatch, ...]

    def relabel_target(self, action):
        """Keep a correlated transform as a shared view until inspected."""
        action = tuple(sorted(dict(action).items()))
        return AAMHierarchyView(self, action) if action else self

    def _materialize_target(self, action):
        """Transport a hierarchy by one correlated target permutation.

        Unmentioned atoms (including augmentation copies) are fixed. Exact
        generators are conjugated, not copied into the old target frame.
        """
        from dataclasses import replace
        action = dict(action)
        extent = max(max(action, default=-1), max(action.values(), default=-1)) + 1
        generators = tuple(g for f in self.fragments for g in (f.target_generators or ()))
        if all(a == b for a, b in action.items()) and all(g.degree >= extent for g in generators):
            return self
        frames = {}
        transported = {}
        def image(atom):
            return int(action.get(atom, atom))
        def domain(item):
            return replace(item, p_atoms=tuple(image(a) for a in item.p_atoms))
        def generator(item):
            if item.images not in transported:
                degree = max(item.degree, extent)
                if degree not in frames:
                    images = tuple(action.get(a, a) for a in range(degree))
                    inverse = [0] * degree
                    for atom, target in enumerate(images):
                        inverse[target] = atom
                    frames[degree] = images, inverse
                images, inverse = frames[degree]
                padded = item.images + tuple(range(item.degree, degree))
                moved = tuple(images[padded[a]] for a in inverse)
                transported[item.images] = (item if moved == item.images else AtomPermutation(moved))
            return transported[item.images]
        return AAMHierarchy(tuple(replace(fragment,
            representative_assignments=tuple((a, image(b)) for a, b in
                                               fragment.representative_assignments),
            symmetry_domains=tuple(domain(d) for d in fragment.symmetry_domains),
            automorph_domains=tuple(domain(d) for d in fragment.automorph_domains),
            target_generators=(None if fragment.target_generators is None else
                               tuple(generator(g) for g in fragment.target_generators)))
            for fragment in self.fragments))

    @classmethod
    def from_record(cls, branch_symmetry):
        fragments = []
        for position, raw in enumerate(
                dict(branch_symmetry or {}).get("fragments") or ()):
            symmetry = raw.get("symmetry") or {}
            raw_generators = symmetry.get("automorph_generators")
            domains = []
            for block in symmetry.get("blocks") or ():
                domains.append(SymmetryDomain(
                    tuple(block.get("r_atoms") or ()),
                    tuple(block.get("p_atoms") or ()),
                    str(block.get("source") or "sym_block"),
                    bool(block.get("extendable", False))))
            automorph_domains = tuple(SymmetryDomain(
                tuple(block.get("r_atoms") or ()),
                tuple(block.get("p_atoms") or ()),
                str(block.get("source") or "exact_automorph_group"),
                bool(block.get("extendable", False)))
                for block in symmetry.get("automorph_blocks") or ())
            fragments.append(FragmentMatch(
                fragment_index=int(raw.get("fragment_index", position)),
                island_index=int(raw.get("island_idx", position)),
                r_atoms=tuple(sorted(map(int, raw.get("fragment") or ()))),
                deferred_edges=tuple(
                    tuple(sorted(map(int, edge)))
                    for edge in raw.get("deferred_edges") or ()),
                symmetry_domains=tuple(domains),
                representative_assignments=tuple(
                    (int(source), int(target))
                    for source, target in dict(
                        symmetry.get("witness") or {}).items()),
                exact_fixed=tuple(map(
                    int, symmetry.get("exact_fixed") or ())),
                multiplicity=int(symmetry.get("multiplicity", 1)),
                automorph_domains=automorph_domains,
                target_generators=(
                    None if raw_generators is None else tuple(
                        AtomPermutation(tuple(map(int, generator)))
                        for generator in raw_generators))))
        return cls(tuple(fragments))

    @property
    def has_complete_exact_target_groups(self):
        return bool(self.fragments) and all(
            fragment.has_exact_target_group for fragment in self.fragments)

    def to_record(self):
        def domain_record(domain):
            return {
                "r_atoms": list(domain.r_atoms),
                "p_atoms": list(domain.p_atoms),
                "source": domain.source,
                "extendable": bool(domain.extendable),
            }

        fragments = []
        for fragment in self.fragments:
            symmetry = {
                "witness": dict(fragment.representative_assignments),
                "blocks": [
                    domain_record(domain)
                    for domain in fragment.symmetry_domains
                ],
                "exact_fixed": list(fragment.exact_fixed),
                "multiplicity": int(fragment.multiplicity),
                "automorph_blocks": [
                    domain_record(domain)
                    for domain in fragment.automorph_domains
                ],
            }
            if fragment.target_generators is not None:
                symmetry["automorph_generators"] = [
                    list(generator.images)
                    for generator in fragment.target_generators
                ]
            fragments.append({
                "fragment_index": int(fragment.fragment_index),
                "island_idx": int(fragment.island_index),
                "fragment": list(fragment.r_atoms),
                "deferred_edges": [
                    list(edge) for edge in fragment.deferred_edges],
                "symmetry": symmetry,
            })
        return {
            "rule": "typed_aam_hierarchy",
            "fragments": fragments,
            "blocks": [],
        }


@dataclass(frozen=True)
class AAMHierarchyView:
    """A base hierarchy plus one coordinate transform, without copied groups."""

    base: AAMHierarchy
    target_action: tuple[tuple[int, int], ...]

    @cached_property
    def materialized(self):
        return self.base._materialize_target(dict(self.target_action))

    @property
    def fragments(self):
        return self.materialized.fragments

    @property
    def has_complete_exact_target_groups(self):
        return self.base.has_complete_exact_target_groups

    def relabel_target(self, action):
        action, prior = dict(action), dict(self.target_action)
        if not action:
            return self
        combined = tuple(sorted((a, action.get(prior.get(a, a), prior.get(a, a)))
                                for a in set(prior) | set(action)))
        return AAMHierarchyView(self.base, combined)

    def to_record(self):
        """Explicitly requested, fully materialized legacy hierarchy record."""
        return self.materialized.to_record()


@dataclass(frozen=True)
class AAMBranch:
    """One unique completed hierarchical branch under a mechanism."""

    representative: AtomBijection
    hierarchy: AAMHierarchy
    encounter_count: int = 1
    cuts: tuple[tuple[int, int], ...] = ()
    covered_path_count: int = 1
    mapping_family: Mapping = field(default_factory=dict)
    path_provenance: tuple[Mapping, ...] = ()
    target_group: PermutationGroup | None = None

    def __post_init__(self):
        object.__setattr__(self, "encounter_count", int(self.encounter_count))
        object.__setattr__(self, "cuts", tuple(sorted(
            tuple(sorted(map(int, cut))) for cut in self.cuts)))
        object.__setattr__(
            self, "covered_path_count", int(self.covered_path_count))
        object.__setattr__(
            self, "mapping_family", dict(self.mapping_family or {}))
        object.__setattr__(self, "path_provenance", tuple(
            dict(item) for item in self.path_provenance))

    @property
    def has_exact_mapping_family(self):
        return self.target_group is not None


@dataclass(frozen=True)
class PostAAMMechanism:
    """Deduplicated AAM mechanism before chirality/RMSD post-processing."""

    mechanism_key: tuple
    representative: AtomBijection
    hierarchy: AAMHierarchy
    endpoint_source_symmetry: PermutationGroup
    endpoint_target_symmetry: PermutationGroup
    branches: tuple[AAMBranch, ...] = ()
    raw_branch_count: int = 1
    metadata: Mapping = field(default_factory=dict)

    def __post_init__(self):
        degree = self.representative.degree
        if (self.endpoint_source_symmetry.degree != degree
                or self.endpoint_target_symmetry.degree != degree):
            raise ValueError("mechanism groups and representative disagree")
        if any(branch.representative.degree != degree for branch in self.branches):
            raise ValueError("AAM branch mapping degree differs from mechanism")

    @classmethod
    def from_pool_entry(cls, mechanism_key, entry, *,
                        endpoint_source_symmetry=None,
                        endpoint_target_symmetry=None):
        """Create the decision object without reading concrete witnesses."""
        mapping = AtomBijection.from_mapping(entry["mapping"])
        symmetry = entry.get("branch_symmetry") or {}
        raw_branches = list(entry.get("branches") or ())
        branches = []
        for raw in raw_branches:
            mapping_family = dict(raw.get("mapping_family") or {})
            raw_group = raw.get("target_group_generators")
            if raw_group is None:
                raw_group = mapping_family.get("target_generators")
            family_representative = (
                mapping_family.get("representative_mapping")
                or raw.get("mapping") or entry["mapping"])
            branches.append(AAMBranch(
                representative=AtomBijection.from_mapping(
                    family_representative, mapping.degree),
                hierarchy=AAMHierarchy.from_record(
                    raw.get("hierarchy") or symmetry),
                encounter_count=int(raw.get("encounter_count", 1)),
                cuts=tuple(tuple(map(int, cut))
                           for cut in raw.get("cuts") or ()),
                covered_path_count=int(raw.get("covered_path_count", 1)),
                mapping_family=mapping_family,
                path_provenance=tuple(
                    dict(record)
                    for record in raw.get("path_provenance") or ()),
                target_group=(
                    None if raw_group is None else
                    PermutationGroup.from_generator_mappings(
                        mapping.degree, raw_group))))
        branches = tuple(branches)
        if not branches:
            branches = (AAMBranch(
                representative=mapping,
                hierarchy=AAMHierarchy.from_record(symmetry),
                encounter_count=int(entry.get("dedup_count", 1))),)
        return cls(
            mechanism_key=tuple(mechanism_key),
            representative=mapping,
            hierarchy=AAMHierarchy.from_record(symmetry),
            endpoint_source_symmetry=(
                endpoint_source_symmetry
                or PermutationGroup.trivial(mapping.degree)),
            endpoint_target_symmetry=(
                endpoint_target_symmetry
                or PermutationGroup.trivial(mapping.degree)),
            branches=branches,
            raw_branch_count=int(entry.get("dedup_count", 1)),
            metadata={"has_no_cut": bool(entry.get("has_no_cut", False))},
        )

    @classmethod
    def from_aam_graphs(cls, mechanism_key, entry, graph_R, graph_P, *,
                        symmetry_wbo_tolerance):
        """Construct exact endpoint groups directly from graph automorphisms."""
        from ..matcher.orbits import _nauty_atom_generators

        degree = len(entry["mapping"])
        source = PermutationGroup.from_generator_mappings(
            degree, _nauty_atom_generators(
                graph_R, wbo_tol=float(symmetry_wbo_tolerance)))
        target = PermutationGroup.from_generator_mappings(
            degree, _nauty_atom_generators(
                graph_P, wbo_tol=float(symmetry_wbo_tolerance)))
        return cls.from_pool_entry(
            mechanism_key, entry,
            endpoint_source_symmetry=source,
            endpoint_target_symmetry=target)

    def symmetry_record(self):
        """Serialize the typed analytical hierarchy and group availability."""
        return {
            "mechanism_key": self.mechanism_key,
            "representative_mapping": self.representative.as_dict(),
            "endpoint_source_generators": [
                [list(cycle) for cycle in generator.cycles()]
                for generator in self.endpoint_source_symmetry.generators
            ],
            "endpoint_target_generators": [
                [list(cycle) for cycle in generator.cycles()]
                for generator in self.endpoint_target_symmetry.generators
            ],
            "endpoint_source_orbits": [
                list(orbit)
                for orbit in self.endpoint_source_symmetry.orbits()],
            "endpoint_target_orbits": [
                list(orbit)
                for orbit in self.endpoint_target_symmetry.orbits()],
            "endpoint_symmetry_role": (
                "auxiliary_graph_symmetry_not_free_mapping_candidates"),
            "analytical_branches": [{
                "representative_mapping": branch.representative.as_dict(),
                "encounter_count": branch.encounter_count,
                "covered_path_count": branch.covered_path_count,
                "cuts": [list(cut) for cut in branch.cuts],
                "mapping_family": dict(branch.mapping_family),
                "path_provenance": [dict(record)
                                    for record in branch.path_provenance],
                "exact_mapping_family_available": (
                    branch.has_exact_mapping_family),
                "target_group_generators": (
                    None if branch.target_group is None else [
                        list(generator.images)
                        for generator in branch.target_group.generators
                    ]),
                "fragments": [{
                    "fragment_index": fragment.fragment_index,
                    "island_index": fragment.island_index,
                    "r_atoms": list(fragment.r_atoms),
                    "deferred_edges": [list(edge)
                                       for edge in fragment.deferred_edges],
                    "symmetry_domains": [{
                        "r_atoms": list(domain.r_atoms),
                        "p_atoms": list(domain.p_atoms),
                        "source": domain.source,
                    } for domain in fragment.symmetry_domains],
                    "exact_target_group_available": (
                        fragment.has_exact_target_group),
                    "target_generators": (
                        None if fragment.target_generators is None else [
                            list(generator.images)
                            for generator in fragment.target_generators
                        ]),
                } for fragment in branch.hierarchy.fragments],
            } for branch in self.branches],
            "fragments": [{
                "fragment_index": fragment.fragment_index,
                "island_index": fragment.island_index,
                "r_atoms": list(fragment.r_atoms),
                "deferred_edges": [list(edge) for edge in fragment.deferred_edges],
                "symmetry_domains": [{
                    "r_atoms": list(domain.r_atoms),
                    "p_atoms": list(domain.p_atoms),
                    "source": domain.source,
                } for domain in fragment.symmetry_domains],
                "exact_target_group_available": (
                    fragment.has_exact_target_group),
                "target_generators": (
                    None if fragment.target_generators is None else [
                        list(generator.images)
                        for generator in fragment.target_generators
                    ]),
            } for fragment in self.hierarchy.fragments],
            "raw_branch_count": self.raw_branch_count,
        }


@dataclass(frozen=True)
class ConstraintEvaluation:
    name: str
    valid: bool
    violations: int = 0
    metadata: Mapping = field(default_factory=dict)


class MappingConstraint(Protocol):
    name: str

    def evaluate(self, mechanism: PostAAMMechanism,
                 mapping: AtomBijection) -> ConstraintEvaluation: ...


class MappingObjective(Protocol):
    name: str

    def score(self, mechanism: PostAAMMechanism,
              mapping: AtomBijection) -> float: ...


@dataclass(frozen=True)
class OrientedSimplex:
    """One signed affine orientation constraint in R index order."""

    center_r: int
    neighbors_r: tuple[int, int, int]
    expected_sign: int
    degeneracy_tolerance: float = 0.1

    def __post_init__(self):
        object.__setattr__(self, "center_r", int(self.center_r))
        object.__setattr__(
            self, "neighbors_r", tuple(map(int, self.neighbors_r)))
        if len(self.neighbors_r) != 3:
            raise ValueError("oriented simplex requires exactly three vectors")
        if int(self.expected_sign) not in (-1, 1):
            raise ValueError("oriented simplex sign must be +1 or -1")
        object.__setattr__(self, "expected_sign", int(self.expected_sign))


def _simplex_sign(coords, origin, points, tolerance):
    xyz = np.asarray(coords, dtype=np.longdouble)
    vectors = np.stack([
        xyz[int(point)] - xyz[int(origin)] for point in points
    ])
    lengths = np.sqrt(np.sum(vectors * vectors, axis=1))
    denominator = np.prod(lengths)
    if denominator == 0:
        return 0
    normalized = np.linalg.det(np.asarray(vectors, dtype=float)) / float(
        denominator)
    if abs(normalized) <= float(tolerance):
        return 0
    return 1 if normalized > 0 else -1


@dataclass(frozen=True)
class AffineChiralityConstraint:
    """Mapping constraint independent of search branches and RMSD."""

    product_coords: np.ndarray
    simplices: tuple[OrientedSimplex, ...]
    name: str = "affine_index_chirality"

    def evaluate(self, mechanism: PostAAMMechanism,
                 mapping: AtomBijection):
        product = np.asarray(self.product_coords, dtype=float)
        violations = []
        for simplex in self.simplices:
            actual = _simplex_sign(
                product,
                mapping.images[simplex.center_r],
                tuple(mapping.images[r] for r in simplex.neighbors_r),
                simplex.degeneracy_tolerance)
            if actual != simplex.expected_sign:
                violations.append({
                    "center_r": simplex.center_r,
                    "neighbors_r": list(simplex.neighbors_r),
                    "expected_sign": simplex.expected_sign,
                    "actual_sign": actual,
                })
        return ConstraintEvaluation(
            name=self.name,
            valid=not violations,
            violations=len(violations),
            metadata={"violations": violations})


@dataclass(frozen=True)
class FixedMappingRMSD:
    """Proper-fit RMSD objective with immutable atom correspondence."""

    reactant_coords: np.ndarray
    product_coords: np.ndarray
    name: str = "fixed_mapping_proper_fit_rmsd"

    def score(self, mechanism: PostAAMMechanism,
              mapping: AtomBijection) -> float:
        reactant = np.asarray(self.reactant_coords, dtype=float)
        product = mapping.product_in_reactant_order(self.product_coords)
        if reactant.shape != product.shape:
            raise ValueError("RMSD endpoint coordinate shapes differ")
        centered_R = reactant - reactant.mean(axis=0)
        centered_P = product - product.mean(axis=0)
        u, _singular, vt = np.linalg.svd(centered_P.T @ centered_R)
        rotation = u @ vt
        if np.linalg.det(rotation) < 0.0:
            u[:, -1] *= -1.0
            rotation = u @ vt
        aligned = centered_P @ rotation
        return float(np.sqrt(np.mean(np.sum(
            (aligned - centered_R) ** 2, axis=1))))


@dataclass(frozen=True)
class MappingEvaluation:
    mapping: AtomBijection
    constraints: tuple[ConstraintEvaluation, ...]
    objective_name: str
    objective_value: float

    @property
    def valid(self):
        return all(result.valid for result in self.constraints)


@dataclass(frozen=True)
class PostAAMSelectionProblem:
    mechanism: PostAAMMechanism
    constraints: tuple[MappingConstraint, ...]
    objective: MappingObjective

    def evaluate(self, mapping: AtomBijection):
        if mapping.degree != self.mechanism.representative.degree:
            raise ValueError("candidate mapping has the wrong degree")
        results = tuple(
            constraint.evaluate(self.mechanism, mapping)
            for constraint in self.constraints)
        return MappingEvaluation(
            mapping=mapping,
            constraints=results,
            objective_name=self.objective.name,
            objective_value=self.objective.score(self.mechanism, mapping),
        )
