"""Experimental symmetry closure; can add admissible mappings. Not pure deduplication."""
from collections import defaultdict
from itertools import combinations
from rxn_core.final_branches import FinalBranchCatalogue, FinalFamily, canonical_actions

class _TargetPartitionCanonicalizer:
    """Align equivalent final product partitions using exact target symmetry."""
    def __init__(self, saved, policy):
        import numpy as np
        self.saved, self.policy = saved, policy
        self.n = saved.problem.target_atom_count
        whole = (tuple(range(saved.problem.source_atom_count)), tuple(range(self.n)))
        self.symmetry = saved.rebuild_symmetry(policy, whole)
        self.labels = saved._labels[policy, whole]
        self.generators = self.symmetry['generators']
        self.elements = saved.problem.product.elements
        self.adjacency = {i: set() for i in range(self.n)}
        self.bond_colors = defaultdict(set)
        next_vertex = self.n
        for a, b in combinations(range(self.n), 2):
            label = int(self.labels[a, b])
            if not label:
                continue
            self.bond_colors[('bond', label)].add(next_vertex)
            self.adjacency[next_vertex] = {a, b}
            self.adjacency[a].add(next_vertex)
            self.adjacency[b].add(next_vertex)
            next_vertex += 1
        self.adjacency = {a: sorted(v) for a, v in self.adjacency.items()}
        self.vertices = next_vertex
        self.references = {}
        self.cache = {}
        self.membership = {}

    def contains(self, factor):
        import numpy as np
        if factor in self.membership:
            return self.membership[factor]
        kind, data = factor
        if kind == 'group':
            generators = data
        else:
            generators = []
            for a, b in zip(data, data[1:]):
                images = list(range(self.n));images[a], images[b] = b, a
                generators.append(images)
        good = all(all(self.elements[a] == self.elements[b] for a, b in enumerate(g))
                   and np.array_equal(self.labels, self.labels[np.ix_(g, g)]) for g in generators)
        self.membership[factor] = good
        return good

    def align(self, pairs):
        import numpy as np
        import pynauty
        if pairs in self.cache:
            return self.cache[pairs]
        colors = defaultdict(set)
        for i, (rs, ps) in enumerate(pairs):
            for p in ps:
                colors[('atom', self.elements[p], i)].add(p)
        colors.update(self.bond_colors)
        color_keys = tuple(sorted(colors, key=repr))
        graph = pynauty.Graph(self.vertices, adjacency_dict=self.adjacency,
                             vertex_coloring=[colors[k] for k in color_keys])
        key = tuple(rs for rs, ps in pairs), color_keys, pynauty.certificate(graph)
        order = pynauty.canon_label(graph)
        if key not in self.references:
            self.references[key] = pairs, order
        reference_pairs, reference_order = self.references[key]
        full_map = dict(zip(order, reference_order))
        images = tuple(full_map[p] for p in range(self.n))
        assert all(p < self.n for p in images)
        assert np.array_equal(self.labels, self.labels[np.ix_(images, images)])
        assert tuple((rs, tuple(sorted(images[p] for p in ps))) for rs, ps in pairs) == reference_pairs
        result = reference_pairs, images
        self.cache[pairs] = result
        return result


def rebuild_intrinsic_catalogue(saved):
    """Canonical final pairs, intrinsic shuffles, and common target symmetry.

    First align equivalent product placements under whole-target matching
    symmetry. Then use each final pair's intrinsic symmetry, followed by the
    common whole-target action. Required bonds remain enforced. Any old family
    not certified covered is retained explicitly, closed under the same safe
    whole-target symmetry. The reconstructed union contains every saved family
    and can add further matching-admissible shuffles. No permutations or group
    elements are enumerated.
    """
    import numpy as np
    result = FinalBranchCatalogue(saved.problem)
    result.mode = 'intrinsic_final_pair_reconstruction'
    result.path_count, result.incomplete_paths = saved.path_count, saved.incomplete_paths
    canonicalizers, local_cache, orbit_ids, fallback_ids = {}, {}, {}, {}
    fallback_inputs, edge_fallback = set(), set()
    n = saved.problem.target_atom_count
    for (policy, pairs), family_ids in saved.branches.items():
        if policy not in canonicalizers:
            canonicalizers[policy] = _TargetPartitionCanonicalizer(saved, policy)
        target = canonicalizers[policy]
        reference_pairs, alignment = target.align(pairs)
        branch_key = policy, reference_pairs
        branch = result.branches.setdefault(branch_key, set())
        inverse = tuple(alignment.index(p) for p in range(n))
        owner = {r: i for i, (rs, ps) in enumerate(pairs) for r in rs}
        factors = []
        indices = []
        for pair in reference_pairs:
            symmetry = saved.rebuild_symmetry(policy, pair)
            indices.append({p: i for i, p in enumerate(pair[1])})
            generators = []
            for local in symmetry['generators']:
                images = list(range(n))
                for i, p in enumerate(pair[1]):
                    images[p] = pair[1][local[i]]
                generators.append(tuple(images))
            if generators:
                factors.append(result._intern(result._factors, ('group', tuple(sorted(generators)))))
        local_program = canonical_actions(factors)
        global_factor = ('group', target.generators) if target.generators else None

        def with_global(program):
            program = list(program)
            # Adjacent subgroups of the final global action are redundant.
            if global_factor:
                while program and target.contains(program[-1]):
                    program.pop()
                program.append(global_factor)
            return result._intern(result._programs, canonical_actions(program))

        rebuilt_program = with_global(local_program)

        def conjugate(program):
            factors = []
            for kind, data in program:
                if kind == 'pool':
                    factor = kind, tuple(sorted(alignment[p] for p in data))
                else:
                    factor = kind, tuple(sorted(tuple(alignment[g[inverse[p]]] for p in range(n))
                                                for g in data))
                factors.append(result._intern(result._factors, factor))
            return canonical_actions(factors)

        def locally_covered(program):
            cache_key = branch_key, program
            if cache_key in local_cache:
                return local_cache[cache_key]
            def preserves(images):
                for pair, local in zip(reference_pairs, indices):
                    if any(images[p] not in local for p in pair[1]):
                        return False
                    permutation = [local[images[p]] for p in pair[1]]
                    labels = saved._labels[policy, pair]
                    if not np.array_equal(labels, labels[np.ix_(permutation, permutation)]):
                        return False
                return True
            good = True
            for kind, data in program:
                if kind == 'group':
                    generators = data
                else:
                    generators = []
                    for a, b in zip(data, data[1:]):
                        images = list(range(n));images[a], images[b] = b, a
                        generators.append(images)
                if not all(preserves(g) for g in generators):
                    good = False
                    break
            local_cache[cache_key] = good
            return good

        for family_id in sorted(family_ids):
            family = saved.families[family_id]
            mapping = tuple((r, alignment[p]) for r, p in family.mapping)
            by_r = dict(mapping)
            cross_edges = any(owner[a] != owner[b] for a, b in family.required_edges)
            if not cross_edges:
                pullbacks = []
                for pair, local in zip(reference_pairs, indices):
                    permutation = [local[by_r[r]] for r in pair[0]]
                    labels = saved._labels[policy, pair]
                    pullbacks.append(labels[np.ix_(permutation, permutation)].tobytes())
                key = branch_key, family.required_edges, tuple(pullbacks)
                if key not in orbit_ids:
                    orbit_ids[key] = len(result.families)
                    result.families.append(FinalFamily(mapping, family.required_edges,
                                                       rebuilt_program, policy))
                out_id = orbit_ids[key]
                branch.add(out_id)
                result.families[out_id].provenance.extend(family.provenance)
                if all(target.contains(factor) for factor in family.actions):
                    continue
                aligned_program = conjugate(family.actions)
                if locally_covered(aligned_program):
                    continue
            else:
                edge_fallback.add(family_id)
                aligned_program = conjugate(family.actions)
            fallback_inputs.add(family_id)
            program = with_global(aligned_program)
            key = policy, mapping, family.required_edges, program
            if key not in fallback_ids:
                fallback_ids[key] = len(result.families)
                result.families.append(FinalFamily(mapping, family.required_edges, program, policy))
            out_id = fallback_ids[key]
            result.families[out_id].provenance.extend(family.provenance)
            # An unproved cross-pair action is an explicit coupled branch, not
            # silently treated as an independent fixed fragment combination.
            result.coupled_branches.setdefault(out_id, set()).add(reference_pairs)
    result.branches = {k: v for k, v in result.branches.items() if v}
    result.symmetries = saved.symmetries
    result._labels = saved._labels
    result.reconstruction = dict(saved_families=len(saved.families),
        rebuilt_orbits=len(orbit_ids), retained_original_families=len(fallback_inputs),
        canonical_fallback_families=len(fallback_ids),
        cross_fragment_constraint_fallbacks=len(edge_fallback),
        coupled_branches=len(result.coupled_branches),
        whole_target_partition_classes=sum(len(c.references) for c in canonicalizers.values()),
        union_contains_every_saved_family=True,
        permits_additional_matching_admissible_shuffles=True)
    return result
