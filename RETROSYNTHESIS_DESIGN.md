# Single-Step Retrosynthesis: Reference Design

This document is the normative algorithm and abstraction guide.  It separates
structural search from catalog I/O, recommendation policy, chemistry scoring,
and visualization.

## 1. Contract

Input:

- one target product graph `P`, represented with explicit hydrogens;
- a catalog of precursor graphs `R_i`; and
- bounded search policies.

Output:

- ranked precursor-set recommendations;
- atom ownership from every selected precursor into `P`;
- proposed broken and formed bonds; and
- explicit completeness/cap diagnostics at every bounded stage.

The output is a recommendation, not proof that a reaction is feasible.

## 2. Clean end-to-end pseudocode

```text
RETROSYNTHESIS(P, catalog, policies):
    target = PREPARE_TARGET(P, explicit_hydrogens = true)

    candidate_store = empty persistent shard set

    parallel for catalog record in catalog:
        precursor = PREPARE_PRECURSOR(record, explicit_hydrogens = true)

        if NECESSARY_COMPOSITION_FILTER_REJECTS(precursor, target):
            persist rejection reason
            continue

        result = DETECT_FRAGMENTS(
            source = precursor,
            target,
            policies.detection,
        )

        persist result, including cap diagnostics, immediately

    index = BUILD_CANDIDATE_INDEX(candidate_store, target)

    patterns = ENUMERATE_COVERAGE_PATTERNS(
        target,
        index.coverage_groups,
        policies.coverage_enumeration,
    )

    recommendations = empty bounded ranking

    for pattern in patterns:
        for precursor_set in ENUMERATE_PATTERN_ASSIGNMENTS(pattern, index):
            assembly = VALIDATE_ASSEMBLY(target, precursor_set)
            if assembly is invalid:
                continue

            score = STRUCTURAL_SCORE(assembly)
            recommendations.retain(assembly, score)

    diverse = SELECT_PATTERN_DIVERSE_RECOMMENDATIONS(recommendations)

    for finalist in diverse:
        refined = FINALIST_REFINEMENT(
            finalist,
            chirality = true,
            cut_sweep = true,
            chemical_checks = enabled,
        )
        persist refined result

    return RANK_REFINED_RECOMMENDATIONS(diverse)
```

The initial catalog pass should normally use chirality-insensitive, capped AAM
for recall.  Chirality-sensitive cut sweep belongs in finalist refinement.

## Part I: Detection

Detection is an independent graph capability.  Given any source graph and
target graph, it generates coherent source-fragment placements in the target.
It does not know about catalogs, precursor combinations, construction
patterns, recommendation ranking, or retrosynthesis.

Public package:

```text
rxn_core.fragment_matching
```

Public operation:

```text
detect_fragments(source, target, source_id, config) -> FragmentDetectionResult
```

### 3. Fragment candidate generation

```text
DETECT_FRAGMENTS(source, target, policy):
    diagnostics = new SearchDiagnostics(policy limits)
    initial_placements = empty deduplicated set

    for seed in ORDER_SEEDS(source):
        placements = GROW_CONNECTED_ISLAND(
            query = source,
            target = target,
            seed = seed,
            tolerance = policy.isomorphism_tolerance,
            branch_cap = policy.branch_cap,
        )

        diagnostics.observe(placements)
        if this seed hit its cap:
            diagnostics.mark_incomplete(stage = initial_growth)
            continue

        retain symmetry-distinct placements in initial_placements
        if candidate cap is reached:
            diagnostics.mark_incomplete(stage = candidate_collection)
            break

    best_initial_size = maximum connected size in initial_placements

    for initial in placements of best_initial_size:
        partition = CUT_INITIAL_FRAGMENT(source, initial.source_atoms)

        if partition violates configured fragmentation bounds:
            continue

        validations = COMPETITIVE_AUGMENTED_MATCH(
            source,
            target,
            initial,
            partition,
            policy,
        )

        diagnostics.observe(validations)

        for validation in validations:
            candidate = PROJECT_CANDIDATE(
                source, target, partition, validation)
            if candidate violates configured leftover bounds:
                continue
            emit candidate unless structurally equivalent to one already emitted

    return SearchResult(candidates, diagnostics)
```

The maximum applies only to the initial connected island.  Competitive
augmentation may discover additional target-owned fragments from the same
precursor.

### 4. Competitive augmented matching

```text
COMPETITIVE_AUGMENTED_MATCH(R, P, initial, partition, policy):
    cut_query = copy R
    remove every initial-fragment boundary bond from cut_query

    residual_copies = disconnected copies of atoms outside the initial fragment
    augmented_target = disjoint_union(P, residual_copies)

    matches = COMPLETE_SUBGRAPH_MATCH(
        query = cut_query,
        target = augmented_target,
        fixed_mapping = initial R-to-P mapping,
        tolerance = policy.isomorphism_tolerance,
        branch_cap = policy.branch_cap,
    )

    valid = []
    for match in matches:
        if match changes the fixed initial mapping:
            continue
        if a cut boundary has both endpoints in P but they are not adjacent in P:
            continue
        valid.append(match)

    best_target_ownership = maximum number of R atoms mapped into original P
    return all valid matches having best_target_ownership
```

Residual copies are competitors, never anchors.  A residual component may map
into unused atoms of `P` when it fits; otherwise it remains on its appended
copy.  This is why carbon dioxide can contribute a carbonyl fragment plus
a second oxygen fragment to a carboxyl group at tolerance 0.5.

### 5. Candidate projection

```text
PROJECT_CANDIDATE(R, P, partition, augmented_match):
    target_mapping = mappings whose image is an original atom of P
    copied_residual_mapping = mappings whose image is an appended copy

    target_owned_atoms_R = domain(target_mapping)
    retained_fragments = connected components of target_owned_atoms_R
                         after the recorded boundary cuts
    leftover_fragments = connected components of all other R atoms
                         after the recorded boundary cuts

    return FragmentCandidate(
        source identity,
        target mapping,
        retained fragments,
        leftover fragments,
        boundary cuts,
        attachment atoms,
        copied residual mapping,
    )
```

A candidate is one source contribution and may contain multiple retained
fragments.  It is not an assembly and it is not identified by coverage alone.

## Part II: Assembly

Assembly is the retrosynthesis-specific layer.  It consumes detected fragment
candidates, assigns candidates to precursor records, constructs complete
non-overlapping target covers, and ranks precursor-set recommendations.

Public package:

```text
rxn_core.retrosynthesis
```

Detection never imports or calls assembly.  Assembly depends on detection's
`FragmentCandidate` record.

Neither API is re-exported from the root `rxn_core` namespace.  Callers import
the component they use explicitly:

```text
from rxn_core.fragment_matching import detect_fragments
from rxn_core.retrosynthesis import assemble_fragment_cover
```

### 6. Candidate indexing

```text
BUILD_CANDIDATE_INDEX(shards, P):
    validate schema and target identity for every shard
    canonicalize precursor structure identity
    optionally create controlled ownership variants at attachment atoms
    calculate explicit-H and heavy-atom retention statistics
    calculate chirality diagnostics only when requested

    group candidates by exact target coverage mask
    retain a bounded, ranked, structurally diverse pool per mask

    return CandidateIndex(groups, scan diagnostics, cap diagnostics)
```

Catalog parsing, gzip/JSON serialization, and multiprocessing belong outside
the chemistry/search core.

### 7. Coverage-pattern enumeration

```text
ENUMERATE_COVERAGE_PATTERNS(P, groups, policy):
    full = mask containing every explicit atom of P
    states = { empty coverage: empty pattern }

    repeat up to policy.maximum_precursor_copies times:
        next_states = empty map from coverage to diverse partial patterns

        for state in states:
            pivot = choose one uncovered target atom
            for disjoint candidate mask containing pivot:
                expanded = state plus candidate mask

                if expanded covers full:
                    emit its construction-pattern signature
                else if remaining masks can still complete expanded:
                    retain bounded diverse paths for this coverage state

        if global state cap is exceeded:
            prune by optimistic structural rank
            mark enumeration incomplete

        states = next_states
```

Repeated copies of the same precursor are allowed by policy.  The objective is
not to minimize copy count independently; it first favors fewer unique
precursor structures and then high set retention.

### 8. Assembly validation

```text
VALIDATE_ASSEMBLY(P, selected_candidates):
    owner = empty map from target atom to selected candidate

    for candidate in selected_candidates:
        reject if it contributes no target atom
        reject if any target atom already has an owner
        assign ownership for its entire candidate coverage

    reject unless every explicit atom of P has exactly one owner

    formed_bonds = every P bond whose endpoints have different owners
    if strict attachment validation is enabled:
        reject a formed bond unsupported by both candidates' attachment atoms

    broken_bonds = union of selected candidates' recorded boundary cuts
    return Assembly(selected candidates, formed_bonds, broken_bonds)
```

Overlap is handled by generating explicit ownership variants before assembly,
not by allowing two selected candidates to own one product atom.

### 9. Ranking and diversity

Structural ranking is lexicographic and transparent:

```text
STRUCTURAL_SCORE(assembly):
    return (
        chirality violations,              # finalist stage when enabled
        number of unique precursor structures,
        negative all-atom set retention,   # explicit H included
        number of capped precursor results,
        number of broken precursor bonds,
        number of leftover precursor atoms,
        number of formed product bonds,
        stable identity tie-breaker,
    )
```

Recommendations are diversified by construction pattern before taking the
top variants.  A construction pattern is the isomorphism-aware partition of
target atoms into precursor-owned modules, not merely a list of precursor IDs.

### 10. Finalist refinement

```text
FINALIST_REFINEMENT(assembly):
    rerun careful AAM with chirality enabled
    run cut sweep with separately reported caps
    validate valence, charge, and stereochemical conservation
    estimate bond-edit count and reaction-center coherence
    optionally apply reaction-family or learned feasibility models
    keep original coarse result and refined result linked by stable IDs
```

Refinement must never overwrite the coarse search artifact.  Every expensive
stage consumes saved intermediates and writes a new versioned artifact.

## 11. Module boundaries

```text
rxn_core.fragment_matching.models
    detection config, fragment candidates, results, cap diagnostics

rxn_core.fragment_matching.detection
    connected-island growth and fragment candidate generation

rxn_core.fragment_matching.augmentation
    boundary cutting, augmented target construction, candidate projection

rxn_core.fragment_matching.serialization
    strict fragment-detection artifact conversion

rxn_core.fragment_matching.rdkit_adapter
    optional RDKit-to-weighted-graph conversion

rxn_core.retrosynthesis.models
    assembly result records

rxn_core.retrosynthesis.coverage
    in-memory assembly validation for domain candidates

rxn_core.retrosynthesis.catalog_index
    persisted-record normalization, controlled ownership variants, mask index

rxn_core.retrosynthesis.enumeration
    bounded exact, modular, and recommendation coverage-pattern search

rxn_core.retrosynthesis.ranking
    candidate and assembly score construction; no I/O

tools/search_mcule_retro.py
    catalog reader and parallel adapter around fragment detection

tools/merge_retro_catalog.py
    CLI parsing and composition of coverage/ranking services

viewer tools
    read saved artifacts only; never infer or manufacture mappings
```

Dependencies point inward: CLIs and viewers depend on the domain package;
domain algorithms never depend on a CLI, RDKit file format, JSON, or viewer.

## 12. Non-negotiable invariants

- Hydrogens are explicit in search, coverage, retention, and visualization.
- A target atom has at most one owner in a completed assembly.
- Appended residual copies are competitors, not hard anchors.
- Target-mapped residual fragments respect the original cut attachment.
- Branch-cap hits always propagate as `complete = false` with stage details.
- Repeated precursor copies are representable without duplicating structure
  identity or colors.
- Chirality-insensitive screening cannot reject a candidate solely because of
  stereochemistry; chirality-sensitive refinement may rerank or reject it.
- Saved mappings are the sole source of viewer colors and atom provenance.
- No expensive result is discarded before its versioned intermediate is
  written.
