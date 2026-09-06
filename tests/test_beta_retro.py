"""The beta is an opt-in scheduling policy, not a change to AAM semantics."""
from dataclasses import replace

from rxn_core.fragment_matching import FragmentCandidate, FragmentDetectionConfig
from rxn_core.fragment_matching.connected import find_connected_fragments
from rxn_core.retrosynthesis.beta import (
    BetaPlacement, FragmentQueryBank, recommend_big_blocks,
)
from rxn_core.smiles import smiles_to_weighted_graph


def graph(smiles):
    return smiles_to_weighted_graph(smiles, expand_hydrogens=True)


def block(source, atoms, *, refined=False):
    mapping = tuple(enumerate(atoms))
    candidate = FragmentCandidate(source, mapping, tuple(range(len(atoms))),
        tuple(atoms), (), (), (), (), (), 4,
        retained_fragments=(tuple(range(len(atoms))),))
    return BetaPlacement(candidate, tuple(range(4)), refined)


class FakeBank:
    target = range(4)
    capped_searches = 0

    def __init__(self):
        self.events = []

    ordered_query = FragmentQueryBank.ordered_query


def test_selected_source_is_refined_before_gap_bank_search():
    class Bank(FakeBank):
        def query(self, region):
            self.events.append(('query', frozenset(region)))
            return (block('R', (0, 1, 2)),)
        def refine(self, selected):
            self.events.append(('refine', selected.candidate.source_id))
            return (block('R', (0, 1, 2, 3), refined=True),)
    bank = Bank()
    result = recommend_big_blocks(bank)
    assert bank.events == [('query', frozenset(range(4))), ('refine', 'R')]
    assert len(result.recommendations[0].placements) == 1
    assert not result.exhaustive


def test_gap_queries_and_repeated_reactant_copies_with_overlap():
    class Bank(FakeBank):
        def query(self, region):
            self.events.append(frozenset(region))
            return (block('R', (0, 1, 2)),) if len(region) == 4 else (block('R', (3,)),)
        def refine(self, selected):
            if 3 in selected.covered_atoms:
                return (block('R', (2, 3), refined=True),)
            return (replace(selected, refined=True),)
    bank = Bank()
    result = recommend_big_blocks(bank)
    assert bank.events == [frozenset(range(4)), frozenset({3})]
    chosen = result.recommendations[0].placements
    assert len(chosen) == 2
    assert {p.candidate.source_id for p in chosen} == {'R'}
    assert chosen[0].covered_atoms & chosen[1].covered_atoms == {2}


def test_stalled_largest_choice_does_not_delete_alternatives():
    class Bank(FakeBank):
        def query(self, region):
            return (block('large', (0, 1, 2)), block('small', (0, 1)))
        def refine(self, selected):
            self.events.append(selected.candidate.source_id)
            return () if selected.candidate.source_id == 'large' else (
                block('small', (0, 1, 2, 3), refined=True),)
    bank = Bank()
    result = recommend_big_blocks(bank)
    assert bank.events == ['large', 'small']
    assert result.recommendations[0].placements[0].candidate.source_id == 'small'


def test_no_cover_is_an_explicit_partial_result():
    class Bank(FakeBank):
        def query(self, region):
            return (block('R', (0, 1, 2)),) if len(region) == 4 else ()
        def refine(self, selected):
            return (replace(selected, refined=True),)
    result = recommend_big_blocks(Bank())
    assert not result.recommendations
    assert result.best_partial.uncovered_target_atoms == (3,)


def test_connected_stage_never_augments_and_preserves_search_evidence(monkeypatch):
    import rxn_core.fragment_matching.detection as detection
    def forbidden(*args, **kwargs):
        raise AssertionError('bank-wide connected scan must not augment')
    monkeypatch.setattr(detection, '_augment_initial_family', forbidden)
    result = find_connected_fragments(graph('CCBr'), graph('CCO'), source_id='ethyl',
        config=FragmentDetectionConfig(iso_tolerance=1.0))
    assert result.candidates
    assert result.search_graphs
    assert all(len(c.retained_fragments) == 1 and c.derivations for c in result.candidates)


def test_real_explicit_h_methanol_cover_and_checkpointing():
    saved = []
    bank = FragmentQueryBank([('methane', graph('C')), ('water', graph('O'))], graph('CO'),
        checkpoint=lambda event, result: saved.append((event, result)))
    result = recommend_big_blocks(bank)
    assert result.recommendations
    chosen = result.recommendations[0].placements
    assert set().union(*(p.covered_atoms for p in chosen)) == set(range(6))
    assert {p.candidate.source_id for p in chosen} == {'methane', 'water'}
    assert all(p.refined for p in chosen)
    assert saved and all(evidence.search_graphs for _, evidence in saved)
    for placed in chosen:
        mapping = dict(placed.mapping)
        assert len(set(mapping.values())) == len(mapping)
        source = bank.sources[placed.candidate.source_id]
        for a, p in mapping.items():
            assert source.nodes[a]['element'] == bank.target.nodes[p]['element']
        for a, b in placed.candidate.preserved_source_bonds:
            assert bank.target.has_edge(mapping[a], mapping[b])


def test_gap_local_indices_and_hydrogen_are_preserved():
    bank = FragmentQueryBank([('water', graph('O'))], graph('CO'))
    hydrogen = next(a for a, data in bank.target.nodes(data=True) if data['element'] == 'H')
    options = bank.query({hydrogen})
    assert options and all(p.covered_atoms == {hydrogen} for p in options)
    assert all(p.target_atoms == (hydrogen,) for p in options)
    event_count = len(bank.events)
    assert bank.query({hydrogen}) is options
    assert len(bank.events) == event_count


def test_caps_are_not_reported_as_complete(monkeypatch):
    import rxn_core.fragment_matching.connected as connected
    monkeypatch.setattr(connected, '_initial_fragment_placements',
        lambda *args, **kwargs: ((), 1, 101, False, False, 1, 0, False, ()))
    result = find_connected_fragments(graph('C'), graph('C'))
    assert result.capped_seed_count == 1
    assert not result.complete


def test_parallel_connected_scan_preserves_occupations():
    sources = [('methane', graph('C')), ('water', graph('O'))]
    serial = FragmentQueryBank(sources, graph('CO'))
    parallel = FragmentQueryBank(sources, graph('CO'), workers=2)
    assert {p.key for p in serial.query(range(6))} == {
        p.key for p in parallel.query(range(6))}


def test_sorted_disk_stream_is_consumed_lazily_without_pruning_alternatives():
    class Bank(FakeBank):
        def ordered_query(self, region, placements):
            for i in range(1000):
                self.events.append(i)
                yield block(f'R{i:04d}', (0, 1, 2))
        def refine(self, selected):
            return (block(selected.source_id, (0, 1, 2, 3), refined=True),)
    bank = Bank()
    result = recommend_big_blocks(bank)
    assert result.recommendations
    assert bank.events == [0, 1]
