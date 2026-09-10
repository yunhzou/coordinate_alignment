"""Verify an existing representative certificate without rescoring its whole DAG."""
import time

import pynauty

from golden_evaluation import colored_graph, project
from publication_analysis import certificate_id


def check_ranked_reference(classes, plan, reference, search_row):
    """First stage of the existing representative-then-symbolic verification.

None means a compressed-family query is still needed, not non-recovery.
Class keys contain the exact saved full-atom witness, not an approximate score.
"""
    start = time.perf_counter()
    features, wanted = reference['features'], reference['mapping']
    expected = certificate_id(features,wanted)
    hit = next((row for row in classes if row['id']==expected),None)
    if hit is None:
        return None
    mapping = dict(hit['key'][3])
    assert pynauty.certificate(colored_graph(features,project(mapping,features))) == \
           pynauty.certificate(colored_graph(features,project(wanted,features)))
    heavy_target = len(features[1]['heavy'])
    top = classes[0]
    return dict(reference_pairs=len(wanted),reference_annotation_complete=len(wanted)==heavy_target,
        top1_correct=top['id']==expected,representative_recovery=True,reference_recovery='recovered',
        top_terminal=top['terminal'],witness_terminal=hit['terminal'],
        candidate_terminals=search_row['terminals'],unique_representative_chemistries=len(classes),
        best_target_heavy_coverage=max(-row['key'][0] for row in classes)/max(1,heavy_target),
        best_target_all_atom_coverage=max(-row['key'][1] for row in classes)/plan.input_problem.target_atom_count,
        capped=search_row['capped'],symbolic_queries=0,unknown_queries=0,
        evaluation_seconds=time.perf_counter()-start,top_events=dict(
            broken=top['events']['broken'],formed=top['events']['formed'],
            bond_order_changed=top['events']['order_changed']),
        input_orientation_witness=sorted(mapping.items()),search_direction=plan.direction,
        mapping_space='Witness is input R/P; terminal ID refers to the saved search DAG',
        coverage_space='original input product',verifier='saved_ranked_full_atom_witness_v1',
        reference_semantics='exact heavy-atom relation including unmatched atoms, modulo endpoint chemical symmetry')
