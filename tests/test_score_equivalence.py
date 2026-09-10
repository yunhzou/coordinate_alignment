import random
import numpy as np
import pynauty
from golden_evaluation import colored_graph
from compare_elementary_outputs import event_counts
from score_equivalence import element_pair_features


def test_pair_response_generators_preserve_all_event_counts_for_element_maps():
    elements=['C','C','C','H','H','H']
    r=np.zeros((6,6));p=np.zeros((6,6))
    for matrix,weights in ((r,[1.47,1.48,1.49,.97,.98,.99]),(p,[1.49,1.47,1.48,.99,.97,.98])):
        for (a,b),w in zip(((0,1),(1,2),(2,0),(0,3),(1,4),(2,5)),weights):matrix[a,b]=matrix[b,a]=w
    raw=dict(reactant=dict(elements=elements,wbo=r),product=dict(elements=elements,wbo=p))
    features=element_pair_features(raw)
    generators=[tuple(tuple(g[:6]) for g in pynauty.autgrp(colored_graph([feature]))[0]) for feature in features]
    assert any(any(g[i]!=i for i in range(3)) for g in generators[0])
    rng=random.Random(812)
    for _ in range(40):
        carbon=[0,1,2];hydrogen=[3,4,5];rng.shuffle(carbon);rng.shuffle(hydrogen)
        mapping=carbon+hydrogen
        score=event_counts(r,p,[mapping])[0]
        for g in generators[0]:np.testing.assert_array_equal(event_counts(r,p,[[mapping[g[i]] for i in range(6)]])[0],score)
        for g in generators[1]:np.testing.assert_array_equal(event_counts(r,p,[[g[mapping[i]] for i in range(6)]])[0],score)
