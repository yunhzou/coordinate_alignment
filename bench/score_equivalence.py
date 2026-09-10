"""Endpoint automorphisms that preserve event scoring for admissible atom maps."""
from collections import defaultdict
import numpy as np


def element_pair_features(raw, *, bond_floor=.2, event_tolerance=.5):
    """Compare a bond only with opposite-endpoint bonds of its element pair.

    Element-preserving AAM cannot compare C-C against C-H. Graph adjacency and
    every possible event response within each admissible pair are retained.
    These are score-equivalence colors, not exact WBO or geometric identity.
    """
    ends=[raw['reactant'],raw['product']]
    weights=[]
    for endpoint in ends:
        w=np.asarray(endpoint['wbo']);elements=endpoint['elements'];groups=defaultdict(set)
        for a,b in zip(*np.where(np.triu(w,1)>bond_floor)):
            groups[tuple(sorted((elements[a],elements[b])))].add(float(w[a,b]))
        weights.append(groups)
    result=[]
    for side,endpoint in enumerate(ends):
        w=np.asarray(endpoint['wbo']);elements=endpoint['elements'];bonds=[]
        for a,b in zip(*np.where(np.triu(w,1)>bond_floor)):
            values=sorted(weights[1-side][tuple(sorted((elements[a],elements[b])))])
            bonds.append((int(a),int(b),tuple(int(abs(w[a,b]-value)>event_tolerance) for value in values)))
        result.append(dict(colors=[(e,) for e in elements],bonds=bonds))
    return result
