"""Deterministic orthographic camera selection; coordinates are never altered.

Choose a view separating projected atom centers. This is only a figure-layout
operation and does not use mapping/reference labels or affect search.
"""
from functools import lru_cache
import numpy as np

@lru_cache(maxsize=8)
def _basis(coords):
    xyz=np.asarray(coords,dtype=float)
    xyz-=xyz.mean(axis=0)
    rng=np.random.default_rng(20260910)
    best=-np.inf;chosen=np.eye(3)
    ii,jj=np.triu_indices(len(xyz),1)
    for _ in range(1600):
        q,_=np.linalg.qr(rng.normal(size=(3,3)))
        q[:,2]=np.cross(q[:,0],q[:,1])
        xy=xyz@q[:,:2]
        scale=np.max(np.linalg.norm(xy,axis=1))
        d=np.linalg.norm(xy[ii]-xy[jj],axis=1)/scale
        score=np.min(d)+.08*np.mean(np.sort(d)[:len(xyz)])
        if score>best:best=score;chosen=q
    return chosen

def view_basis(ep):
    return _basis(tuple(tuple(x) for x in ep['coordinates']))

def projection(ep):
    xyz=np.asarray(ep['coordinates'],dtype=float)
    return (xyz-xyz.mean(axis=0))@view_basis(ep)[:,:2]
