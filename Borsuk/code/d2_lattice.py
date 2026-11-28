#!/usr/bin/env python3
"""D2: lattice subsets.  For a lattice L and a squared distance s, a subset S of L with diam(S)^2 = s may be
translated to contain 0; then S lies in L_s = {x in L : |x|^2 <= s} and its diameter graph is an induced subgraph
of H_s = (L_s, pairs at squared distance exactly s).  Step 1: 5-colour H_s (tabu; SAT).  Step 2 if needed: CEGAR
over subsets S of L_s containing 0 with all pairs at squared distance <= s (as in d1_cegar.py).  Exact integers.
Lattices: Z4; D4 (even coordinate sum); A4 (Z^5, coordinate sum 0); A4* (5x - (sum x)1, x in Z^5, scaled by 5).
Usage: python d2_lattice.py LATTICE SMAX
"""
import sys, json, time, itertools
import numpy as np
from pysat.solvers import Solver
sys.path.insert(0, 'Borsuk/code')
from tabu import colour_tabu
sys.path.insert(0, 'Borsuk/code'); import cegar as cg

def points(name, smax):
    if name in ('Z4', 'D4'):
        b = int(np.floor(np.sqrt(smax)))
        P = np.array(list(itertools.product(range(-b, b + 1), repeat=4)), np.int64)
        if name == 'D4': P = P[P.sum(1) % 2 == 0]
    elif name == 'A4':
        b = int(np.floor(np.sqrt(smax)))
        Q = np.array(list(itertools.product(range(-b, b + 1), repeat=4)), np.int64)
        P = np.hstack([Q, -Q.sum(1, keepdims=True)])
    elif name == 'A4*':     # 5x - (sum x) 1 for x in Z^5, x5 = 0 representatives cover Z^5 / Z(1,...,1)
        b = int(np.floor(np.sqrt(smax) / 4)) + 2
        Q = np.array(list(itertools.product(range(-b, b + 1), repeat=4)), np.int64)
        Q = np.hstack([Q, np.zeros((len(Q), 1), np.int64)])
        P = 5 * Q - Q.sum(1, keepdims=True)
        P = np.unique(P, axis=0)
    n2 = (P * P).sum(1); P = P[n2 <= smax]; n2 = n2[n2 <= smax]
    o = np.argsort(n2, kind='stable'); return P[o]           # P[0] = origin

def sqd(A, B): return (A * A).sum(1)[:, None] + (B * B).sum(1)[None] - 2 * sum(A[:, k, None] * B[None, :, k] for k in range(A.shape[1]))

def solve(P, s, tl=3600):
    n2 = (P * P).sum(1); L = P[n2 <= s]; m = len(L)
    D = np.empty((m, m), np.int16)
    for i0 in range(0, m, 256): D[i0:i0 + 256] = sqd(L[i0:i0 + 256], L)
    E = np.triu(D == s, 1); U, V = [a.astype(np.int32) for a in np.nonzero(E)]
    if len(U) == 0 or np.bincount(np.concatenate([U, V]), minlength=m).max() < 5: return {'s': int(s), 'L': m, 'res': 'deg<5'}
    for sd in range(1, 4):
        col, c = colour_tabu(m, U, V, 5, iters=400000 * sd, seed=sd)
        if c == 0:
            assert (col[U] != col[V]).all(); return {'s': int(s), 'L': m, 'edges': len(U), 'res': 'tabu'}
    BAD = D > s
    res = cg.cegar(m, U, V, BAD, tl=tl, root=0)
    out = {'s': int(s), 'L': m, 'edges': len(U), 'res': 'CEGAR ' + res[0]}
    if res[0] == 'CANDIDATE': out['S'] = L[res[1]].tolist()
    return out

if __name__ == '__main__':
    name, smax = sys.argv[1], int(sys.argv[2]); P = points(name, smax)
    vals = sorted(set(((P * P).sum(1)).tolist()) - {0})   # lattice: squared distances = norms
    out = open('Borsuk/.scratch/d2_%s.jsonl' % name.replace('*', 's'), 'w')
    for s in vals:
        r = {'lattice': name, **solve(P, s)}; print(json.dumps(r), flush=True); out.write(json.dumps(r) + '\n'); out.flush()
