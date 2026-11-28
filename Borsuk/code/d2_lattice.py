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
sys.path.insert(0, 'Borsuk/.scratch')
from tabu import colour_tabu

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

def solve(P, s, tl=1200):
    n2 = (P * P).sum(1); L = P[n2 <= s]; m = len(L)
    D = np.empty((m, m), np.int16)
    for i0 in range(0, m, 256): D[i0:i0 + 256] = sqd(L[i0:i0 + 256], L)
    E = np.triu(D == s, 1); U, V = [a.astype(np.int32) for a in np.nonzero(E)]
    if len(U) == 0 or np.bincount(np.concatenate([U, V]), minlength=m).max() < 5: return {'s': int(s), 'L': m, 'res': 'deg<5'}
    for sd in range(1, 4):
        col, c = colour_tabu(m, U, V, 5, iters=400000 * sd, seed=sd)
        if c == 0:
            assert (col[U] != col[V]).all(); return {'s': int(s), 'L': m, 'edges': len(U), 'res': 'tabu'}
    # CEGAR (selection of compatible S containing 0 = index 0)
    bad = np.triu(D > s, 1); A = Solver(name='g4'); A.add_clause([1])
    for a, b in zip(*np.nonzero(bad)): A.add_clause([-(int(a) + 1), -(int(b) + 1)])
    compat = ~(bad | bad.T); adj = [[] for _ in range(m)]
    for a, b in zip(U.tolist(), V.tolist()): adj[a].append(b); adj[b].append(a)
    nxt = m + 1; t0 = time.time(); it = 0
    while True:
        it += 1
        if time.time() - t0 > tl: return {'s': int(s), 'L': m, 'res': 'TIMEOUT', 'iters': it}
        if not A.solve(): return {'s': int(s), 'L': m, 'edges': len(U), 'res': 'CEGAR UNSAT (dead)', 'iters': it}
        mod = A.get_model(); inS = np.zeros(m, bool); inS[[v for v in range(m) if mod[v] > 0]] = True
        for v in np.random.default_rng(it).permutation(m):
            if not inS[v] and compat[v, inS].all(): inS[v] = True
        S = np.flatnonzero(inS); idx = -np.ones(m, int); idx[S] = np.arange(len(S))
        sel = inS[U] & inS[V]; Us, Vs = idx[U[sel]].astype(np.int32), idx[V[sel]].astype(np.int32)
        col = np.zeros(len(S), int) if len(Us) == 0 else None
        for sd in range(1, 4):
            if col is not None: break
            c_, conf = colour_tabu(len(S), Us, Vs, 5, iters=300000 * sd, seed=sd)
            if conf == 0: col = c_
        if col is None:
            return {'s': int(s), 'L': m, 'res': 'CANDIDATE (tabu failed on compatible S)', 'S': L[S].tolist()}
        full = -np.ones(m, int); full[S] = col
        for v in range(m):
            if full[v] < 0:
                cnt = np.zeros(5, int)
                for w in adj[v]:
                    if full[w] >= 0: cnt[full[w]] += 1
                full[v] = int(np.argmin(cnt))
        cl = []
        for e in np.flatnonzero(full[U] == full[V]):
            y = nxt; nxt += 1; A.add_clause([-y, int(U[e]) + 1]); A.add_clause([-y, int(V[e]) + 1]); cl.append(y)
        A.add_clause(cl)

if __name__ == '__main__':
    name, smax = sys.argv[1], int(sys.argv[2]); P = points(name, smax)
    vals = sorted(set(((P * P).sum(1)).tolist()) - {0})   # lattice: squared distances = norms
    out = open('Borsuk/.scratch/d2_%s.jsonl' % name.replace('*', 's'), 'w')
    for s in vals:
        r = {'lattice': name, **solve(P, s)}; print(json.dumps(r), flush=True); out.write(json.dumps(r) + '\n'); out.flush()
