#!/usr/bin/env python3
"""D1 refinement (v2) for thresholds whose local graph H_s[L_s] was not 5-coloured by tabu search.

CEGAR over subsets S of L_s containing x0 whose pairs are all at squared distance <= s: a SAT solver proposes S
(extended to a maximal compatible set); tabu search (SAT fallback) 5-colours H_s[S]; the colouring, extended
greedily to L_s, blocks every S on which it is proper.  UNSAT of the proposer = no subset of P with diameter
sqrt(s) and a non-5-colourable diameter graph.  Exact keys throughout; colourings checked edge by edge.
"""
import sys, json, time
import numpy as np
from pysat.solvers import Solver
sys.path.insert(0, 'Borsuk/code'); sys.path.insert(0, 'Borsuk/.scratch')
import d1
from tabu import colour_tabu

def cegar(name, X, key, fv, tl=1800):
    n = len(X); DA, DB = d1.sqd_block(X, 0, 1); kk = d1.keys_of(DA[0], DB[0])
    sv = fv(key); vals = np.array([fv(int(t)) for t in kk]); tol = 1e-9 * max(1.0, abs(sv))
    L = np.flatnonzero(vals <= sv + tol); Y = X[L]; m = len(L)
    DA, DB = d1.sqd_block(Y, 0, m); K2 = d1.keys_of(DA, DB)
    F = np.vectorize(lambda t: fv(int(t)))(K2)
    E = np.triu(K2 == key, 1); U, V = [a.astype(np.int32) for a in np.nonzero(E)]
    bad = np.triu(F > sv + tol, 1)
    assert not np.any(np.triu(np.abs(F - sv) <= tol, 1) & ~E), 'float tie'
    x0 = int(np.flatnonzero(L == 0)[0])
    A = Solver(name='g4'); A.add_clause([x0 + 1])
    for a, b in zip(*np.nonzero(bad)): A.add_clause([-(int(a) + 1), -(int(b) + 1)])
    compat = ~(bad | bad.T); adj = [[] for _ in range(m)]
    for a, b in zip(U.tolist(), V.tolist()): adj[a].append(b); adj[b].append(a)
    nxt = m + 1; t0 = time.time(); it = 0
    while True:
        it += 1
        if time.time() - t0 > tl: return {'set': name, 's': round(sv, 6), 'L': m, 'result': 'TIMEOUT', 'iters': it}
        if not A.solve(): return {'set': name, 's': round(sv, 6), 'L': m, 'result': 'UNSAT (dead)', 'iters': it, 'sec': round(time.time() - t0, 1)}
        mod = A.get_model(); inS = np.zeros(m, bool); inS[[v for v in range(m) if mod[v] > 0]] = True
        for v in np.random.default_rng(it).permutation(m):
            if not inS[v] and compat[v, inS].all(): inS[v] = True
        S = np.flatnonzero(inS); idx = -np.ones(m, int); idx[S] = np.arange(len(S))
        sel = inS[U] & inS[V]; Us, Vs = idx[U[sel]].astype(np.int32), idx[V[sel]].astype(np.int32)
        col = None
        if len(Us) == 0: col = np.zeros(len(S), int)
        else:
            for sd in range(1, 4):
                c_, conf = colour_tabu(len(S), Us, Vs, 5, iters=200000 * sd, seed=sd)
                if conf == 0: col = c_; break
        if col is None:
            st, col = d1.colour_check(len(S), Us, Vs, 5, budget=20_000_000)
            if col is None:
                return {'set': name, 's': round(sv, 6), 'L': m, 'result': 'CANDIDATE ' + st, 'S': L[S].tolist()}
        col = np.asarray(col); assert (col[Us] != col[Vs]).all()
        full = -np.ones(m, int); full[S] = col
        for v in range(m):
            if full[v] < 0:
                cnt = np.zeros(5, int)
                for w in adj[v]:
                    if full[w] >= 0: cnt[full[w]] += 1
                full[v] = int(np.argmin(cnt))
        mono = np.flatnonzero(full[U] == full[V]); cl = []
        for e in mono:
            y = nxt; nxt += 1; A.add_clause([-y, int(U[e]) + 1]); A.add_clause([-y, int(V[e]) + 1]); cl.append(y)
        A.add_clause(cl)

if __name__ == '__main__':
    src = sys.argv[1] if len(sys.argv) > 1 else 'Borsuk/.scratch/d1v2_results.jsonl'
    todo = {}
    for line in open(src):
        r = json.loads(line)
        for b in r['bad']: todo.setdefault(r['set'], []).append(b['key'])
    for g in sorted({s.split(':')[0] for s in todo}):
        for name, field, pts, vt in d1.sets_for(g):
            if name not in todo: continue
            d1.set_field(field); X = d1.to_int(pts); wf = d1.WF[field]
            fv = lambda key, wf=wf: (lambda ab: ab[0] + ab[1] * wf)(d1.unkey(key))
            for key in todo[name]:
                print(json.dumps(cegar(name, X, key, fv)), flush=True)
