#!/usr/bin/env python3
"""D1 refinement: thresholds s where the global graph H_s(P) was not shown 5-colourable.

P is vertex-transitive (an orbit of a reflection group), so a subset S with diam^2 S = s and a
non-5-colourable diameter graph may be moved by an isometry of P to contain the base point x0; then S lies in
L = {p : d^2(x0, p) <= s}.  Step 1: try to 5-colour H_s[L].  Step 2 (if that fails): CEGAR over subsets S of L
containing x0 with all pairs at squared distance <= s: a SAT solver proposes a maximal such S, a second solver
5-colours H_s[S]; each colouring (extended greedily to L) blocks every S on which it is proper.  The loop ends
with UNSAT (no counterexample at s) or with a non-5-colourable S (a counterexample candidate, printed).
Usage: python d1_local.py   (reads Borsuk/.scratch/d1_results.jsonl)
"""
import sys, json, time
import numpy as np
from pysat.solvers import Solver
sys.path.insert(0, 'Borsuk/code')
import d1

def local(name, X, key, k=5):
    n = len(X); DA, DB = d1.sqd_block(X, 0, 1); kk = d1.keys_of(DA[0], DB[0])
    wf = d1.WF[d1.FIELD]; fl = lambda key_: (lambda ab: ab[0] + ab[1] * wf)(d1.unkey(key_))
    sv = fl(key)
    vals = np.array([fl(int(t)) for t in kk])
    L = [0] + [j for j in range(1, n) if vals[j] <= sv + 1e-9 * max(1.0, abs(sv))]
    # exactness: no value within the tolerance band other than s itself
    near = [j for j in L if abs(vals[j] - sv) <= 1e-9 * max(1.0, abs(sv)) and int(kk[j]) != key]
    assert not near, 'float tie'
    Y = X[L]; m = len(L)
    DA, DB = d1.sqd_block(Y, 0, m); K2 = d1.keys_of(DA, DB)
    F = np.vectorize(lambda t: fl(int(t)))(K2) if m <= 400 else None
    ii, jj = np.nonzero(np.triu(K2 == key, 1)); U, V = ii.astype(np.int32), jj.astype(np.int32)
    st, col = d1.colour_check(m, U, V, k, budget=2_000_000)
    if col is not None:
        assert d1.verify(U, V, col, k); return {'set': name, 'L': m, 'edges': len(U), 'local': st}
    # CEGAR
    if F is None:
        F = np.zeros((m, m))
        for a in range(m):
            F[a] = [fl(int(t)) for t in K2[a]]
    bad = np.triu(F > sv + 1e-9 * max(1.0, abs(sv)), 1)
    A = Solver(name='g4'); A.add_clause([1])
    for a, b in zip(*np.nonzero(bad)): A.add_clause([-(int(a) + 1), -(int(b) + 1)])
    adj = [[] for _ in range(m)]
    for a, b in zip(U.tolist(), V.tolist()): adj[a].append(b); adj[b].append(a)
    compat = ~(bad | bad.T)
    nxt = m + 1; it = 0; t0 = time.time()
    while True:
        it += 1
        if not A.solve(): A.delete(); return {'set': name, 'L': m, 'edges': len(U), 'local': st, 'cegar': 'UNSAT', 'iters': it, 'sec': round(time.time() - t0, 1)}
        mod = A.get_model(); S = [v for v in range(m) if mod[v] > 0]
        inS = np.zeros(m, bool); inS[S] = True
        for v in range(m):           # extend to a maximal compatible set
            if not inS[v] and compat[v, inS].all(): inS[v] = True
        S = np.flatnonzero(inS).tolist(); idx = {v: i for i, v in enumerate(S)}
        es = [(a, b) for a, b in zip(U.tolist(), V.tolist()) if a in idx and b in idx]
        Us = np.array([idx[a] for a, b in es], np.int32); Vs = np.array([idx[b] for a, b in es], np.int32)
        st2, c2 = d1.colour_check(len(S), Us, Vs, k, budget=5_000_000)
        if c2 is None:
            return {'set': name, 'L': m, 'local': st, 'cegar': 'CANDIDATE ' + st2, 'S': [L[v] for v in S], 'iters': it}
        col = [-1] * m
        for v in S: col[v] = c2[idx[v]]
        for v in range(m):
            if col[v] < 0:
                cnt = [0] * k
                for w in adj[v]:
                    if col[w] >= 0: cnt[col[w]] += 1
                col[v] = int(np.argmin(cnt))
        mono = [(a, b) for a, b in zip(U.tolist(), V.tolist()) if col[a] == col[b]]
        cl = []
        for a, b in mono:
            y = nxt; nxt += 1; A.add_clause([-y, a + 1]); A.add_clause([-y, b + 1]); cl.append(y)
        A.add_clause(cl)

if __name__ == '__main__':
    todo = {}
    for line in open('Borsuk/.scratch/d1_results.jsonl'):
        r = json.loads(line)
        for b in r['bad']: todo.setdefault(r['set'], []).append(b['key'])
    groups = sorted({s.split(':')[0] for s in todo})
    for g in groups:
        for name, field, pts, vt in d1.sets_for(g):
            if name not in todo: continue
            assert vt
            d1.set_field(field); X = d1.to_int(pts)
            for key in todo[name]:
                print(json.dumps(local(name, X, key)), flush=True)
