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
sys.path.insert(0, 'Borsuk/code'); import cegar as cg

def run(name, X, key, fv, tl=600):
    n = len(X); DA, DB = d1.sqd_block(X, 0, 1); kk = d1.keys_of(DA[0], DB[0])
    sv = fv(key); vals = np.array([fv(int(t)) for t in kk]); tol = 1e-9 * max(1.0, abs(sv))
    L = np.flatnonzero(vals <= sv + tol); Y = X[L]; m = len(L)
    EQ = np.zeros((m, m), bool); BAD = np.zeros((m, m), bool)
    for i0 in range(0, m, 128):
        DA, DB = d1.sqd_block(Y, i0, min(m, i0 + 128)); K2 = d1.keys_of(DA, DB)
        F = np.vectorize(lambda t: fv(int(t)))(K2)
        EQ[i0:i0 + 128] = K2 == key; BAD[i0:i0 + 128] = F > sv + tol
        assert not np.any((np.abs(F - sv) <= tol) & (K2 != key)), 'float tie'
    U, V = [a.astype(np.int32) for a in np.nonzero(np.triu(EQ, 1))]
    root = int(np.flatnonzero(L == 0)[0])
    res = cg.cegar(m, U, V, BAD, tl=tl, root=root)
    out = {'set': name, 's': round(sv, 6), 'L': m, 'result': res[0], 'info': res[1] if res[0] != 'CANDIDATE' else None}
    if res[0] == 'CANDIDATE': out['S'] = [int(L[v]) for v in res[1]]
    return out

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
                print(json.dumps(run(name, X, key, fv)), flush=True)
