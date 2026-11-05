#!/usr/bin/env python3
"""ext.py - try to extend a valid colouring of K_n0 (file) by t new vertices via SAT.
usage: ext.py <colouring file> t [out file]"""
import sys, os, itertools, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sym import write_colouring, check_fast
def read(fn):
    rows = [l.strip() for l in open(fn) if l.strip() and not l.startswith('#')]
    n = int(rows[0]); rows = rows[1:1 + n]
    return n, [[-1 if i == j else int(rows[i][j]) for j in range(n)] for i in range(n)]
def extend(M, t, solver='cadical153', extra_clauses=None, return_all=False):
    from pysat.solvers import Solver
    n0 = len(M); N = n0 + t
    # edge variables for new edges: (a,b) with b >= n0
    var = {}
    def ev(a, b, c):
        if a > b: a, b = b, a
        key = (a, b)
        if key not in var: var[key] = 3 * len(var) + 1
        return var[key] + c
    newedges = [(a, b) for b in range(n0, N) for a in range(b)]
    for (a, b) in newedges: ev(a, b, 0)
    cl = []
    for (a, b) in newedges:
        cl.append([ev(a, b, 0), ev(a, b, 1), ev(a, b, 2)])
        for x, y in ((0, 1), (0, 2), (1, 2)): cl.append([-ev(a, b, x), -ev(a, b, y)])
    news = list(range(n0, N))
    # precompute colour masks of base subsets incrementally
    nb = 0
    for j in range(1, min(t, 6) + 1):
        s = 6 - j
        for J in itertools.combinations(news, j):
            Jedges = list(itertools.combinations(J, 2))
            for S in itertools.combinations(range(n0), s):
                m = 0
                for x, y in itertools.combinations(S, 2): m |= 1 << M[x][y]
                if m == 7: continue
                for c in range(3):
                    if m >> c & 1: continue
                    lits = [ev(a, b, c) for a in S for b in J] + [ev(a, b, c) for a, b in Jedges]
                    cl.append(lits); nb += 1
    if extra_clauses: cl += extra_clauses(ev)
    with Solver(name=solver, bootstrap_with=cl) as sv:
        ok = sv.solve()
        if not ok: return None, nb
        model = set(x for x in sv.get_model() if x > 0)
        M2 = [[-1] * N for _ in range(N)]
        for a in range(n0):
            for b in range(n0): M2[a][b] = M[a][b]
        for (a, b) in newedges:
            c = next(c for c in range(3) if ev(a, b, c) in model)
            M2[a][b] = M2[b][a] = c
        return M2, nb
if __name__ == '__main__':
    fn = sys.argv[1]; t = int(sys.argv[2]); out = sys.argv[3] if len(sys.argv) > 3 else None
    n0, M = read(fn)
    t0 = time.time()
    M2, nb = extend(M, t)
    print(f"base n0={n0} t={t} clauses={nb} result={'SAT' if M2 else 'UNSAT'} time={time.time()-t0:.1f}s")
    if M2:
        print('check', check_fast(M2, 6))
        if out: write_colouring(out, M2, note=f'extension of {os.path.basename(fn)} by {t}')

