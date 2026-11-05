#!/usr/bin/env python3
"""fullsat.py - SAT for a 3-colouring of K_n with every K_k 3-coloured, no symmetry assumed.
Symmetry breaking: colours of edges (0,1..n-1) non-decreasing; colour 0 is the most frequent at vertex 0
and colour 1 at least as frequent as colour 2 at vertex 0 (encoded via sorted prefix: c(0,j) pattern 0..0 1..1 2..2
with #0 >= #1 >= #2).  usage: fullsat.py n k [solver]"""
import sys, os, itertools, time, json
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from sym import write_colouring, check_fast
n, k = int(sys.argv[1]), int(sys.argv[2]); solver = sys.argv[3] if len(sys.argv) > 3 else 'kissat404'
E = {}
for a in range(n):
    for b in range(a + 1, n): E[(a, b)] = len(E)
v = lambda e, c: 3 * e + c + 1
cl = []
for e in range(len(E)):
    cl.append([v(e, 0), v(e, 1), v(e, 2)])
    for x, y in ((0, 1), (0, 2), (1, 2)): cl.append([-v(e, x), -v(e, y)])
for S in itertools.combinations(range(n), k):
    es = [E[p] for p in itertools.combinations(S, 2)]
    for c in range(3): cl.append([v(e, c) for e in es])
# vertex-0 row sorted: colour(0,j) <= colour(0,j+1)
row = [E[(0, j)] for j in range(1, n)]
for i in range(len(row) - 1):
    a, b = row[i], row[i + 1]
    # forbid a=1,b=0 ; a=2,b=0 ; a=2,b=1
    cl.append([-v(a, 1), -v(b, 0)]); cl.append([-v(a, 2), -v(b, 0)]); cl.append([-v(a, 2), -v(b, 1)])
# counts: with sorted row, #0 >= #1 >= #2 means: let t = n-1. position-based: #0>=#1 iff row[i] for i < ceil(t/2)... encode via
# for each i: if row[i]==1 then row[2*i - ...]; simpler: #0 >= #1  <=>  for all j: (row[j]==1) -> (row[j - #?])... use direct:
# #0 >= #1: for every i, if row[i] is the first 1 at index i (i.e. #0 = i), then row[2i] is not 1 (so #1 <= i).
t = n - 1
for i in range(t):
    # first1 at i : row[i]==1 and (i==0 or row[i-1]==0)
    if 2 * i < t:
        lits = [-v(row[i], 1)] + ([-v(row[i - 1], 0)] if i > 0 else []) + [-v(row[2 * i], 1)]
        cl.append(lits)
    # if i == 0 then row[0]==1 means #0 = 0 -> must have #1 = 0, contradiction unless... forbid row[0]==1 with any 1 after: covered by 2*i=0 -> lits [-row0=1, -row0=1] ok
# #1 >= #2: first 2 at index j with first 1 at i: #1 = j - i, #2 = t - j  => need t - j <= j - i  -> 2j - i >= t
for i in range(t):
    for j in range(i + 1, t):
        if 2 * j - i < t:
            lits = [-v(row[i], 1)] + ([-v(row[i - 1], 0)] if i > 0 else []) + [-v(row[j], 2), -v(row[j - 1], 1)]
            cl.append(lits)
# if no 1 at all but 2s exist: #1=0 < #2 -> forbid: row[j]==2 and row[j-1]==0
for j in range(1, t):
    cl.append([-v(row[j], 2), -v(row[j - 1], 0)])
cl.append([-v(row[0], 2)])
print(json.dumps(dict(n=n, k=k, vars=3 * len(E), clauses=len(cl), solver=solver)), flush=True)
from pysat.solvers import Solver
t0 = time.time()
with Solver(name=solver, bootstrap_with=cl) as s:
    ok = s.solve()
    res = dict(n=n, k=k, res='SAT' if ok else 'UNSAT', t=round(time.time() - t0, 1), solver=solver)
    if ok:
        mod = set(x for x in s.get_model() if x > 0)
        M = [[-1] * n for _ in range(n)]
        for (a, b), e in E.items():
            c = next(c for c in range(3) if v(e, c) in mod); M[a][b] = M[b][a] = c
        res['check'] = check_fast(M, k)[0]
        write_colouring(os.path.join(HERE, 'sol', f'full_n{n}_k{k}.txt'), M, note=f'full SAT n={n} k={k}')
print(json.dumps(res), flush=True)
with open(os.path.join(HERE, 'results.jsonl'), 'a') as fh: fh.write(json.dumps(res) + '\n')

