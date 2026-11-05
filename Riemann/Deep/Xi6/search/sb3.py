#!/usr/bin/env python3
"""sb3.py - SAT on plain P:m:3:f instance with symmetry breaking (colour order along W_0, orbit-shift canonical X_01[0], X_02[0]).
usage: sb3.py <inst tag e.g. f_P_12_3_0> <m> [solver]"""
import sys, os, pickle, time, json
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from sym import write_colouring, check_fast
tag, m = sys.argv[1], int(sys.argv[2]); solver = sys.argv[3] if len(sys.argv) > 3 else 'cadical153'
n, Cp, M = pickle.load(open(os.path.join(HERE, 'inst', tag + '.pkl'), 'rb'))
E = [list(map(int, l.split())) for l in open(os.path.join(HERE, 'inst', tag + '.hyp')).read().strip().split('\n')[1:]]
v = lambda k, c: 3 * k + c + 1
cl = []
for k in range(M):
    cl.append([v(k, 0), v(k, 1), v(k, 2)])
    for x, y in ((0, 1), (0, 2), (1, 2)): cl.append([-v(k, x), -v(k, y)])
for T in E:
    for c in range(3): cl.append([v(k, c) for k in T])
seq = [Cp[0][d] for d in range(1, m // 2 + 1)]      # W_0 classes d = 1..m/2
cl.append([v(seq[0], 0)])
for p in range(1, len(seq)):
    cl.append([-v(seq[p], 2)] + [v(seq[q], 1) for q in range(p)])
for j in (1, 2):
    k0 = Cp[0][j * m]
    for c in (1, 2):
        for d in range(1, m):
            for c2 in range(c):
                cl.append([-v(k0, c), -v(Cp[0][j * m + d], c2)])
from pysat.solvers import Solver
t0 = time.time()
print('built', len(cl), flush=True)
with Solver(name=solver, bootstrap_with=cl) as s:
    ok = s.solve()
    res = dict(inst=tag, sb='colour-order W0 + shift-canon X01,X02', solver=solver, res='SAT' if ok else 'UNSAT', t=round(time.time() - t0, 1))
    if ok:
        mod = set(x for x in s.get_model() if x > 0)
        col = [next(c for c in range(3) if v(k, c) in mod) for k in range(M)]
        Mx = [[-1 if a == b else col[Cp[a][b]] for b in range(n)] for a in range(n)]
        res['check'] = check_fast(Mx, 6)[0]
        write_colouring(os.path.join(HERE, 'sol', 'sb3_' + tag + '.txt'), Mx, note='sb3 ' + tag)
print(json.dumps(res), flush=True)
with open(os.path.join(HERE, 'results.jsonl'), 'a') as fh: fh.write(json.dumps(res) + '\n')

