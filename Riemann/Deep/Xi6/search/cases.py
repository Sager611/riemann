#!/usr/bin/env python3
"""cases.py - SAT on a dumped plain instance (inst/f_<tag>.hyp/.pkl) with extra unit clauses.
usage: cases.py <tag> <case-json: {"vertexpairs":[[a,b,colour],...]}> [nosym]"""
import sys, os, json, pickle, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from sym import write_colouring, check_fast
tag, case = sys.argv[1], json.loads(sys.argv[2])
hp = os.path.join(HERE, 'inst', 'f_' + tag)
n, Cp, m = pickle.load(open(hp + '.pkl', 'rb'))
lines = open(hp + '.hyp').read().strip().split('\n')
E = [list(map(int, l.split())) for l in lines[1:]]
from pysat.solvers import Solver
v = lambda k, c: 3 * k + c + 1
cl = []
for k in range(m):
    cl.append([v(k, 0), v(k, 1), v(k, 2)])
    for x, y in ((0, 1), (0, 2), (1, 2)): cl.append([-v(k, x), -v(k, y)])
for T in E:
    for c in range(3): cl.append([v(k, c) for k in T])
for a, b, c in case.get('vertexpairs', []):
    cl.append([v(Cp[a][b], c)])
t0 = time.time()
with Solver(name=os.environ.get('XI6_SOLVER','cadical153'), bootstrap_with=cl) as s:
    ok = s.solve()
    res = 'SAT' if ok else 'UNSAT'
    out = dict(tag=tag, case=case, res=res, t=round(time.time() - t0, 1))
    if ok:
        mod = set(x for x in s.get_model() if x > 0)
        col = [next(c for c in range(3) if v(k, c) in mod) for k in range(m)]
        M = [[-1 if a == b else col[Cp[a][b]] for b in range(n)] for a in range(n)]
        out['check'] = check_fast(M, 6)[0]
        write_colouring(os.path.join(HERE, 'sol', f'case_{tag}_{int(time.time())}.txt'), M, note=f'{tag} case {case}')
print(json.dumps(out), flush=True)
with open(os.path.join(HERE, 'results.jsonl'), 'a') as fh: fh.write(json.dumps(out) + '\n')

