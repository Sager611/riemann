#!/usr/bin/env python3
"""cube34.py - exhaust P:11:3:1 case (fixed point colours (0,1,2)) by cubes on (W0 rep, W1, W2)."""
import sys, os, pickle, itertools, time, json
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from sym import write_colouring, check_fast
m = 11
n, Cp, M = pickle.load(open(os.path.join(HERE, 'inst', 'f_P_11_3_1.pkl'), 'rb'))
lines = open(os.path.join(HERE, 'inst', 'f_P_11_3_1.hyp')).read().strip().split('\n')
E = [list(map(int, l.split())) for l in lines[1:]]
v = lambda k, c: 3 * k + c + 1
cl = []
for k in range(M):
    cl.append([v(k, 0), v(k, 1), v(k, 2)])
    for x, y in ((0, 1), (0, 2), (1, 2)): cl.append([-v(k, x), -v(k, y)])
for T in E:
    for c in range(3): cl.append([v(k, c) for k in T])
fx = [Cp[33][i * m] for i in range(3)]
for i in range(3): cl.append([v(fx[i], i)])
W = lambda i, d: Cp[i * m][i * m + d]          # d = 1..5
X = lambda i, j, d: Cp[i * m][j * m + d]       # i<j, d = 0..10
# shift canonicalisation: X_01[0] and X_02[0] hold the minimum colour present
for (i, j) in ((0, 1), (0, 2)):
    k0 = X(i, j, 0)
    for c in (1, 2):
        for d in range(1, m):
            for c2 in range(c):
                cl.append([-v(k0, c), -v(X(i, j, d), c2)])
# valid W lists (computed independently)
def has_clique(adj, nn, k):
    def ext(cand, size):
        if size == k: return True
        while cand:
            u = cand.bit_length() - 1; cand &= ~(1 << u)
            if size + 1 + bin(cand & adj[u]).count('1') >= k and ext(cand & adj[u], size + 1): return True
            if size + bin(cand).count('1') < k: return False
        return False
    return ext((1 << nn) - 1, 0)
def circ_adj(Wt, cols):
    adj = [0] * m
    for x in range(m):
        for y in range(m):
            if x != y:
                d = (y - x) % m; d = min(d, m - d)
                if Wt[d - 1] in cols: adj[x] |= 1 << y
    return adj
valid = {0: [], 1: [], 2: []}
for Wt in itertools.product(range(3), repeat=5):
    for i in range(3):
        oth = [c for c in range(3) if c != i]
        if any(has_clique(circ_adj(Wt, {i, x}), m, 5) for x in oth): continue
        if has_clique(circ_adj(Wt, set(oth)), m, 6): continue
        valid[i].append(Wt)
def mult(Wt, a): return tuple(Wt[min((a * d) % m, m - (a * d) % m) - 1] for d in range(1, 6))
reps0 = sorted({min(mult(Wt, a) for a in range(1, m)) for Wt in valid[0]})
assump = lambda i, Wt: [v(W(i, d), Wt[d - 1]) for d in range(1, 6)]
from pysat.solvers import Solver
s = Solver(name='cadical153', bootstrap_with=cl)
t0 = time.time()
p01 = {(a, b) for a in reps0 for b in valid[1] if s.solve(assumptions=assump(0, a) + assump(1, b))}
p02 = {(a, c) for a in reps0 for c in valid[2] if s.solve(assumptions=assump(0, a) + assump(2, c))}
print('pairs01', len(p01), 'of', len(reps0) * 30, 'pairs02', len(p02), 't', round(time.time() - t0, 1), flush=True)
p12 = {}
cnt = 0; found = None
for a in reps0:
    B = [b for b in valid[1] if (a, b) in p01]; Cc = [c for c in valid[2] if (a, c) in p02]
    for b in B:
        for c in Cc:
            cnt += 1
            if s.solve(assumptions=assump(0, a) + assump(1, b) + assump(2, c)):
                found = (a, b, c); break
        if found: break
    print('W0', a, 'B', len(B), 'C', len(Cc), 'cubes so far', cnt, 't', round(time.time() - t0, 1), flush=True)
    if found: break
res = dict(case='P:11:3:1 fixed (0,1,2)', cubes=cnt, res='SAT' if found else 'UNSAT', t=round(time.time() - t0, 1))
if found:
    mod = set(x for x in s.get_model() if x > 0)
    col = [next(c for c in range(3) if v(k, c) in mod) for k in range(M)]
    Mx = [[-1 if a_ == b_ else col[Cp[a_][b_]] for b_ in range(n)] for a_ in range(n)]
    res['check'] = check_fast(Mx, 6)[0]
    write_colouring(os.path.join(HERE, 'sol', 'cube34.txt'), Mx, note='P:11:3:1 case (0,1,2)')
print(json.dumps(res), flush=True)
with open(os.path.join(HERE, 'results.jsonl'), 'a') as fh: fh.write(json.dumps(res) + '\n')

