#!/usr/bin/env python3
# controls for xi_check.py + structural check of K33 (Z_11 x 3 description in Riemann/Deep/Xi6/REPORT.md section 2)
import sys, random, subprocess, os
sys.path.insert(0, os.path.dirname(__file__))
import xi_check as X
def write(fn, n, col):
    with open(fn, 'w') as f:
        f.write('%d\n' % n)
        for i in range(n): f.write(''.join('.' if i == j else str(col[i][j]) for j in range(n)) + '\n')
def check(n, col, k):
    ok = True; ws = []
    for a, b in [(0, 1), (0, 2), (1, 2)]:
        adj = [sum(1 << j for j in range(n) if j != i and col[i][j] in (a, b)) for i in range(n)]
        w = X.clique_number(adj, n); ws.append(w); ok &= w <= k - 1
    return ok, ws
# (a) structure: rebuild K33 from W/X tables
W = [[0,0,1,2,0],[2,2,1,2,0],[1,0,0,2,0]]
Xt = {(0,1):'10020100122', (0,2):'12111202220', (1,2):'20210111020'}
n = 33; col = [[None]*n for _ in range(n)]
for i in range(3):
    for x in range(11):
        for j in range(3):
            for y in range(11):
                u, v = 11*i + x, 11*j + y
                if u == v: continue
                d = (y - x) % 11
                if i == j: c = W[i][min(d, 11 - d) - 1]
                elif i < j: c = int(Xt[(i, j)][d])
                else: c = int(Xt[(j, i)][(x - y) % 11])
                col[u][v] = c
n2, col2 = X.parse(sys.argv[1])
print('structure (W/X tables reproduce the file):', all(col[u][v] == col2[u][v] for u in range(n) for v in range(n) if u != v))
# Z_11 shift invariance
sh = lambda u: 11*(u // 11) + (u % 11 + 1) % 11
print('Z_11 shift invariant:', all(col2[sh(u)][sh(v)] == col2[u][v] for u in range(n) for v in range(n) if u != v))
# (b) control: cubic-residue colouring of K31 (arXiv 2611.03878 Thm 3.2) must FAIL for k=6, PASS for k=7
H = {1,30,2,29,4,27,8,23,15,16}
cls = lambda d: 0 if d in H else (1 if d in {(3*h) % 31 for h in H} else 2)
assert len({(3*h) % 31 for h in H} | H | {(9*h) % 31 for h in H}) == 30
c31 = [[None if i == j else cls((i - j) % 31) for j in range(31)] for i in range(31)]
print('K31 cubic residues, k=6:', check(31, c31, 6), ' k=7:', check(31, c31, 7))
# (c) control: random colourings of K33 must fail
random.seed(1); fails = 0
for t in range(20):
    r = [[None]*33 for _ in range(33)]
    for i in range(33):
        for j in range(i+1, 33): r[i][j] = r[j][i] = random.randrange(3)
    fails += not check(33, r, 6)[0]
print('random K33 colourings rejected: %d/20' % fails)
# (d) single-edge perturbations of K33: how many single recolourings keep it valid?
keep = 0; tot = 0
for u in range(33):
    for v in range(u+1, 33):
        for c in range(3):
            if c == col2[u][v]: continue
            tot += 1
            old = col2[u][v]; col2[u][v] = col2[v][u] = c
            keep += check(33, col2, 6)[0]
            col2[u][v] = col2[v][u] = old
print('single-edge recolourings that stay valid: %d of %d' % (keep, tot))
