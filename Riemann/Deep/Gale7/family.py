#!/usr/bin/env python3
"""Exact check of the simplex-product family V_n for Gale's conjecture (Python standard library only).

The family.  R^n = R^m1 (+) R^m2, m1 + m2 = n.  u_0..u_m1 in R^m1 and v_0..v_m2 in R^m2 are the vertices of
regular simplices: |u_i| = |v_j| = 1, u_i.u_k = -1/m1 (i != k), v_j.v_l = -1/m2 (j != l).  For parameters
p = (rA, rB, a, b, c, d, e, f) the set V_n(p) has (m1+1) + (m2+1) + 3(m1+1)(m2+1) points:
    P_i = rA u_i,   Q_j = rB v_j,   A_ij = a u_i - b v_j,   B_ij = -c u_i + d v_j,   C_ij = e u_i + f v_j.
Every squared distance depends only on the two orbits and on which indices agree, so it is one of 36 explicit
rational functions of p (a "class"); with rational p every squared distance is an exact rational number.

Lemma (proved in NOTE.md, section 3; checked by exhaustive search below for small m1, m2).  Let H be the graph on
V_n whose edges are the pairs in the 17 classes
    PP(i!=k)  QQ(j!=l)  P_k A_ij(k!=i)  Q_l B_ij(l!=j)  A_ij A_kl(i!=k)  B_ij B_kl(j!=l)  P_k C_ij(k!=i)
    A_ij C_kl(i!=k)  Q_l C_ij(l!=j)  B_ij C_kl(j!=l)  C_ij C_kl(i!=k, j!=l)  P_i B_ij  Q_j A_ij
(AA, BB, AC, BC count as two classes each: same or different second index).  Then H has no proper colouring
with n + 1 colours.
Consequence.  If every pair in these 17 classes has squared distance >= t, then every partition of V_n into n+1
parts has a part containing two points at distance >= sqrt(t); so beta_{n+1}(V_n) >= sqrt(t) and
beta_{n+1}(V_n)/D(V_n) >= sqrt(t / D^2), D^2 = the largest of the 36 class values.  If t/D^2 > sigma_n^2, V_n is
a counterexample to Gale's conjecture in R^n (sigma_n^2 = (n+1)/(n+2) for even n, 1/2 + sqrt((n-1)/(n+3))/2 for
odd n; the ball satisfies beta_{n+1}(B^n) <= sigma_n by the regular-simplex partition).

Usage:
  python3 family.py                    check every row of family_params.txt (exact), print a summary
  python3 family.py --n 7              check one dimension, print all 36 classes
  python3 family.py --lemma            also re-check the colouring lemma by exhaustive search for small (m1, m2)
"""
import os, sys, time
from fractions import Fraction as F
from decimal import Decimal, getcontext

ORB = ['P', 'Q', 'A', 'B', 'C']
HAS_A = {'P': True, 'Q': False, 'A': True, 'B': True, 'C': True}
HAS_B = {'P': False, 'Q': True, 'A': True, 'B': True, 'C': True}
# the 17 edge classes used by the lemma: (orbit, orbit, relation of first indices, relation of second indices)
EDGES = [('P','P','n','-'), ('Q','Q','-','n'), ('P','A','n','-'), ('Q','B','-','n'), ('A','A','n','e'), ('A','A','n','n'),
         ('B','B','e','n'), ('B','B','n','n'), ('P','C','n','-'), ('A','C','n','e'), ('A','C','n','n'), ('Q','C','-','n'),
         ('B','C','e','n'), ('B','C','n','n'), ('C','C','n','n'), ('P','B','e','-'), ('Q','A','-','e')]

def coeffs(p):
    rA, rB, a, b, c, d, e, f = p
    return {'P': (rA, 0), 'Q': (0, rB), 'A': (a, -b), 'B': (-c, d), 'C': (e, f)}

def classes():
    out = []
    for x in range(5):
        for y in range(x, 5):
            X, Y = ORB[x], ORB[y]
            r1 = ['e', 'n'] if HAS_A[X] and HAS_A[Y] else ['-']
            r2 = ['e', 'n'] if HAS_B[X] and HAS_B[Y] else ['-']
            for s1 in r1:
                for s2 in r2:
                    if X == Y and s1 != 'n' and s2 != 'n': continue     # the same point
                    out.append((X, Y, s1, s2))
    return out

def part(x, y, rel, m):
    """|x w_i - y w_k|^2 for unit simplex vectors w (dot 1 if rel == 'e', -1/m if 'n'); '-' : one coefficient is 0"""
    if rel == '-': return x*x + y*y
    g = F(1) if rel == 'e' else F(-1, m)
    return x*x + y*y - 2*x*y*g

def class_d2(m1, m2, p, cl):
    X, Y, s1, s2 = cl
    cf = coeffs(p)
    (a1, b1), (a2, b2) = cf[X], cf[Y]
    return part(a1, a2, s1, m1) + part(b1, b2, s2, m2)

def sigma2_above(r, n):
    """exact test r > sigma_n^2"""
    if n % 2 == 0: return r > F(n + 1, n + 2)
    t = 2*r - 1
    return t > 0 and t*t > F(n - 1, n + 3)

def sigma_dec(n):
    getcontext().prec = 40
    if n % 2 == 0: return (Decimal(n + 1)/Decimal(n + 2)).sqrt()
    return (Decimal(1)/2 + (Decimal(n - 1)/Decimal(n + 3)).sqrt()/2).sqrt()

def check(n, m1, m2, p, verbose=False):
    assert m1 + m2 == n and m1 >= 1 and m2 >= 1
    cls = classes()
    assert len(cls) == 36 and all(e in cls for e in EDGES)
    vals = {cl: class_d2(m1, m2, p, cl) for cl in cls}
    assert all(v > 0 for v in vals.values()), 'two orbits share a point'
    D2 = max(vals.values()); t = min(vals[e] for e in EDGES)
    r = t / D2
    if verbose:
        for cl in sorted(cls, key=lambda c: -vals[c]):
            tag = 'edge' if cl in EDGES else ''
            print('  %s%s %s%s  d2/D2 = %.10f %s' % (cl[0], cl[1], cl[2], cl[3], float(vals[cl]/D2), tag))
    return r, sigma2_above(r, n)

def read_params(path):
    rows = []
    for line in open(path):
        line = line.split('#')[0].strip()
        if not line: continue
        w = line.split()
        rows.append((int(w[0]), int(w[1]), int(w[2]), tuple(F(x) for x in w[3:11])))
    return rows

# ---------------- independent exhaustive check of the lemma for small (m1, m2) ----------------
def lemma_graph(m1, m2):
    I, J = range(m1 + 1), range(m2 + 1)
    V = [('P', i, None) for i in I] + [('Q', None, j) for j in J] + [(t, i, j) for t in 'ABC' for i in I for j in J]
    idx = {v: k for k, v in enumerate(V)}
    adj = [set() for _ in V]
    def rel(x, y):
        return '-' if x is None or y is None else ('e' if x == y else 'n')
    for s in range(len(V)):
        for t in range(s + 1, len(V)):
            X, Y = V[s], V[t]
            if ORB.index(X[0]) > ORB.index(Y[0]): X, Y = Y, X
            if (X[0], Y[0], rel(X[1], Y[1]), rel(X[2], Y[2])) in EDGES:
                adj[s].add(t); adj[t].add(s)
    return len(V), adj

def colourable(n, adj, k):
    """exhaustive DSATUR backtracking with forward checking; colours of a greedy clique fixed (w.l.o.g.)."""
    order = sorted(range(n), key=lambda v: -len(adj[v]))
    cl = [order[0]]; cand = set(adj[order[0]])
    while cand and len(cl) < k:
        v = max(cand, key=lambda u: len(adj[u] & cand)); cl.append(v); cand &= adj[v]
    col = [-1]*n; dom = [(1 << k) - 1]*n
    def assign(v, c, trail):
        col[v] = c; trail.append((v, None))
        for w in adj[v]:
            if col[w] == -1 and dom[w] >> c & 1:
                trail.append((w, dom[w])); dom[w] &= ~(1 << c)
                if dom[w] == 0: return False
        return True
    def undo(trail, mark):
        while len(trail) > mark:
            v, d = trail.pop()
            if d is None: col[v] = -1
            else: dom[v] = d
    trail = []
    for c, v in enumerate(cl):
        if not (dom[v] >> c & 1) or not assign(v, c, trail): return False
    def rec(used):
        best, bs = -1, 99
        for v in range(n):
            if col[v] == -1:
                s = bin(dom[v]).count('1')
                if s < bs: best, bs = v, s
        if best == -1: return True
        for c in range(min(k, used + 1)):
            if dom[best] >> c & 1:
                mark = len(trail)
                if assign(best, c, trail) and rec(max(used, c + 1)): return True
                undo(trail, mark)
        return False
    return rec(len(cl))

def main():
    args = sys.argv[1:]
    here = os.path.dirname(os.path.abspath(__file__))
    path = os.path.join(here, 'family_params.txt')
    if '--lemma' in args:
        t0 = time.time()
        for m1, m2 in [(1,1), (1,2), (2,2), (1,3), (2,3), (3,3), (2,4), (3,4)]:
            nv, adj = lemma_graph(m1, m2)
            ok = not colourable(nv, adj, m1 + m2 + 1)
            print('lemma graph H(%d,%d): %d vertices, %d edges, not %d-colourable: %s' % (m1, m2, nv, sum(map(len, adj))//2, m1+m2+1, ok))
            assert ok
        print('lemma checks: %.1f s' % (time.time() - t0))
    rows = read_params(path)
    if '--n' in args:
        n = int(args[args.index('--n') + 1])
        rows = [r for r in rows if r[0] == n]
        for n, m1, m2, p in rows:
            r, ok = check(n, m1, m2, p, verbose=True)
            getcontext().prec = 40
            val = (Decimal(r.numerator)/Decimal(r.denominator)).sqrt()
            print('n = %d, split (%d,%d), %d points: beta_%d/D >= sqrt(%s/%s)' % (n, m1, m2, (m1+1)+(m2+1)+3*(m1+1)*(m2+1), n+1, r.numerator, r.denominator))
            print('   = %s...  > sigma_%d = %s...: %s' % (str(val)[:14], n, str(sigma_dec(n))[:14], ok))
        return 0
    bad = []; worst = None
    for n, m1, m2, p in rows:
        r, ok = check(n, m1, m2, p)
        if not ok: bad.append(n)
        getcontext().prec = 40
        exc = (Decimal(r.numerator)/Decimal(r.denominator)).sqrt() - sigma_dec(n)
        if worst is None or exc < worst[1]: worst = (n, exc)
    ns = [r[0] for r in rows]
    print('checked %d dimensions n = %d..%d (exact rational arithmetic)' % (len(rows), min(ns), max(ns)))
    print('beta_{n+1}/D > sigma_n for all of them: %s' % (not bad))
    print('smallest excess over sigma_n: %.4e at n = %d' % (float(worst[1]), worst[0]))
    return 0 if not bad else 1

if __name__ == '__main__':
    sys.exit(main())

