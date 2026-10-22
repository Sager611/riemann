#!/usr/bin/env python3
"""Exact certificate check: beta_4(V)/diam(V) for a finite point set V in R^3 (Gale's conjecture, k = 4 parts).

Usage:  python3 verify.py [points.txt] [--cnf FILE]

points.txt : one point per line, three exact rationals "x y z" (integers or p/q); '#' lines are comments.
Only the Python standard library is used; all arithmetic is exact (integers / fractions).

What is checked:
  1. All squared distances are computed exactly (coordinates scaled to a common denominator);
     D^2 = max squared distance.
  2. For a threshold s (a squared distance), G(>= s) is the graph on V joining pairs at squared distance >= s.
     The script finds the critical value dc^2: the largest squared distance s such that G(>= s) is NOT
     4-colourable (binary search over the distinct squared distances).
  3. Non-4-colourability of G(>= dc^2) is proved by exhaustive backtracking: the four vertices of one K4 get
     colours 0,1,2,3 (w.l.o.g., colour symmetry), then DSATUR-style branching with forward checking; every
     branch ends in a domain wipe-out.  The node count is printed.  A second, independent exhaustive search
     (plain DSATUR backtracking, no propagation, a different K4) is run as a cross-check.
  4. An explicit proper 4-colouring of G(> dc^2) is found and checked edge by edge, so beta_4(V) = dc exactly.
  5. beta_4(V)/D = sqrt(dc^2/D^2) is printed as an exact reduced fraction under the root and in decimal, and
     compared EXACTLY with (a) the ball value sqrt((3+sqrt 3)/6) and (b) the previous record
     sqrt(39541040058709/50000010510015) = 0.8892809651... (101-point set of arXiv 2610.28828).
Optional: --cnf FILE writes the DIMACS CNF "G(>= dc^2) is 4-colourable" (variable 4*i+c+1 = vertex i gets
colour c; no symmetry breaking), for an independent SAT/DRAT check (e.g. cadical --proof + drat-trim).

Why this proves the bound: any partition of V into 4 parts is a 4-colouring of G(>= dc^2); since none is
proper, some part contains two points at distance >= dc, i.e. has diameter >= dc.  Hence beta_4(V) >= dc,
and the colouring of G(> dc^2) shows equality.
"""
import os, sys, time
from fractions import Fraction
from math import gcd, isqrt
from decimal import Decimal, getcontext

PREV = Fraction(39541040058709, 50000010510015)   # (beta_4/D)^2 of the 101-point set of arXiv 2610.28828

def read_points(path):
    pts = []
    for line in open(path):
        line = line.split('#')[0].strip()
        if not line: continue
        xs = line.split()
        assert len(xs) == 3, line
        pts.append(tuple(Fraction(s) for s in xs))
    assert len(set(pts)) == len(pts), 'repeated point'
    return pts

def integer_coords(pts):
    den = 1
    for p in pts:
        for c in p:
            den = den * c.denominator // gcd(den, c.denominator)
    return [tuple(int(c * den) for c in p) for p in pts], den

def find_k4(n, adj, skip=()):
    order = sorted(range(n), key=lambda v: -len(adj[v]))
    for a in order:
        for b in sorted(adj[a]):
            if b <= a: continue
            ab = adj[a] & adj[b]
            for c in sorted(ab):
                if c <= b: continue
                for d in sorted(ab & adj[c]):
                    if d > c and (a, b, c, d) not in skip:
                        return (a, b, c, d)
    return None

def colour_fc(n, adj, k4):
    """Backtracking with forward checking. Returns (colouring or None, nodes)."""
    dom = [15] * n; col = [-1] * n; nodes = [0]
    nbrs = [sorted(adj[v]) for v in range(n)]
    trail = []
    def assign(v, c):
        stack = [(v, c)]
        while stack:
            u, cu = stack.pop()
            if col[u] != -1:
                if col[u] != cu: return False
                continue
            if not (dom[u] >> cu) & 1: return False
            trail.append((u, dom[u], True)); col[u] = cu; dom[u] = 1 << cu
            bit = 1 << cu
            for w in nbrs[u]:
                if dom[w] & bit:
                    if col[w] != -1: return False
                    trail.append((w, dom[w], False)); dom[w] &= ~bit
                    if dom[w] == 0: return False
                    if dom[w] & (dom[w] - 1) == 0:
                        stack.append((w, dom[w].bit_length() - 1))
        return True
    def undo(mark):
        while len(trail) > mark:
            u, d, was_assign = trail.pop()
            dom[u] = d
            if was_assign: col[u] = -1
    def rec():
        nodes[0] += 1
        best, bs, bd = -1, 5, -1
        for v in range(n):
            if col[v] == -1:
                s = bin(dom[v]).count('1')
                if s < bs or (s == bs and len(nbrs[v]) > bd):
                    best, bs, bd = v, s, len(nbrs[v])
        if best == -1: return True
        for c in range(4):
            if (dom[best] >> c) & 1:
                mark = len(trail)
                if assign(best, c) and rec(): return True
                undo(mark)
        return False
    if k4 is not None:
        for i, v in enumerate(k4):
            if not assign(v, i): return None, 1
    ok = rec()
    return (list(col) if ok else None), nodes[0]

def colour_plain(n, adj, k4, limit=10**7):
    """Plain DSATUR backtracking (no propagation): colour = smallest-saturation-first, check only
    already-coloured neighbours. Returns (colouring or None, nodes)."""
    col = [-1] * n; nodes = [0]
    nbrs = [sorted(adj[v]) for v in range(n)]
    for i, v in enumerate(k4 or ()):
        col[v] = i
    def rec():
        nodes[0] += 1
        if nodes[0] > limit: raise RuntimeError('node limit')
        best, bsat, bdeg = -1, -1, -1
        for v in range(n):
            if col[v] == -1:
                sat = len({col[w] for w in nbrs[v] if col[w] != -1})
                if sat > bsat or (sat == bsat and len(nbrs[v]) > bdeg):
                    best, bsat, bdeg = v, sat, len(nbrs[v])
        if best == -1: return True
        used = {col[w] for w in nbrs[best] if col[w] != -1}
        for c in range(4):
            if c not in used:
                col[best] = c
                if rec(): return True
                col[best] = -1
        return False
    ok = rec()
    return (list(col) if ok else None), nodes[0]

def graph(n, d2, s, strict):
    adj = [set() for _ in range(n)]; m = 0
    for (i, j), d in d2.items():
        if d > s or (not strict and d == s):
            adj[i].add(j); adj[j].add(i); m += 1
    return adj, m

def main():
    args = sys.argv[1:]
    cnf = None
    if '--cnf' in args:
        k = args.index('--cnf'); cnf = args[k + 1]; del args[k:k + 2]
    path = args[0] if args else os.path.join(os.path.dirname(os.path.abspath(__file__)), 'points.txt')
    t0 = time.time()
    sys.setrecursionlimit(20000)
    pts = read_points(path); n = len(pts)
    P, den = integer_coords(pts)
    d2 = {}
    for i in range(n):
        for j in range(i + 1, n):
            d2[(i, j)] = (P[i][0]-P[j][0])**2 + (P[i][1]-P[j][1])**2 + (P[i][2]-P[j][2])**2
    D2 = max(d2.values())
    dist = sorted(set(d2.values()), reverse=True)
    def noncol(idx):
        adj, _ = graph(n, d2, dist[idx], strict=False)
        return colour_fc(n, adj, find_k4(n, adj))[0] is None
    hi = next(i for i, d in enumerate(dist) if 10000 * d < 7225 * D2)      # first value below (0.85 D)^2
    assert noncol(hi)
    lo = 0
    while lo < hi:
        mid = (lo + hi) // 2
        if noncol(mid): hi = mid
        else: lo = mid + 1
    dc2 = dist[lo]
    adj, m = graph(n, d2, dc2, strict=False)
    k4 = find_k4(n, adj)
    c1, nodes1 = colour_fc(n, adj, k4)
    assert c1 is None
    k4b = find_k4(n, adj, skip={k4})
    c1b, nodes1b = colour_plain(n, adj, k4b)
    assert c1b is None
    adj2, m2 = graph(n, d2, dc2, strict=True)
    c2, _ = colour_fc(n, adj2, find_k4(n, adj2))
    assert c2 is not None and all(c2[v] != c2[w] for v in range(n) for w in adj2[v])
    r = Fraction(dc2, D2)
    getcontext().prec = 50
    val = (Decimal(r.numerator) / Decimal(r.denominator)).sqrt()
    # exact comparisons: r > PREV ; r > (3+sqrt3)/6  <=>  6r-3 > 0 and (6r-3)^2 > 3
    above_ball = (6 * r - 3 > 0) and (6 * r - 3) ** 2 > 3
    print('file: %s   points: %d   common denominator: %d' % (path, n, den))
    print('D^2           = %d / %d^2' % (D2, den))
    print('dc^2          = %d / %d^2   (critical squared distance)' % (dc2, den))
    print('G(>= dc): %d edges; NOT 4-colourable: forward-checking search %d nodes (K4 %s fixed), '
          'plain DSATUR search %d nodes (K4 %s fixed)' % (m, nodes1, k4, nodes1b, k4b))
    print('G(>  dc): %d edges; explicit proper 4-colouring found and checked => beta_4(V) = dc exactly' % m2)
    print('beta_4(V)/D(V) = sqrt(%d/%d)' % (r.numerator, r.denominator))
    print('               = %s...' % str(val)[:16])
    print('> ball value sqrt((3+sqrt3)/6) = 0.8880738...: %s' % above_ball)
    print('> previous record sqrt(39541040058709/50000010510015) = 0.8892809651...: %s (gain %.3e)'
          % (r > PREV, float(val) - float(Decimal(PREV.numerator) / Decimal(PREV.denominator)) ** 0.5))
    if cnf:
        with open(cnf, 'w') as f:
            f.write('p cnf %d %d\n' % (4 * n, n + 4 * m))
            for i in range(n):
                f.write('%d %d %d %d 0\n' % tuple(4 * i + c + 1 for c in range(4)))
            for (i, j), d in sorted(d2.items()):
                if d >= dc2:
                    for c in range(4):
                        f.write('-%d -%d 0\n' % (4 * i + c + 1, 4 * j + c + 1))
        print('CNF written to %s' % cnf)
    print('time %.1f s' % (time.time() - t0))
    if r > PREV and above_ball:
        print('RESULT: certified beta_4(V)/D(V) = sqrt(%d/%d) = %s... > previous record 0.8892809651...'
              % (r.numerator, r.denominator, str(val)[:12]))
    else:
        print('RESULT: no improvement over the previous record')
    return 0 if r > PREV else 1

if __name__ == '__main__':
    sys.exit(main())
