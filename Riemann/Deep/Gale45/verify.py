#!/usr/bin/env python3
"""Exact certificate check: beta_k(V)/diam(V) with k = d+1 parts for a finite point set V in R^d
(Gale's conjecture; d = 4, 5 here, any d works).

Usage:  python3 verify.py POINTS [POINTS ...] [--cnf FILE]

POINTS : one point per line, d exact rationals (integers, p/q or finite decimals), separated by spaces or
commas; '#' starts a comment.  Only the Python standard library is used; all arithmetic is exact.

What is checked (for each file):
  1. All squared distances are computed exactly as integers (coordinates scaled to a common denominator);
     D^2 = max squared distance.
  2. For a squared distance s, G(>= s) is the graph on V joining pairs at squared distance >= s.  The
     critical value dc^2 = the largest squared distance s with G(>= s) NOT k-colourable is found by binary
     search over the distinct squared distances.
  3. G(>= dc^2) is proved not k-colourable by exhaustive backtracking: one k-clique gets colours 0..k-1
     (w.l.o.g., colour symmetry), then DSATUR branching with forward checking; every branch ends in a
     domain wipe-out.  A second, independent exhaustive search (plain DSATUR backtracking without
     propagation, a different k-clique fixed) is run as a cross-check.
  4. An explicit proper k-colouring of G(> dc^2) is found and checked edge by edge: beta_k(V) = dc exactly.
  5. beta_k(V)/D = sqrt(dc^2/D^2) is compared EXACTLY with the ball value sigma_d of the regular-simplex
     partition (sigma_d^2 = (d+1)/(d+2) for even d, 1/2 + sqrt((d-1)/(d+3))/2 for odd d; Gale's conjecture
     says no set beats the ball, and beta_{d+1}(B^d) <= sigma_d) and with the previous record for that d
     (github.com/vsvor/Gale-conjecture, best_known/, as listed in arXiv 2610.28828, Table 2).
Why this proves the bound: any partition of V into k parts is a k-colouring of G(>= dc^2); none is proper,
so some part contains two points at distance >= dc, i.e. has diameter >= dc.  Hence beta_k(V) >= dc, and the
colouring of G(> dc^2) gives equality.
Optional: --cnf FILE writes the DIMACS CNF "G(>= dc^2) is k-colourable" (variable k*i+c+1 = vertex i gets
colour c) for an independent SAT/DRAT check (e.g. cadical + drat-trim); with --sym the CNF also contains unit
clauses giving the vertices of the clique used in step 3 the colours 0,1,...,q-1 (w.l.o.g.), which makes the
DRAT proof much shorter.
"""
import os, sys, time
from fractions import Fraction
from math import gcd
from decimal import Decimal, getcontext

# (beta_k/D)^2 of the previous records, exact (recomputed from points.csv of the repository)
PREV = {
    4: (Fraction(8338402506083942692737062451903818, 10000000000000000120293888772747518),
        '61-point set best_known/4d_61p_091314 (0.913148536990)'),
    5: (Fraction(8570383072970710989222622174759524, 9999999999999999525977876802951354),
        '87-point set best_known/5d_87p_092576 (0.925763634680)'),
}

SYM = False   # --sym: the CNF also fixes the colours of one clique (unit clauses)

def read_points(path):
    pts = []
    for line in open(path):
        line = line.split('#')[0].strip()
        if not line: continue
        pts.append(tuple(Fraction(s) for s in line.replace(',', ' ').split()))
    d = len(pts[0])
    assert all(len(p) == d for p in pts), 'ragged coordinates'
    assert len(set(pts)) == len(pts), 'repeated point'
    return pts, d

def integer_coords(pts):
    den = 1
    for p in pts:
        for c in p:
            den = den * c.denominator // gcd(den, c.denominator)
    return [tuple(int(c * den) for c in p) for p in pts], den

def above_ball(r, d):
    """exact test r > sigma_d^2"""
    if d % 2 == 0:
        return r > Fraction(d + 1, d + 2)
    t = 2 * r - 1
    return t > 0 and t * t > Fraction(d - 1, d + 3)

def sigma_decimal(d):
    getcontext().prec = 50
    if d % 2 == 0:
        return (Decimal(d + 1) / Decimal(d + 2)).sqrt()
    return (Decimal(1) / 2 + (Decimal(d - 1) / Decimal(d + 3)).sqrt() / 2).sqrt()

def find_clique(n, adj, k, skip=()):
    """a large clique (size <= k) found greedily from the 40 highest-degree start vertices, avoiding the
    cliques in skip.  Any clique is fine: it is only used to fix colours (symmetry breaking)."""
    order = sorted(range(n), key=lambda v: (-len(adj[v]), v))
    best = None
    for a in order[:40]:
        cl = [a]; cand = set(adj[a])
        while cand and len(cl) < k:
            v = max(cand, key=lambda u: (len(adj[u] & cand), -u))
            cl.append(v); cand &= adj[v]
        t = tuple(sorted(cl))
        if t in skip: continue
        if best is None or len(t) > len(best): best = t
        if len(best) == k: break
    return best

def colour_fc(n, adj, clique, k):
    """Exhaustive backtracking with forward checking.  Colour symmetry is broken by fixing the clique to
    colours 0,1,...,q-1 and by never trying more than one so-far-unused colour (unused colours are
    interchangeable; the used colours always form a prefix 0..nu-1).  Returns (colouring or None, nodes)."""
    full = (1 << k) - 1
    dom = [full] * n; col = [-1] * n; nodes = [0]; nu = [0]
    nbrs = [sorted(adj[v]) for v in range(n)]
    pc = [bin(x).count('1') for x in range(1 << k)]
    trail = []
    def assign(v, c):
        stack = [(v, c)]
        while stack:
            u, cu = stack.pop()
            if col[u] != -1:
                if col[u] != cu: return False
                continue
            if not (dom[u] >> cu) & 1: return False
            assert cu <= nu[0]
            trail.append((u, dom[u], 1)); col[u] = cu; dom[u] = 1 << cu
            if cu == nu[0]:
                trail.append((-1, nu[0], 2)); nu[0] += 1
            bit = 1 << cu
            for w in nbrs[u]:
                if dom[w] & bit:
                    if col[w] != -1: return False
                    trail.append((w, dom[w], 0)); dom[w] &= ~bit
                    if dom[w] == 0: return False
                    if dom[w] & (dom[w] - 1) == 0:
                        stack.append((w, dom[w].bit_length() - 1))
        return True
    def undo(mark):
        while len(trail) > mark:
            u, dd, kind = trail.pop()
            if kind == 2: nu[0] = dd; continue
            dom[u] = dd
            if kind == 1: col[u] = -1
    def rec():
        nodes[0] += 1
        best, bs, bd = -1, k + 1, -1
        for v in range(n):
            if col[v] == -1:
                s = pc[dom[v]]
                if s < bs or (s == bs and len(nbrs[v]) > bd):
                    best, bs, bd = v, s, len(nbrs[v])
        if best == -1: return True
        for c in range(min(k, nu[0] + 1)):
            if (dom[best] >> c) & 1:
                mark = len(trail)
                if assign(best, c) and rec(): return True
                undo(mark)
        return False
    if clique is not None:
        for i, v in enumerate(clique):
            if not assign(v, i): return None, 1
    ok = rec()
    return (list(col) if ok else None), nodes[0]

def colour_plain(n, adj, clique, k, limit=10**7):
    """Plain DSATUR backtracking (no propagation): most-saturated uncoloured vertex first (ties: higher
    degree), try every colour not used by a coloured neighbour; the clique is fixed to colours 0..q-1 and at
    most one unused colour is tried.  Saturation is kept incrementally.  Returns (colouring or None, nodes)."""
    col = [-1] * n; nodes = [0]
    nbrs = [sorted(adj[v]) for v in range(n)]
    deg = [len(x) for x in nbrs]
    cnt = [[0] * k for _ in range(n)]   # cnt[v][c] = number of neighbours of v with colour c
    sat = [0] * n
    def put(v, c):
        col[v] = c
        for w in nbrs[v]:
            cnt[w][c] += 1
            if cnt[w][c] == 1: sat[w] += 1
    def take(v):
        c = col[v]; col[v] = -1
        for w in nbrs[v]:
            cnt[w][c] -= 1
            if cnt[w][c] == 0: sat[w] -= 1
    for i, v in enumerate(clique or ()):
        assert col[v] == -1 and cnt[v][i] == 0
        put(v, i)
    nu0 = len(clique or ())
    def rec(nu):
        nodes[0] += 1
        if nodes[0] > limit: raise RuntimeError('node limit')
        best, bsat, bdeg = -1, -1, -1
        for v in range(n):
            if col[v] == -1 and (sat[v] > bsat or (sat[v] == bsat and deg[v] > bdeg)):
                best, bsat, bdeg = v, sat[v], deg[v]
        if best == -1: return True
        cb = cnt[best]
        for c in range(min(k, nu + 1)):
            if cb[c] == 0:
                put(best, c)
                if rec(max(nu, c + 1)): return True
                take(best)
        return False
    ok = rec(nu0)
    return (list(col) if ok else None), nodes[0]

def graph(n, d2, s, strict):
    adj = [set() for _ in range(n)]; m = 0
    for (i, j), dd in d2.items():
        if dd > s or (not strict and dd == s):
            adj[i].add(j); adj[j].add(i); m += 1
    return adj, m

def check(path, cnf=None):
    t0 = time.time()
    pts, d = read_points(path); n = len(pts); k = d + 1
    P, den = integer_coords(pts)
    d2 = {}
    for i in range(n):
        Pi = P[i]
        for j in range(i + 1, n):
            Pj = P[j]
            d2[(i, j)] = sum((a - b) * (a - b) for a, b in zip(Pi, Pj))
    D2 = max(d2.values())
    dist = sorted(set(d2.values()), reverse=True)
    def noncol(idx):
        adj, _ = graph(n, d2, dist[idx], strict=False)
        return colour_fc(n, adj, find_clique(n, adj, k), k)[0] is None
    hi = next(i for i, x in enumerate(dist) if 10000 * x < 7225 * D2)      # first value below (0.85 D)^2
    assert noncol(hi)
    lo = 0
    while lo < hi:
        mid = (lo + hi) // 2
        if noncol(mid): hi = mid
        else: lo = mid + 1
    dc2 = dist[lo]
    adj, m = graph(n, d2, dc2, strict=False)
    cl = find_clique(n, adj, k)
    c1, nodes1 = colour_fc(n, adj, cl, k)
    assert c1 is None
    cl2 = find_clique(n, adj, k, skip={cl})
    c1b, nodes1b = colour_plain(n, adj, cl2, k)
    assert c1b is None
    adj2, m2 = graph(n, d2, dc2, strict=True)
    c2, _ = colour_fc(n, adj2, find_clique(n, adj2, k), k)
    assert c2 is not None and all(0 <= c2[v] < k for v in range(n))
    assert all(c2[v] != c2[w] for v in range(n) for w in adj2[v])
    r = Fraction(dc2, D2)
    getcontext().prec = 50
    val = (Decimal(r.numerator) / Decimal(r.denominator)).sqrt()
    ball = above_ball(r, d)
    print('file: %s   d = %d, k = %d parts, points: %d, common denominator: %d' % (path, d, k, n, den))
    print('D^2  = %d / %d^2' % (D2, den))
    print('dc^2 = %d / %d^2   (critical squared distance)' % (dc2, den))
    print('G(>= dc): %d edges; NOT %d-colourable: forward checking %d nodes (K%d %s fixed), '
          'plain DSATUR %d nodes (K%d %s fixed)' % (m, k, nodes1, len(cl), cl, nodes1b, len(cl2), cl2))
    print('G(>  dc): %d edges; explicit proper %d-colouring found and checked => beta_%d(V) = dc exactly'
          % (m2, k, k))
    print('beta_%d(V)/D(V) = sqrt(%d/%d)' % (k, r.numerator, r.denominator))
    print('               = %s...' % str(val)[:16])
    print('> ball value sigma_%d = %s...: %s' % (d, str(sigma_decimal(d))[:12], ball))
    rec = False
    if d in PREV:
        pr, label = PREV[d]
        rec = r > pr
        pv = (Decimal(pr.numerator) / Decimal(pr.denominator)).sqrt()
        print('> previous record %s: %s (gain %.4e)' % (label, rec, float(val - pv)))
    if cnf:
        units = [(v, c) for c, v in enumerate(cl)] if SYM else []
        with open(cnf, 'w') as f:
            if SYM:
                f.write('c clique %s fixed to colours 0..%d (w.l.o.g. by colour symmetry)\n' % (list(cl), len(cl) - 1))
            f.write('p cnf %d %d\n' % (k * n, n + k * m + len(units)))
            for v, c in units:
                f.write('%d 0\n' % (k * v + c + 1))
            for i in range(n):
                f.write(' '.join(str(k * i + c + 1) for c in range(k)) + ' 0\n')
            for (i, j), dd in sorted(d2.items()):
                if dd >= dc2:
                    for c in range(k):
                        f.write('-%d -%d 0\n' % (k * i + c + 1, k * j + c + 1))
        print('CNF written to %s' % cnf)
    print('time %.1f s' % (time.time() - t0))
    if rec and ball:
        print('RESULT: certified beta_%d(V)/D(V) = %s... > previous record and > sigma_%d' % (k, str(val)[:14], d))
    elif ball:
        print('RESULT: certified counterexample to Gale (above sigma_%d), not above the previous record' % d)
    else:
        print('RESULT: not above the ball value')
    return rec, ball

def main():
    args = sys.argv[1:]
    global SYM
    cnf = None
    if '--sym' in args:
        SYM = True; args.remove('--sym')
    if '--cnf' in args:
        i = args.index('--cnf'); cnf = args[i + 1]; del args[i:i + 2]
    sys.setrecursionlimit(20000)
    here = os.path.dirname(os.path.abspath(__file__))
    paths = args or [os.path.join(here, f) for f in ('points4.txt', 'points5.txt') if os.path.exists(os.path.join(here, f))]
    ok = True
    for p in paths:
        rec, ball = check(p, cnf if len(paths) == 1 else None)
        ok = ok and rec and ball
        print()
    return 0 if ok else 1

if __name__ == '__main__':
    sys.exit(main())

