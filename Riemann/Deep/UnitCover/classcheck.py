#!/usr/bin/env python3
"""classcheck.py -- exhaustive check of the moment lemma's tuple bound (DERIVATION sec 3.3) on small groups.
For G = (Z/qZ)^x (dual group = product of cyclic groups) and every k-tuple a of nonprincipal characters:
  (1) R(a) = {e in {-1,0,1}^k : sum e_i a_i = 0} is one of the closed sets enumerated by cert.closed_sets(k);
  (2) #{a : R(a) = R} <= |ker| = H^(k-rho) prod_i |G[d_i]|  (exact kernel size, Smith form);
  (3) the bound |G[d]| <= 2 d^omega(q) and |G[2]| <= tau(q) used in sec 3.4;
  (4) for k = 2, 3 the hand classification of DERIVATION sec 3.2 (stochastic domination of the pi values).
Run: uv run --no-project --with python-flint --with mpmath python3 classcheck.py   (pure integers; about 2 min)"""
import itertools, math, sys, os, time
from fractions import Fraction as Fr
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cert

def factor(n):
    f = {}; d = 2
    while d * d <= n:
        while n % d == 0: f[d] = f.get(d, 0) + 1; n //= d
        d += 1
    if n > 1: f[n] = f.get(n, 0) + 1
    return f

def cyc(q):
    ns = []
    for p, a in factor(q).items():
        if p == 2:
            if a == 2: ns.append(2)
            elif a >= 3: ns += [2, 2 ** (a - 2)]
        else: ns.append((p - 1) * p ** (a - 1))
    return [n for n in ns if n > 1]

def tau(q):
    t = 1
    for a in factor(q).values(): t *= a + 1
    return t

def Gd(ns, d):
    t = 1
    for n in ns: t *= math.gcd(n, d)
    return t

HAND = {2: lambda H, T2, T3: [(Fr(1, 4), (H - 1) ** 2), (Fr(3, 8), 2 * (H - T2)), (Fr(1, 2), T2 - 1)],
        3: lambda H, T2, T3: [(Fr(1, 8), (H - 1) ** 3), (Fr(5, 32), 4 * (H - 1) ** 2), (Fr(3, 16), 6 * (H - 1) ** 2),
                              (Fr(1, 4), 12 * (T2 - 1) * H + 12 * (H - 1) + 3 * (T2 - 1) * (H - 1)),
                              (Fr(5, 16), 4 * (H - T2)), (Fr(11, 32), 4 * (T3 - 1)), (Fr(1, 2), T2 - 1)]}

def check(q, k, closed):
    ns = cyc(q)
    G = list(itertools.product(*[range(n) for n in ns])); H = len(G)
    z = tuple(0 for _ in ns)
    nz = [x for x in G if x != z]
    cube = list(itertools.product((-1, 0, 1), repeat=k))
    index = {R: (pi, rho, ed) for (pi, rho, ed, R) in closed}
    cnt = {}; vals = []
    om = len(factor(q))
    T2 = Gd(ns, 2); T3 = Gd(ns, 3)
    assert T2 <= tau(q), (q, T2)
    for d in range(2, 13): assert Gd(ns, d) <= 2 * d ** om, (q, d)
    for a in itertools.product(nz, repeat=k):
        R = frozenset(e for e in cube if all(sum(ei * ai[j] for ei, ai in zip(e, a)) % ns[j] == 0 for j in range(len(ns))))
        if R not in index: return "R(a) not enumerated: q=%d a=%s" % (q, a)
        pi = index[R][0]
        pi_direct = sum(Fr(1, 2 ** (k + sum(1 for x in e if x))) for e in R)
        if pi != pi_direct: return "pi mismatch"
        cnt[R] = cnt.get(R, 0) + 1; vals.append(pi)
    for R, c in cnt.items():
        pi, rho, ed = index[R]
        ker = H ** (k - rho)
        for d in ed: ker *= Gd(ns, d)
        if c > ker: return "count %d > ker %d (q=%d, R rank %d, ed %s)" % (c, ker, q, rho, ed)
    if k in HAND:
        slots = HAND[k](H, T2, T3)
        for p in sorted(set(vals)):
            have = sum(1 for v in vals if v >= p); cap = sum(c for (v, c) in slots if v >= p)
            if have > cap: return "hand classes fail q=%d k=%d p=%s" % (q, k, p)
    return None

if __name__ == "__main__":
    t0 = time.time()
    plan = {2: [q for q in range(3, 121)], 3: [q for q in range(3, 61)] + [63, 65, 72, 80, 91, 104, 105, 120],
            4: [q for q in range(3, 41) if len(list(itertools.product(*[range(n) for n in cyc(q)]))) <= 12] + [16, 24, 40, 48, 60, 120]}
    bad = 0
    for k in (2, 3, 4):
        closed = cert.closed_sets(k)
        ok = 0
        for q in plan[k]:
            H = 1
            for n in cyc(q): H *= n
            if k == 4 and H > 32: continue
            msg = check(q, k, closed)
            if msg: print("FAIL", msg); bad += 1
            else: ok += 1
        print("classcheck k=%d: %d closed sets; %d moduli checked exhaustively (all tuples), %s" % (k, len(closed), ok, "ok" if bad == 0 else "FAILURES"), flush=True)
    print("classcheck: %s in %.0f s" % ("ALL OK" if bad == 0 else "%d FAILURES" % bad, time.time() - t0))
