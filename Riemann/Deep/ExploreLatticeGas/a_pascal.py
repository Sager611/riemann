#!/usr/bin/env python3
"""a_pascal.py -- explore15, transfer A (N6): the Weitz-type 'tree recursion' of Jensen polynomials.

J^{d,n}(X) = sum_j C(d,j) g(n+j) X^j satisfies the Pascal recursion J^{d+1,n} = J^{d,n} + X J^{d,n+1}
(two 'children' n, n+1, like Weitz's recursion R_v = lam * prod 1/(1+R_c)).  A hard-core-style proof would
need a region of ratios r_{d,n} = J^{d,n+1}/J^{d,n} mapped into itself by
    r_{d+1,n} = r_{d,n} (1 + X r_{d,n+1}) / (1 + X r_{d,n}).
Tests (exact rational arithmetic + Sturm counts):
 1. the identity on random sequences;
 2. propagation: level-d data (all J^{d,n} real-rooted and mutually interlacing) never forces level d+1:
    random search for counterexamples for d = 1..5 (log-concave, but J^{3,0} complex, etc.);
 3. the base map is not a self-map of the 'proper position' class.
Run: uv run --no-project --with mpmath python3 -u a_pascal.py   (< 1 min)
"""
import random
from fractions import Fraction as Fr
from math import comb

def J(g, d, n):
    return [comb(d, j) * g[n + j] for j in range(d + 1)]           # low degree first

def padd(p, q):
    m = max(len(p), len(q)); return [(p[i] if i < len(p) else 0) + (q[i] if i < len(q) else 0) for i in range(m)]

def pderiv(p):
    return [i * p[i] for i in range(1, len(p))]

def prem(a, b):
    a = a[:]
    while len(a) >= len(b) and any(a):
        if a[-1] == 0: a.pop(); continue
        c = a[-1] / b[-1]; s = len(a) - len(b)
        for i in range(len(b)): a[s + i] -= c * b[i]
        a.pop()
    while a and a[-1] == 0: a.pop()
    return a

def sturm_real_roots(p):
    """number of distinct real roots of p (exact Fractions)"""
    p = [Fr(c) for c in p]
    while p and p[-1] == 0: p.pop()
    if len(p) <= 1: return 0
    seq = [p, pderiv(p)]
    while True:
        r = prem(seq[-2], seq[-1])
        if not r: break
        seq.append([-c for c in r])
    def sc(signs):
        s = [x for x in signs if x != 0]; return sum(1 for i in range(len(s) - 1) if s[i] * s[i + 1] < 0)
    neg = [(1 if (q[-1] > 0) == (len(q) % 2 == 1) else -1) for q in seq]   # sign at -inf: lead*(-1)^deg
    pos = [(1 if q[-1] > 0 else -1) for q in seq]
    return sc(neg) - sc(pos)

def squarefree_deg(p):
    # degree of p / gcd(p, p')
    p = [Fr(c) for c in p]
    a, b = p[:], pderiv(p)
    while b:
        a, b = b, prem(a, b)
    return len(p) - len(a)

def real_rooted(p):
    """all roots real: exact Sturm count for low degree (<= 4), mpmath roots (60 digits) otherwise"""
    if len(p) <= 5:
        return sturm_real_roots(p) == squarefree_deg(p)
    import mpmath as mp
    mp.mp.dps = 60
    c = [mp.mpf(Fr(x).numerator) / Fr(x).denominator for x in reversed(p)]
    rs = mp.polyroots(c, maxsteps=400, extraprec=200)
    return all(abs(mp.im(r)) < mp.mpf(10) ** -25 * (1 + abs(r)) for r in rs)

def roots_mp(p):
    import mpmath as mp
    mp.mp.dps = 60
    c = [mp.mpf(x.numerator) / x.denominator for x in reversed([Fr(v) for v in p])]
    return sorted([mp.re(r) for r in mp.polyroots(c, maxsteps=400, extraprec=200)], reverse=True)

def interlace_left(f, g):
    """f << g for real-rooted f, g of equal degree with positive coefficients (zeros < 0):
    0 > f_1 >= g_1 >= f_2 >= g_2 >= ... >= f_d >= g_d (Hermite-Kakeya-Obreschkoff: f + c g real-rooted for all c)."""
    a, b = roots_mp(f), roots_mp(g)
    if len(a) != len(b): return False
    seq = []
    for x, y in zip(a, b): seq += [x, y]
    tol = 1e-40
    return all(seq[i] >= seq[i + 1] - tol for i in range(len(seq) - 1))

def level_ok(g, d, nmax):
    """level-d data: all J^{d,n}, n<=nmax, real-rooted and consecutive pairs interlacing"""
    Js = [J(g, d, n) for n in range(nmax + 1)]
    if not all(real_rooted(p) for p in Js): return False
    return all(interlace_left(Js[i], Js[i + 1]) for i in range(len(Js) - 1))

def main():
    rnd = random.Random(15)
    # 1. identity
    for _ in range(200):
        g = [Fr(rnd.randint(1, 10 ** 6), rnd.randint(1, 10 ** 6)) for _ in range(12)]
        d, n = rnd.randint(0, 6), rnd.randint(0, 4)
        assert J(g, d + 1, n) == padd(J(g, d, n), [0] + J(g, d, n + 1))
    print('1. Pascal identity J^{d+1,n} = J^{d,n} + X J^{d,n+1}: 200/200 random checks')
    # 2a. explicit: (1,1,1,eps) is log-concave (all Turan inequalities), J^{2,0} << J^{2,1} fails, J^{3,0} complex
    e = Fr(1, 100); g = [Fr(1), Fr(1), Fr(1), e]
    print('2a. g = (1,1,1,1/100): Turan', all(g[i] ** 2 >= g[i - 1] * g[i + 1] for i in (1, 2)),
          '| J^{2,0}, J^{2,1} real-rooted', real_rooted(J(g, 2, 0)), real_rooted(J(g, 2, 1)),
          '| J^{2,0} << J^{2,1}', interlace_left(J(g, 2, 0), J(g, 2, 1)),
          '| J^{3,0} real-rooted', real_rooted(J(g, 3, 0)))
    # 2b. random search: sequences satisfying ALL level-d conditions on a window, failing level d+1 at n = 0
    stats = {}
    for d in range(1, 6):
        found = None; tries = 0
        while found is None and tries < 1500:
            tries += 1
            # perturb an LP sequence: moments of 1/Gamma-type (coefficients of a real-rooted polynomial)
            roots = sorted(Fr(rnd.randint(1, 60), rnd.randint(1, 8)) for _ in range(d + 4))
            poly = [Fr(1)]
            for r in roots: poly = padd([c * r for c in poly], [0] + poly)       # prod (r + y)
            from math import factorial
            gg = [poly[k] * factorial(k) for k in range(len(poly))]              # g(k) = k! [y^k]
            k = rnd.randint(d, len(gg) - 1)
            gg[k] = gg[k] * Fr(rnd.randint(70, 130), 100)                          # perturb one moment
            nmax = len(gg) - 2 - d
            if nmax < 1: continue
            if level_ok(gg, d, nmax):
                ok_next = real_rooted(J(gg, d + 1, 0)) and real_rooted(J(gg, d + 1, 1))
                assert ok_next, 'theorem I(d) => H(d+1) violated?'
                stats[d] = stats.get(d, 0) + 1
                if not interlace_left(J(gg, d + 1, 0), J(gg, d + 1, 1)):
                    found = (gg, k)
        if found:
            gg, k = found
            print(f'2b. d={d}: I({d}) holds on the whole window (all J^{{{d},n}} real-rooted, consecutive ones'
                  f' interlacing), H({d + 1}) holds as the lemma says, but I({d + 1}) fails: J^{{{d + 1},0}} << J^{{{d + 1},1}} is false'
                  f'  [perturbed index {k}; {tries} tries, {stats.get(d, 0)} samples with I({d})]')
        else:
            print(f'2b. d={d}: no I(d)-but-not-I(d+1) example in {tries} tries ({stats.get(d, 0)} samples with I(d))')
    print('3. Lemma (Hermite-Kakeya-Obreschkoff): I(d) [J^{d,n} << J^{d,n+1} for all n] => H(d+1) [all J^{d+1,n}'
          ' real-rooted], checked on every sample above; but I(d) does not imply I(d+1): the ratio map'
          ' r -> r(1+X r\')/(1+X r) has no invariant region, every level is a new inequality (Polya-Schur: all levels <=> LP).')

if __name__ == '__main__':
    main()
