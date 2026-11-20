#!/usr/bin/env python3
# alpaut probe P1b: LemmaX with squarefree composite odd f and with l | f (both allowed by XData,
# AV/LevelCohomology/UnramifiedCocycleIsCoboundaryOfBernoulliUnit.lean:14-22). phi = prod over p | f of phi_p (order m_p | p-1,
# nontrivial => primitive mod p), ord phi = lcm m_p | l^d - 1 (so gcd(m_p, l) = 1).  B = P(zeta_M)/N, N = f*2^n.
# l not dividing f: violation iff gcd(P, Phi_M) != 1 mod l.  l | f (v_l(N) = 1): ||B|| = 1 for all embeddings iff v_lambda(P) = 1 for
# every lambda | l: check (a) Phi_M | P mod l and (b) v_l(Res(Phi_M, P)) = phi(M) (sum of residue degrees).
# usage: lemmax_probe2.py LMAX PMAX NPRIMES EXTRA
import sys, math, itertools
import numpy as np
from flint import nmod_poly, fmpz_poly, fmpz
from sympy import primerange, primitive_root, n_order, totient, divisors

def v2(x):
    c = 0
    while x % 2 == 0: x //= 2; c += 1
    return c

def vl(x, l):
    x = abs(int(x)); c = 0
    if x == 0: return 10**9
    while x % l == 0: x //= l; c += 1
    return c

def coeffs_P(ps, ms, n, s):
    f = math.prod(ps); N = f * 2**n; q = 2**(n-2); m = math.lcm(*ms) if ms else 1; M = math.lcm(m, q, 2)
    lg = np.zeros(2**n, dtype=np.int64); sg = np.zeros(2**n, dtype=np.int64); x = 1
    for t in range(q):
        lg[x] = t; lg[(-x) % 2**n] = t; sg[(-x) % 2**n] = 1; x = x * 5 % 2**n
    a = np.arange(1, N + 1, dtype=np.int64)
    keep = (a % 2 == 1)
    k = (M // q) * lg[a % 2**n]
    if s == -1: k = k + (M // 2) * sg[a % 2**n]
    for p, mp in zip(ps, ms):
        # p may be a prime power q = r^k (cyclic unit group, phi(q) elements)
        g = primitive_root(p); lf = np.zeros(p, dtype=np.int64); y = 1
        for t in range(int(totient(p))): lf[y] = t; y = y * g % p
        r0 = min(x for x in range(3, p + 1) if p % x == 0)
        keep &= (a % r0 != 0)
        k = k + (M // mp) * lf[a % p]
    k = k[keep] % M; w = a[keep]
    return np.bincount(k, weights=w.astype(np.float64), minlength=M), M, N

def parity(ps, ms):
    s = 1
    for p, mp in zip(ps, ms):
        # -1 = g^(phi(p)/2) for cyclic (Z/p)^x, phi_p(g^t) = zeta_mp^t
        s *= (1 if (int(totient(p)) // 2) % mp == 0 else -1)
    return s

def run_pp(Lmax, extra):
    # prime-power conductors f = r^2, r^3 (r odd prime != l); phi primitive mod r^k <=> r | ord phi; need ord phi | l^d - 1
    tests = viol = 0
    for l in primerange(3, Lmax):
        for q in [9, 27, 25, 49, 121, 169]:
            r0 = min(x for x in range(3, q + 1) if q % x == 0)
            if r0 == l: continue
            for mq in [m for m in divisors(int(totient(q))) if m % r0 == 0 and math.gcd(m, l) == 1][:2]:
                d = 1
                while (l**d - 1) % mq or (l**d - 1) % 4: d += 1
                v = v2(l**d - 1); nmin = max(3, math.ceil(math.log2(q * 4**v) - 1e-12))
                for n in range(nmin, nmin + extra):
                    if 2**(n-2) * mq > 1 << 12 or q * 2**n > 1 << 19: continue
                    s = -parity([q], [mq])
                    cf, M, N = coeffs_P([q], [mq], n, s); tests += 1
                    C = fmpz_poly.cyclotomic(M)
                    g = nmod_poly([int(c) % l for c in cf], l).gcd(nmod_poly(C.coeffs(), l))
                    if g.degree() > 0:
                        viol += 1; print('VIOLATION pp l=%d f=%d ord=%d n=%d M=%d gcd_deg=%d' % (l, q, mq, n, M, g.degree()), flush=True)
    print('prime-power tests', tests, 'violations', viol, flush=True)

def run(Lmax, Pmax, npr, extra):
    tests = viol = 0
    for l in primerange(3, Lmax):
        print('l', l, 'tests', tests, flush=True)
        pool = [p for p in primerange(3, Pmax)]
        for r in range(1, npr + 1):
            for ps in itertools.combinations(pool, r):
                if math.prod(ps) > 300: continue
                choices = [[mp for mp in divisors(p - 1) if mp > 1 and math.gcd(mp, l) == 1] for p in ps]
                if any(not c for c in choices): continue
                for ms in itertools.product(*[c[-1:] + c[:1] for c in choices]):
                    ms = list(ms); m = math.lcm(*ms); f = math.prod(ps)
                    d = 1
                    while (l**d - 1) % m or (l**d - 1) % 4: d += 1
                    v = v2(l**d - 1); nmin = max(3, math.ceil(math.log2(f * 4**v) - 1e-12))
                    for n in range(nmin, nmin + extra):
                        if 2**(n-2) * m > 1 << 12 or f * 2**n > 1 << 19: continue
                        s = -parity(ps, ms)
                        cf, M, N = coeffs_P(ps, ms, n, s); tests += 1
                        C = fmpz_poly.cyclotomic(M)
                        if l not in ps:
                            g = nmod_poly([int(c) % l for c in cf], l).gcd(nmod_poly(C.coeffs(), l))
                            bad = g.degree() > 0; info = 'gcd_deg=%d' % g.degree()
                        else:
                            Pm = nmod_poly([int(c) % l for c in cf], l); Cm = nmod_poly(C.coeffs(), l)
                            a_ok = (Pm % Cm).degree() < 0 if hasattr(Pm % Cm, 'degree') else False
                            P = fmpz_poly([int(c) for c in cf]); R = C.resultant(P)
                            b = vl(R, l); phiM = int(totient(M))
                            bad = (not a_ok) or b != phiM; info = 'PhiM|P mod l=%s v_l(Res)=%d phi(M)=%d' % (a_ok, b, phiM)
                        if bad:
                            viol += 1
                            print('VIOLATION l=%d d=%d v=%d f=%s orders=%s n=%d nmin=%d xi(-1)=%d M=%d %s' % (l, d, v, ps, ms, n, nmin, s, M, info), flush=True)
                        elif l in ps and tests % 10 == 0:
                            print('  ok l|f: l=%d f=%s orders=%s n=%d %s' % (l, ps, ms, n, info), flush=True)
    print('tests', tests, 'violations', viol, flush=True)

if __name__ == '__main__':
    if sys.argv[1] == 'pp': run_pp(int(sys.argv[2]), int(sys.argv[3]))
    else: run(int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]))
