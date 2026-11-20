#!/usr/bin/env python3
# alpaut probe (vectorized): test the archive's LemmaX (AV/LevelCohomology/UnramifiedCocycleIsCoboundaryOfBernoulliUnit.lean:14-27)
# XData l d phi xi: phi primitive mod odd f, ord phi | l^d-1, 4 | l^d-1, xi primitive mod 2^n, f*4^{v2(l^d-1)} <= 2^n,
# phi(-1)xi(-1) = -1  ==>  ||B_{1,phi xi}||_l = 1 (every embedding).  Since l does not divide N = f*2^n, a prime lambda | l divides
# B_{1,chi} iff it divides P(zeta_M) = N*B_{1,chi} = sum_{a<=N,(a,N)=1} a*zeta_M^{k(a)}; all lambda at once: gcd(P, Phi_M) over GF(l).
# Galois conjugates (phi^k, xi^k) also satisfy XData, so gcd != 1 means LemmaX fails for some embedding.
# usage: lemmax_probe.py MODE LMAX FMAX EXTRA   (MODE thr: n = nmin..nmin+EXTRA-1 ; below: n = nmin-EXTRA..nmin-1)
import sys, math
import numpy as np
from flint import nmod_poly, fmpz_poly
from sympy import primerange, primitive_root, n_order

def v2(x):
    c = 0
    while x % 2 == 0: x //= 2; c += 1
    return c

def test(l, f, m, j, n, s):
    N = f * 2**n; q = 2**(n-2); M = math.lcm(m, q, 2)
    lg = np.zeros(2**n, dtype=np.int64); sg = np.zeros(2**n, dtype=np.int64); x = 1
    for t in range(q):
        lg[x] = t; lg[(-x) % 2**n] = t; sg[(-x) % 2**n] = 1; x = x * 5 % 2**n
    a = np.arange(1, N + 1, dtype=np.int64)
    keep = (a % 2 == 1)
    k = (M // q) * lg[a % 2**n]
    if s == -1: k = k + (M // 2) * sg[a % 2**n]
    if f > 1:
        g = primitive_root(f); lf = np.zeros(f, dtype=np.int64); y = 1
        for t in range(f - 1): lf[y] = t; y = y * g % f
        keep &= (a % f != 0)
        k = k + (M // m) * j * lf[a % f]
    k = k[keep] % M; w = a[keep] % l
    coeffs = np.bincount(k, weights=w.astype(np.float64), minlength=M)
    P = nmod_poly([int(c) % l for c in coeffs], l)
    C = nmod_poly(fmpz_poly.cyclotomic(M).coeffs(), l)
    return P.gcd(C).degree(), M, n_order(l, M)

def run(mode, Lmax, Fmax, extra):
    hits = 0; tests = 0
    for l in primerange(3, Lmax):
        print(mode, 'l', l, 'tests so far', tests, flush=True)
        for f in [1] + list(primerange(3, Fmax)):
            if f == l: continue
            ms = [1] if f == 1 else [m for m in range(2, f) if (f - 1) % m == 0 and math.gcd(m, l) == 1]
            for m in ms:
                d = 1
                while (l**d - 1) % m or (l**d - 1) % 4: d += 1
                v = v2(l**d - 1); nmin = max(3, math.ceil(math.log2(f * 4**v) - 1e-12))
                ns = range(nmin, nmin + extra) if mode == 'thr' else range(max(4, nmin - extra), nmin)
                for n in ns:
                    if 2**(n-2) * m > 1 << 13 or f * 2**n > 1 << 19: continue
                    js = [j for j in range(1, m + 1) if math.gcd(j, m) == 1] if m > 1 else [0]
                    for j in js[:2]:
                        pphi = 1 if f == 1 else (1 if (j * ((f - 1) // 2)) % m == 0 else -1)
                        s = -pphi
                        deg, M, F = test(l, f, m, j, n, s); tests += 1
                        if deg: hits += 1
                        if deg:
                            print(mode, 'l=%d d=%d v=%d f=%d m=%d j=%d n=%d nmin=%d xi(-1)=%d M=%d F=%d gcd_deg=%d DIVISIBLE' % (l, d, v, f, m, j, n, nmin, s, M, F, deg), flush=True)
    print(mode, 'Lmax', Lmax, 'Fmax', Fmax, 'extra', extra, 'tests', tests, 'divisible', hits, flush=True)

if __name__ == '__main__':
    run(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]))
