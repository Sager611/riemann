#!/usr/bin/env python3
"""Small-analog tests of Lemma arith:injectivity and (arith:divisibility), arXiv 2610.30476 (TeX l.534-704).
The proofs only need: p0 odd, (d/p0) = -1, m = ceil(sqrt p0), 1 <= s0 < p0. We build the integer matrix of
iota for small p0 and several s0 (incl. the paper's s0 = floor(p0/30) where >= 1) and certify full column rank
by Gaussian elimination modulo two large primes (rank mod l = n  =>  rank over Q = n  =>  injective on Z^n).
Divisibility iota(g_{r,s}^{(p)}) in p^h Z^M is checked with exact integers in Z[sqrt d]."""
import numpy as np, math, sys, itertools
L1, L2 = 2147483647, 2147483629
def leg(d, p): 
    r = pow(d % p, (p - 1) // 2, p)
    return 0 if d % p == 0 else (1 if r == 1 else -1)
def is_prime(n):
    return n > 1 and all(n % k for k in range(2, math.isqrt(n) + 1))
def zmul(a, b, d, mod=None):
    x = a[0] * b[0] + d * a[1] * b[1]; y = a[0] * b[1] + a[1] * b[0]
    return (x % mod, y % mod) if mod else (x, y)
def powers(al, n, d, mod=None):
    out = [(1, 0)]
    for _ in range(n): out.append(zmul(out[-1], al, d, mod))
    return out
def iota_matrix(p0, s0, d, mod):
    m = math.isqrt(p0 - 1) + 1
    cols = [(r, s) for r in range(p0) for s in range(s0)]
    rows = []
    for u in range(m):
        for v in range(1, m):
            A = powers((u, v), p0, d, mod); B = powers((u, -v), s0, d, mod)
            for j in range(s0):
                r0 = []; r1 = []
                for (r, s) in cols:
                    if j > r: r0.append(0); r1.append(0); continue
                    ff = math.perm(r, j) % mod
                    val = zmul(A[r - j], B[s], d, mod)
                    r0.append(ff * val[0] % mod); r1.append(ff * val[1] % mod)
                rows.append(r0); rows.append(r1)
    return np.array(rows, dtype=np.int64), len(cols), len(rows)
def rank_mod(M, mod):
    M = M.copy() % mod; nr, nc = M.shape; rk = 0
    for c in range(nc):
        piv = None
        nz = np.nonzero(M[rk:, c])[0]
        if nz.size == 0: continue
        piv = rk + int(nz[0])
        M[[rk, piv]] = M[[piv, rk]]
        inv = pow(int(M[rk, c]), mod - 2, mod)
        M[rk] = (M[rk] * inv) % mod
        col = M[:, c].copy(); col[rk] = 0
        nzr = np.nonzero(col)[0]
        if nzr.size:
            M[nzr] = (M[nzr] - (col[nzr, None] * M[rk][None, :]) % mod) % mod
        rk += 1
        if rk == nr: break
    return rk
print("(1) injectivity: full column rank of iota (n = p0*s0 unknowns)")
cases = []
for p0 in [7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 61, 67, 97]:
    for d in [-1, -2, -3, 2, 3, 5, -7, 6, -163]:
        if leg(d, p0) != -1: continue
        for s0 in sorted(set([1, 2, 3, max(1, p0 // 30)] + ([p0 - 1] if p0 <= 13 else []))):
            if s0 >= p0 or p0 * s0 > 400: continue
            cases.append((p0, s0, d))
bad = []
for (p0, s0, d) in cases:
    M, n, R = iota_matrix(p0, s0, d, L1)
    rk = rank_mod(M, L1)
    if rk < n:
        M2, _, _ = iota_matrix(p0, s0, d, L2); rk = max(rk, rank_mod(M2, L2))
    if rk < n: bad.append((p0, s0, d, rk, n))
print(f"    {len(cases)} cases (p0<=97, s0 in {{1,2,3,floor(p0/30),p0-1}}, n<=400): rank deficient = {bad if bad else 'none'}")
# also: split primes (proof does not apply) -- informational
spl = []
for p0 in [7, 11, 13, 17, 19, 23, 29, 31]:
    for d in [-1, -2, -3, 2, 3, 5]:
        if leg(d, p0) == 1:
            for s0 in [1, 2, 3]:
                M, n, R = iota_matrix(p0, s0, d, L1); spl.append(rank_mod(M, L1) == n)
print(f"    (info) split p0 analog: full rank in {sum(spl)}/{len(spl)} cases")
# (2) divisibility, exact
print("(2) divisibility iota(g_{r,s}^{(p)}) in p^h Z^M, exact integers")
def iota_poly(poly, p0, s0, d):
    m = math.isqrt(p0 - 1) + 1; out = []
    maxa = max(a for a, b in poly) if poly else 0; maxb = max(b for a, b in poly) if poly else 0
    for u in range(m):
        for v in range(1, m):
            A = powers((u, v), maxa, d); B = powers((u, -v), maxb, d)
            for j in range(s0):
                x = y = 0
                for (a, b), c in poly.items():
                    if a < j: continue
                    ff = math.perm(a, j); val = zmul(A[a - j], B[b], d)
                    x += c * ff * val[0]; y += c * ff * val[1]
                out += [x, y]
    return out
def vp(nv, p):
    if nv == 0: return 10**9
    k = 0
    while nv % p == 0: nv //= p; k += 1
    return k
tot = 0; fails = []; frob = 0; frobn = 0
pairs = []
for d, s0 in [(-1, 4), (-2, 3), (2, 5), (-3, 4), (3, 3), (5, 3), (-7, 3), (-163, 2)]:
    p0 = next(p for p in range(19, 200) if is_prime(p) and leg(d, p) == -1)
    pairs.append((p0, s0, d))
print("    pairs (p0, s0, d):", pairs)
for p0, s0, d in pairs:
    assert leg(d, p0) == -1
    Pset = [p for p in range(3, p0 + 1) if is_prime(p) and leg(d, p) == -1]
    for p in Pset:
        for u in range(5):
            for v in range(1, 5):
                al = powers((u, v), p, d)[p]
                frob += ((al[0] - u) % p == 0 and (al[1] + v) % p == 0); frobn += 1
        for r in range(p0):
            for s in range(s0):
                h = min(r // p, s0 - 1 - s)
                if h == 0: continue
                poly = {}
                for k in range(h + 1):
                    key = (r - p * k, s + k); poly[key] = poly.get(key, 0) + (-1) ** k * math.comb(h, k)
                vec = iota_poly(poly, p0, s0, d)
                mv = min(vp(abs(x), p) for x in vec)
                tot += 1
                if mv < h: fails.append((p0, s0, d, p, r, s, h, mv))
print(f"    {tot} basis vectors with h>=1 checked: failures = {fails if fails else 'none'}; Frobenius rho((u+vt)^p) = u-vt mod p held {frob}/{frobn} times")


# (3) negative control: for SPLIT p (chi(p)=+1) the same divisibility must fail somewhere (test is sensitive)
neg = []
for p0, s0, d in pairs[:4]:
    Sp = [p for p in range(3, p0 + 1) if is_prime(p) and leg(d, p) == 1]
    for p in Sp:
        failed = False
        for r in range(p, p0):
            for s in range(s0):
                h = min(r // p, s0 - 1 - s)
                if h == 0: continue
                poly = {}
                for k in range(h + 1):
                    key = (r - p * k, s + k); poly[key] = poly.get(key, 0) + (-1) ** k * math.comb(h, k)
                if min(vp(abs(x), p) for x in iota_poly(poly, p0, s0, d)) < h: failed = True; break
            if failed: break
        neg.append((d, p, failed))
print("(3) negative control, split p: divisibility fails as it should in", sum(f for _, _, f in neg), "/", len(neg), "cases", neg)
