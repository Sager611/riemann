#!/usr/bin/env python3
"""Sanity checks (floating point, not part of the proof).
  uv run --no-project --with numpy --with mpmath python3 sanity.py > logs/sanity_output.txt
1. least primes in every class of Gal(Q(zeta_q)/Q), 3 <= q <= 2000, against both bounds;
2. least primes in the three classes of the S_3 closures of the 508 cubic fields with |d| <= 3000 (logs/s3_least.txt, PARI);
3. explicit-formula normalisation: |H| S(x) against x Phi(1) for Q(zeta_q)/Q, q = 5, 7, 11, 13, x = 10^7,
   with the error bound of DERIVATION (3.5) at delta = 1/32 and in the GRH shape delta = 1/2;
4. the Hadamard identity (DERIVATION 3.3) for zeta (L = Q) at sigma = 2 with 2000 zeros.
"""
import math, numpy as np, mpmath as mp

D_CT = 1 / 32
B2 = (777 - math.sqrt(1081)) / 852
D_L = 1 - B2
C_CT, C_L = 0.484, 0.62

def sieve(n):
    s = np.ones(n + 1, bool); s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]: s[i * i::i] = False
    return np.nonzero(s)[0]
P = sieve(10_000_000)

def phi(n):
    r = n; m = n; p = 2
    while p * p <= m:
        if m % p == 0:
            while m % p == 0: m //= p
            r -= r // p
        p += 1
    if m > 1: r -= r // m
    return r

def logdisc_cyclo(q):
    f = phi(q); s = f * math.log(q); m = q; p = 2
    while p * p <= m:
        if m % p == 0:
            while m % p == 0: m //= p
            s -= f / (p - 1) * math.log(p)
        p += 1
    if m > 1: s -= f / (m - 1) * math.log(m)
    return s

print('1. Q(zeta_q)/Q, 3 <= q <= 2000, every class a mod q')
worst = (0, None); worstL = (0, None); n_classes = 0
for q in range(3, 2001):
    if q % 4 == 2: continue           # Q(zeta_q) = Q(zeta_{q/2})
    L = logdisc_cyclo(q); ell = L + math.log(3)
    res = P % q; first = {}
    for p, r in zip(P[:400000], res[:400000]):
        if math.gcd(int(r), q) == 1 and r not in first: first[r] = int(p)
        if len(first) == phi(q): break
    assert len(first) == phi(q), q
    pm = max(first.values()); n_classes += phi(q)
    e1 = math.exp(D_CT * math.log(pm)) / ell; e2 = math.exp(D_L * math.log(pm)) / ell
    if e1 > worst[0]: worst = (e1, (q, pm))
    if e2 > worstL[0]: worstL = (e2, (q, pm))
    assert pm <= (C_CT * ell) ** 32 and pm <= (C_L * ell) ** (1 / D_L)
print(f'   {n_classes} classes; all least primes below both bounds.')
print(f'   max p^(1/32)/log(3 Delta) = {worst[0]:.4f} at (q, p) = {worst[1]}  (theorem: 0.484)')
print(f'   max p^(1-B2)/log(3 Delta) = {worstL[0]:.4f} at (q, p) = {worstL[1]}  (theorem: 0.62)')

print('\n2. S_3 closures of cubic fields, |d| <= 3000 (PARI nflist; classes id, transposition, 3-cycle)')
rows = [list(map(int, l.split())) for l in open('logs/s3_least.txt') if l.strip()]
w1 = w2 = (0, None)
for d, DL, a, b, c in rows:
    ell = math.log(3 * DL); pm = max(a, b, c)
    assert pm <= (C_CT * ell) ** 32 and pm <= (C_L * ell) ** (1 / D_L)
    e1 = pm ** D_CT / ell; e2 = pm ** D_L / ell
    if e1 > w1[0]: w1 = (e1, (d, pm))
    if e2 > w2[0]: w2 = (e2, (d, pm))
print(f'   {len(rows)} fields, {3 * len(rows)} classes; all least primes below both bounds.')
print(f'   max p^(1/32)/log(3 Delta_L) = {w1[0]:.4f} at (d, p) = {w1[1]};  max p^(1-B2)/log(3 Delta_L) = {w2[0]:.4f} at {w2[1]}')

print('\n3. explicit formula normalisation, L = Q(zeta_q), K = Q, C = {g}, g a generator, E = Q, |H| = q-1')
h, ss, s0 = 2.0, 0.0, 2.0
U = 3 * h
def Fk(u):
    v = u / h
    p2 = lambda t: np.where(t > 0, t, 0.0) ** 2
    return np.exp(ss * u) * h * h * 0.5 * (p2(v) - 3 * p2(v - 1) + 3 * p2(v - 2) - p2(v - 3))
def Phi(s):
    z = s - ss
    return ((1 - np.exp(-h * z)) / z) ** 3 if abs(z) > 1e-12 else h ** 3
P1, P0 = Phi(1.0).real, Phi(0.0).real
gam = np.linspace(0, 200, 400001)
def Mline(beta):
    vals = np.abs(((1 - np.exp(-h * (beta + 1j * gam - ss))) / (beta + 1j * gam - ss)) ** 3) * ((s0 - beta) ** 2 + gam ** 2) / (s0 - beta)
    return vals.max()
g_s0 = min((math.log(math.pi) - float(mp.digamma(s0 / 2))) / 2, (math.log(2 * math.pi) - float(mp.digamma(s0))) / 2)
x = 1e7
for q in [5, 7, 11, 13]:
    g = int(mp.primitive_root(q)) if hasattr(mp, 'primitive_root') else next(a for a in range(2, q) if len({pow(a, k, q) for k in range(1, q)}) == q - 1)
    n = q - 1; Lam = (q - 2) * math.log(q)
    S = 0.0; good = 0.0
    for p in P[P <= x]:
        p = int(p); lp = math.log(p); pk = p; m = 1
        while pk <= x:
            u = math.log(x / pk)
            if u <= U:
                kv = float(Fk(np.array([u]))[0])
                th = 1.0 / (q - 1) if p == q else (1.0 if pk % q == g else 0.0)
                S += lp * th * kv
                if m == 1 and p != q and th == 1.0: good += lp * kv
            pk *= p; m += 1
    eps = math.exp(U) / (x - math.exp(U))
    B = Lam / 2 + 1 / s0 + 1 / (s0 - 1) - n * g_s0
    triv = (n - 1) * P0 + n * P0 * eps
    err = abs(n * S - x * P1)
    bnd_ct = x ** (1 - D_CT) * max(Mline(1 - D_CT), x ** (2 * D_CT - 1) * Mline(D_CT)) * B + triv
    bnd_grh = x ** 0.5 * Mline(0.5) * B + triv
    print(f'   q={q:2d} g={g}: |H|S(x)/(x Phi(1)) = {n * S / (x * P1):.6f};  |error| = {err:.4g}  <= bound(1/32) = {bnd_ct:.4g}: {err <= bnd_ct};'
          f'  <= GRH-shape bound(1/2) = {bnd_grh:.4g}: {err <= bnd_grh};  |H| Good/(x Phi(1)) = {n * good / (x * P1):.6f}')

print('\n4. Hadamard identity for zeta at sigma = 2: sum_rho Re 1/(2-rho) = 1/2 + 1 + (psi(1) - log pi)/2 + zeta\'/zeta(2)')
mp.mp.dps = 20
rhs = 0.5 + 1 + (mp.digamma(1) - mp.log(mp.pi)) / 2 + mp.zeta(2, derivative=1) / mp.zeta(2)
N = 2000; s = mp.mpf(0); T = 0
for k in range(1, N + 1):
    t = mp.im(mp.zetazero(k)); T = t
    s += 2 * mp.mpf(1.5) / (mp.mpf(2.25) + t * t)
tail = 3 / (2 * mp.pi) * (mp.log(T / (2 * mp.pi)) + 1) / T
print(f'   RHS = {mp.nstr(rhs, 10)};  sum over first {N} zeros = {mp.nstr(s, 10)};  + tail estimate {mp.nstr(tail, 3)} = {mp.nstr(s + tail, 8)}')
