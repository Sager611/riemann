# alpext P3: the windowed Landau endgame (L1; ZeroFreeOfRpowWindow.lean:11-14, CharFamily.lean:24-36) for L(s,rho) = zeta_K/zeta,
# K the cubic field of discriminant -23, rho = Ind psi (psi the cubic class-group character of Q(sqrt-23)) = L(s, eta(z)eta(23z)).
# (1) -L'/L(s,rho) = sum Lambda(n) b(n) n^-s with b(p^k) = tr rho(Frob_p^k) (real), b(23^k) = 1: exact check against the
#     q-expansion of eta(z)eta(23z) (Satake recursion), p < 2000.
# (2) completed Lambda(s) = (sqrt23/2pi)^s Gamma(s) L(s): functional equation with root number +1 (t-independence of the
#     incomplete-gamma series), L(1,rho) = 2 pi log(alpha)/sqrt23 (class number formula of K), L'/L(1,rho).
# (3) the analogue of Gm: Gz + w*G_rho has coefficients Lambda(n)(1 + w tr) >= 0 for w <= 1 (exact, n <= 1e6);
#     G_rho(y) = sum Lambda(n) b(n)(y/n - 1) + y L'/L(1,rho) against the GRH envelope 2 Lambda'/Lambda(1) sqrt(y), y <= 1e6.
import mpmath as mp
from math import isqrt, log, sqrt
mp.mp.dps = 30
NQ = 2000
eta = [0] * (NQ + 1)
k = 0
while True:
    done = True
    for kk in (k, -k) if k else (0,):
        e = kk * (3 * kk - 1) // 2
        if e <= NQ: eta[e] += (-1) ** abs(kk); done = False
    if done and k > 0: break
    k += 1
f = [0] * (NQ + 1)  # q * eta(q) * eta(q^23) / q^(24/24): eta(z)eta(23z) = q prod(1-q^n)(1-q^23n)
for i, ci in enumerate(eta):
    if not ci: continue
    for j in range(0, (NQ - 1 - i) // 23 + 1):
        cj = eta[j]
        if cj and 1 + i + 23 * j <= NQ: f[1 + i + 23 * j] += ci * cj
def leg(a, p): r = pow(a % p, (p - 1) // 2, p); return -1 if r == p - 1 else r
def pmulmod(A, B, p):
    r = [0] * 5
    for i, ai in enumerate(A):
        for j, bj in enumerate(B): r[i + j] += ai * bj
    for kk in (4, 3):
        if r[kk]: r[kk - 2] += r[kk]; r[kk - 3] += r[kk]; r[kk] = 0
    return [v % p for v in r[:3]]
def xpow(p):
    R = [1, 0, 0]; B = [0, 1, 0]; e = p
    while e:
        if e & 1: R = pmulmod(R, B, p)
        B = pmulmod(B, B, p); e >>= 1
    return R
def ftype(p):
    if p == 23: return 'ram'
    if p == 2: return '3c'      # x^3-x-1 irreducible mod 2; (-23/2): 2 splits in F (-23 = 1 mod 8)
    L = leg(-23, p)
    return 'tr' if L == -1 else ('id' if xpow(p) == [0, 1, 0] else '3c')
def trk(t, k):
    if t == 'ram': return 1
    if t == 'id': return 2
    if t == 'tr': return 2 if k % 2 == 0 else 0
    return 2 if k % 3 == 0 else -1
def chi23(p): return 0 if p == 23 else (1 if p == 2 else leg(-23, p))
N = 10 ** 6
sv = bytearray([1]) * (N + 1); sv[0] = sv[1] = 0
for i in range(2, isqrt(N) + 1):
    if sv[i]: sv[i * i::i] = bytearray(len(range(i * i, N + 1, i)))
primes = [p for p in range(2, N + 1) if sv[p]]
FT = {p: ftype(p) for p in primes}
# (1) Satake check: a(p) from eta product, b(p^k) = alpha^k + beta^k via recursion s_k = a s_{k-1} - chi s_{k-2}
bad = 0; nchk = 0
for p in primes:
    if p >= NQ: break
    a = f[p]; c = chi23(p)
    if p != 23:
        s = [2, a]
        for kk in range(2, 7): s.append(a * s[-1] - c * s[-2])
        for kk in range(1, 7):
            nchk += 1
            if s[kk] != trk(FT[p], kk): bad += 1
    else:
        nchk += 1
        if a != 1: bad += 1
# multiplicativity of a(n) from the eta product, n < 2000 (Hecke)
def amult(n):
    r = 1; m = n
    for p in primes:
        if p * p > m: break
        if m % p == 0:
            e = 0
            while m % p == 0: m //= p; e += 1
            a = f[p]; c = chi23(p); s = [1, a]
            for kk in range(2, e + 1): s.append(a * s[-1] - c * s[-2])
            r *= s[e]
    if m > 1: r *= f[m]
    return r
badm = sum(1 for n in range(1, NQ + 1) if f[n] != amult(n))
print(f'(1) b(p^k) = tr rho(Frob_p^k) from eta(z)eta(23z): {nchk} checks, {bad} mismatches; a(n) multiplicative n<=2000: {badm} mismatches')
print('    first a(n):', f[1:30])
# (2) completed L-function
c0 = mp.sqrt(23) / (2 * mp.pi)
A = [mp.mpf(f[n]) for n in range(NQ + 1)]
def Lam(s, t=mp.mpf(1), nmax=120):
    tot = mp.mpf(0)
    for n in range(1, nmax):
        if not A[n]: continue
        tot += A[n] * ((c0 / n) ** s * mp.gammainc(s, n * t / c0) + (c0 / n) ** (1 - s) * mp.gammainc(1 - s, n / (t * c0)))
    return tot
for s in (mp.mpf('0.7'), mp.mpc('0.3', '1.7')):
    d = abs(Lam(s, mp.mpf(1)) - Lam(s, mp.mpf('1.3')))
    print(f'(2) |Lambda(s;t=1) - Lambda(s;t=1.3)| at s={s}: {mp.nstr(d, 5)} (0 <=> functional equation, root number +1)')
L1 = Lam(mp.mpf(1)) / c0
alpha = mp.findroot(lambda z: z ** 3 - z - 1, 1.3)
cnf = 2 * mp.pi * mp.log(alpha) / mp.sqrt(23)
print(f'    L(1,rho) = {mp.nstr(L1, 15)}; class number formula 2 pi log(alpha)/sqrt23 = {mp.nstr(cnf, 15)}')
dLam = mp.diff(lambda s: Lam(s), mp.mpf(1))
LpL1 = dLam / Lam(mp.mpf(1)) - mp.log(c0) + mp.euler
env = 2 * dLam / Lam(mp.mpf(1))
print(f'    Lambda\'/Lambda(1) = {mp.nstr(dLam / Lam(mp.mpf(1)), 12)}; L\'/L(1,rho) = {mp.nstr(LpL1, 12)}; GRH envelope sum 1/|rho|^2 = {mp.nstr(env, 10)}')
LpL1 = float(LpL1); env = float(env)
# (3) positivity of Gz + w G_rho and the GRH-scale size of G_rho
negs = {1.0: 0, 0.2: 0, 1.01: 0}
pp = []  # prime powers (n, log p, tr)
for p in primes:
    q = p; kk = 1
    while q <= N:
        t = trk(FT[p], kk); pp.append((q, log(p), t))
        for w in negs:
            if 1 + w * t < 0: negs[w] += 1
        q *= p; kk += 1
pp.sort()
print('(3) prime powers n <= 1e6 with 1 + w tr rho(n) < 0:', negs, '(w <= 1 keeps Gz + w G_rho a nonnegative combination)')
import bisect
ns = [x[0] for x in pp]
S0 = [0.0]; S1 = [0.0]; T0 = [0.0]; T1 = [0.0]
for (n, lp, t) in pp:
    S0.append(S0[-1] + lp * t); S1.append(S1[-1] + lp * t / n)
    T0.append(T0[-1] + lp); T1.append(T1[-1] + lp / n)
def Grho(y):
    i = bisect.bisect_right(ns, y)
    return y * S1[i] - S0[i] + y * LpL1
def Gz(y):
    i = bisect.bisect_right(ns, y)
    return y * T1[i] - T0[i] - y * log(y) + (1 + 0.5772156649015329) * y
ys = [10 ** (2 + 4 * j / 4000) for j in range(4001)]
r = [Grho(y) / sqrt(y) for y in ys]
print(f'    G_rho(y)/sqrt(y), y in [1e2,1e6] (4001 points): max {max(r):+.4f}, min {min(r):+.4f}; envelope {env:.4f}')
for lo, hi in ((1e2, 1e3), (1e3, 1e4), (1e4, 1e5), (1e5, 1e6)):
    seg = [v for y, v in zip(ys, r) if lo <= y <= hi]
    print(f'      [{lo:.0e},{hi:.0e}]: max |G_rho|/sqrt y = {max(abs(v) for v in seg):.4f}')
# windows [Y, Y^1.5] of the level's Gm (base F form: Gz + G_chi(-23) omitted; base-Q form w = 1/5 as Level.cyc of order 3)
for Y in (1e2, 1e3, 1e4):
    seg = [(Gz(y) + 0.2 * Grho(y)) / y ** 0.5 for y in ys if Y <= y <= Y ** 1.5]
    print(f'    window [{Y:.0e}, Y^1.5]: max (Gz + G_rho/5)/sqrt(y) = {max(seg):+.4f}')
