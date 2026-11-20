# alpext P1: the colour balance (AskMainBalance, Formalism/SizeWindowBrackets/FirstPrimeNeSecondPrime.lean:84-101) with the
# cubic class-group character psi of F = Q(sqrt-23) (L(s,psi) = L(s,rho), rho the 2-dim rep of S3 = Gal(H/Q)) in place of chi.
# Idea copied from Riemann/Deep/AlpChain/code/colour_balance_exact.py (abelian E0 x H, formal log symbols, exact Fractions).
# Host group for psi: G = S3 x C2 (C2 = Gal(K/Q), K = Q(sqrt-7), d = 7 = 7 mod 16), Phi = {e = 0} (split in K).
# A: partner t g = cconj + g + topPl t (PartnersToFactorsAt.lean:24-26) must be well defined on Frobenius CLASSES.
# B: orthogonality with linear characters only (all AddChar of an abelian group) vs with Irr(G) (incl. rho).
# C: Chebotarev independence: any CM host E has E cap H inside F; class functions of the host have zero rho-component.
import itertools, random
from fractions import Fraction as Fr
from math import isqrt
S3 = list(itertools.permutations(range(3)))
comp = lambda a, b: tuple(a[b[i]] for i in range(3))
inv = lambda a: tuple(sorted(range(3), key=lambda i: a[i]))
def sgn(s):
    return 1 if sum(1 for i in range(3) for j in range(i+1, 3) if s[i] > s[j]) % 2 == 0 else -1
def trrho(s):  # character of the 2-dim irreducible rep of S3: fixed points - 1
    return sum(1 for i in range(3) if s[i] == i) - 1
G = [(s, e) for s in S3 for e in (0, 1)]
gmul = lambda g, h: (comp(g[0], h[0]), (g[1] + h[1]) % 2)
ginv = lambda g: (inv(g[0]), g[1])
def gpow(g, k):
    r = ((0, 1, 2), 0)
    for _ in range(k): r = gmul(r, g)
    return r
cls = lambda g: frozenset(gmul(gmul(h, g), ginv(h)) for h in G)
CL = sorted(set(cls(g) for g in G), key=lambda C: (min(C)[1], len(C), sorted(C)))
tau_c = (0, 2, 1)            # complex conjugation on the roots (P0: real root fixed)
cconj = (tau_c, 1)           # acts nontrivially on K and on H
print('A. classes of G = S3 x C2:', len(CL), 'sizes', [len(C) for C in CL])
print('   centre of G:', [g for g in G if all(gmul(g, h) == gmul(h, g) for h in G)], '; cconj central:', all(gmul(cconj, h) == gmul(h, cconj) for h in G))
rnd = random.Random(1); bad = 0; tot = 0
for trial in range(50):
    wC = {C: Fr(rnd.randint(10, 99), 10) for C in CL}
    w = lambda g: wC[cls(g)]
    for C in CL:
        vals = set(w(gmul(cconj, g)) for g in C); tot += 1
        if len(vals) > 1: bad += 1
print(f'   class weights w: g -> w(cconj*g) non-constant on {bad}/{tot} (class, trial) pairs => ASK depends on the Frobenius representative')
# abelian control: replace S3 by C6 (abelian, as in the archive's E0 x H): always constant
C6 = [(k, e) for k in range(6) for e in (0, 1)]
m6 = lambda g, h: ((g[0] + h[0]) % 6, (g[1] + h[1]) % 2)
bad6 = 0
for trial in range(50):
    wv = {g: Fr(rnd.randint(10, 99), 10) for g in C6}
    for g in C6:
        if len(set([wv[m6((3, 1), g)]])) > 1: bad6 += 1
print(f'   abelian control (C6 x C2, every element its own class): non-constant {bad6}')
# B. orthogonality
lin = [(a, b) for a in (0, 1) for b in (0, 1)]
def lin_char(c, g): return (sgn(g[0]) if c[0] else 1) * ((-1) ** g[1] if c[1] else 1)
irr = [('lin', c) for c in lin] + [('rho', b) for b in (0, 1)]
def chi(c, g):
    if c[0] == 'lin': return lin_char(c[1], g)
    return trrho(g[0]) * ((-1) ** g[1] if c[1] else 1)
# rank of linear characters on class functions
import sympy as sp
Mlin = sp.Matrix([[lin_char(c, next(iter(C))) for C in CL] for c in lin])
Mirr = sp.Matrix([[chi(c, next(iter(C))) for C in CL] for c in irr])
print('B. class functions: dim', len(CL), '; rank(linear chars) =', Mlin.rank(), '; rank(Irr incl. rho) =', Mirr.rank())
# real primes: Frobenius in S3 from x^3-x-1 mod p (labelled permutation via roots in F_p / F_p^2 is not needed:
# only the class matters for class functions; pick the canonical representative of the class)
def leg(a, p): r = pow(a % p, (p - 1) // 2, p); return -1 if r == p - 1 else r
def pmulmod(A, B, p):
    r = [0] * 5
    for i, ai in enumerate(A):
        for j, bj in enumerate(B): r[i + j] += ai * bj
    for k in (4, 3):
        if r[k]: r[k - 2] += r[k]; r[k - 3] += r[k]; r[k] = 0
    return [v % p for v in r[:3]]
def xpow(p):
    R = [1, 0, 0]; B = [0, 1, 0]; e = p
    while e:
        if e & 1: R = pmulmod(R, B, p)
        B = pmulmod(B, B, p); e >>= 1
    return R
def frob(p):
    L = leg(-23, p)
    s = (0, 2, 1) if L == -1 else ((0, 1, 2) if xpow(p) == [0, 1, 0] else (1, 2, 0))
    return (s, 0 if leg(-7, p) == 1 else 1)
Ny = 3000
P = [p for p in range(3, Ny + 1) if all(p % r for r in range(2, isqrt(p) + 1)) and p not in (7, 23)]
FR = {p: frob(p) for p in P}
def ramp_by(fun, y):  # formal sum_{p^k<=y, p unramified} log p * fun(Frob_p^k) * (y/p^k - 1), as {p: Fraction}
    out = {}
    for p in P:
        k = 1
        while p ** k <= y:
            v = Fr(fun(gpow(FR[p], k))) * (y / Fr(p ** k) - 1)
            if v: out[p] = out.get(p, 0) + v
            k += 1
    return out
def lincomb(terms):
    out = {}
    for c, d in terms:
        for k, v in d.items(): out[k] = out.get(k, 0) + c * v
    return {k: v for k, v in out.items() if v != 0}
ok_irr = 0; res_lin = []
for y in (Fr(997, 1), Fr(2999, 2), Fr(3000, 1)):
    for C in CL:
        g0 = next(iter(C))
        if g0[1] == 0: continue  # only classes outside Phi (inert in K) enter ASK
        lhs = ramp_by(lambda h: 1 if h in C else 0, y)
        rhs_irr = lincomb([(Fr(len(C), len(G)) * chi(c, g0), ramp_by(lambda h, c=c: chi(c, h), y)) for c in irr])
        rhs_lin = lincomb([(Fr(len(C), len(G)) * lin_char(c, g0), ramp_by(lambda h, c=c: lin_char(c, h), y)) for c in lin])
        d_irr = lincomb([(1, lhs), (-1, rhs_irr)]); d_lin = lincomb([(1, lhs), (-1, rhs_lin)])
        ok_irr += (len(d_irr) == 0)
        res_lin.append((float(y), len(C), trrho(g0[0]), len(d_lin)))
print('   real primes 3 <= p <= 3000 (p != 7, 23): identity with Irr(G) exact in', ok_irr, '/ 9 (y, class) cases')
print('   with linear characters only, nonzero log-p symbols per (y, |C|, tr rho(C)):', res_lin)
# C. Chebotarev independence for CM hosts: E cap H inside F => rho averages to 0 on every fibre of G_E x_{Gal(F/Q)} S3
fib = {+1: [s for s in S3 if sgn(s) == 1], -1: [s for s in S3 if sgn(s) == -1]}
print('C. mean of tr rho on the fibres over Gal(F/Q): +1 ->', Fr(sum(trrho(s) for s in fib[1]), 3), '; -1 ->', Fr(sum(trrho(s) for s in fib[-1]), 3), '; whole S3 ->', Fr(sum(trrho(s) for s in S3), 6))
# empirical: host Q(zeta_23) (abelian CM, contains F): mean of tr rho(Frob_p) per residue class a = p mod 23, p <= 1e6
N = 10 ** 6
sv = bytearray([1]) * (N + 1); sv[0] = sv[1] = 0
for i in range(2, isqrt(N) + 1):
    if sv[i]: sv[i * i::i] = bytearray(len(range(i * i, N + 1, i)))
S = {a: [0, 0] for a in range(1, 23)}; chiS = [0, 0]
for p in range(5, N + 1):
    if not sv[p] or p == 23: continue
    t = trrho(frob(p)[0]); a = p % 23; S[a][0] += t; S[a][1] += 1
print('   mean tr rho(Frob_p) by p mod 23, p <= 1e6 (n = count):')
print('   ' + '; '.join(f'{a}:{S[a][0]/S[a][1]:+.4f}(n={S[a][1]})' for a in range(1, 23)))
mx = max(abs(S[a][0] / S[a][1]) * S[a][1] ** 0.5 for a in S)
print(f'   max over classes |mean| * sqrt(n) = {mx:.3f} (O(1) = no correlation; a Dirichlet-type target would give |mean| ~ 1)')
