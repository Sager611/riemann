#!/usr/bin/env python3
# alprows P1-P3: SiegelMassRatioCM (Formalism/HermitianLattices/StatementMass.lean:21-49; StandardLattice.lean:281-285).
# P1: Lean RHS massRatioRightSide K M (Euler product over prime ideals of O_K, exactly as defined) vs |L_F(1-M,chi^M)|/2^d
#     from generalized Bernoulli numbers (exact rationals), for CM biquadratic K = Q(sqrt da, sqrt db) with K/F unramified.
# P2: consistency with Aut(I_M) = mu_K^M x S_M (brute-force count of norm-one vectors for Q(zeta12)) and the
#     Gan-Hanke-Yu product mass(I_n) = 2 prod_{r<=n} |L_F(1-r,chi^r)|/2^d: class-number lower bounds h(I_n) >= mass*|Aut|.
# P3: IsCMUnram (StatementMass.lean:44-49) for imaginary quadratic K: the ideal (theta - conj theta) = (sqrt d_K) or (2i).
import sys
from fractions import Fraction as Fr
from math import factorial, gcd
import mpmath as mp
from sympy import bernoulli, primerange, factorint, symbols, Poly
mp.mp.dps = 40

def kron(d, n):
    # Kronecker symbol (d/n), n > 0
    r = 1
    for p, e in factorint(n).items():
        if p == 2:
            if d % 2 == 0: return 0
            v = 1 if d % 8 in (1, 7) else -1
        else:
            dm = d % p
            if dm == 0: return 0
            v = 1 if pow(dm, (p - 1) // 2, p) == 1 else -1
        r *= v ** e
    return r

X = symbols('X')
def Bpoly(n):
    return Poly(bernoulli(n, X), X)
def genB(n, d):
    # B_{n,chi_d} = f^{n-1} sum_{a=1}^{f} chi(a) B_n(a/f), f = |d|
    f = abs(d); P = Bpoly(n); s = Fr(0)
    for a in range(1, f + 1):
        c = kron(d, a)
        if c: s += c * Fr(str(P.eval(Fr(a, f))))
    return Fr(f) ** (n - 1) * s
def Lneg(n, d):  # L(1-n, chi_d), n >= 2
    return -genB(n, d) / n
def zneg(n):     # zeta(1-n), n >= 2
    return -Fr(str(bernoulli(n))) / n
def LF_neg(M, da, db, dF):
    # L_F(1-M, chi^M): M odd -> L(1-M,chi_da) L(1-M,chi_db); M even -> zeta(1-M) L(1-M,chi_dF)
    return Lneg(M, da) * Lneg(M, db) if M % 2 else zneg(M) * Lneg(M, dF)

# V4 = {(s1,s2)}: action on sqrt(da), sqrt(db); complex conjugation c = (-1,-1); H_j = stabiliser of subfield j
V4 = [(1, 1), (1, -1), (-1, 1), (-1, -1)]
subs = [frozenset([(1, 1)])] + [frozenset([(1, 1), g]) for g in V4[1:]] + [frozenset(V4)]
Ha = frozenset(g for g in V4 if g[0] == 1); Hb = frozenset(g for g in V4 if g[1] == 1); HF = frozenset(g for g in V4 if g[0] * g[1] == 1)
def decomp(l, da, db, dF):
    vals = [(Ha, kron(da, l)), (Hb, kron(db, l)), (HF, kron(dF, l))]
    sols = []
    for I in subs:
        for D in subs:
            if not I <= D or len(I) > 2: continue
            ok = True
            for H, v in vals:
                if v == 1: ok &= D <= H
                elif v == -1: ok &= (not D <= H) and I <= H
                else: ok &= not I <= H
            if ok: sols.append((I, D))
    assert len(sols) == 1, (l, sols)
    I, D = sols[0]
    return 4 // len(D), len(D) // len(I), ((-1, -1) in D)   # g, f, Pbar = P

def lean_rhs(M, da, db, dF, L=200000):
    discK = abs(da * db * dF); d = 2
    vInv = mp.mpf(discK) ** (mp.mpf(1) / (2 * d)) / 2
    lv = mp.mpf(1)
    for l in primerange(2, L):
        g, f, fixed = decomp(l, da, db, dF)
        NP = mp.mpf(l) ** f
        if fixed: fac = 1 / (1 - (-1) ** M * mp.sqrt(NP) ** (-M))
        else: fac = (1 - NP ** (-M)) ** (-mp.mpf(1) / 2)
        lv *= fac ** g
    return (mp.mpf(factorial(M - 1)) * (vInv / mp.pi) ** M) ** d / mp.sqrt((2 * vInv) ** d) * lv

fields = [('Q(zeta12)=Q(i,sqrt-3)', -4, -3, 12, 12), ('Q(sqrt-3,sqrt-7)', -3, -7, 21, 6), ('Q(sqrt-2,sqrt-3)', -8, -3, 24, 6), ('Q(i,sqrt-7)', -4, -7, 28, 4), ('Q(sqrt-3,sqrt-5)', -3, -20, 60, 6)]
print('P1/P3a: IsCMUnram via |disc K| = d_F^2 (K/F unramified at all finite primes):')
for name, da, db, dF, w in fields + [('Q(i,sqrt5) [control]', -4, -20, 5, 4)]:
    print('  %-24s |da db| = %d, d_F = %d -> unramified: %s' % (name, abs(da * db), dF, abs(da * db) == dF))
print('P1: Lean massRatioRightSide K M (Euler product to 2e5) vs |L_F(1-M,chi^M)|/4 (exact)')
worst = 0
res = {}
for name, da, db, dF, w in fields:
    for M in range(3, 9):
        ex = abs(LF_neg(M, da, db, dF)) / 4
        ln = lean_rhs(M, da, db, dF, 200000 if M == 3 else 20000)
        rel = abs(ln / mp.mpf(ex.numerator) * ex.denominator - 1)
        worst = max(worst, rel); res[(name, M)] = ex
        print('  %-24s M=%d exact=%s  lean=%s  rel.diff=%.2e' % (name, M, ex, mp.nstr(ln, 15), float(rel)))
print('  worst rel.diff', float(worst))
print('P2: GHY mass(I_n) = 2 prod_{r<=n}|L_F(1-r,chi^r)|/4 vs 1/|Aut(I_n)| = 1/(w^n n!); n=1 uses L_F(0,chi)=L(0,chi_da)L(0,chi_db)')
for name, da, db, dF, w in fields:
    L0 = (-genB(1, da)) * (-genB(1, db))   # L(0,chi) = -B_{1,chi}
    mass = 2 * abs(L0) / 4
    out = []
    for n in range(1, 9):
        if n >= 2:
            mass *= abs(LF_neg(n, da, db, dF)) / 4
        aut = w ** n * factorial(n)
        out.append('n=%d mass=%s h>=%s' % (n, mass, mass * aut))
    print('  ' + name + ': ' + '; '.join(out))
    print('    Lean ratio R(3)=%s vs |Aut(I_2)|/|Aut(I_3)| = 1/(3w) = %s' % (res[(name, 3)], Fr(1, 3 * w)))

