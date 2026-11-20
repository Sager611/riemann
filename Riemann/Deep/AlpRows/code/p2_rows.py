#!/usr/bin/env python3
# alprows P4-P6 (exact unless stated)
# P4: the 1/m + m/x balance in LinkDIV_at's error terms (Formalism/SharpRestriction/GlobalTheoremInputs.lean:83-85,
#     B0 = sizeBudget = f*m*b*x^eps, Nc = budgetNc (Formalism/LocalRingCardinality/DataArith.lean:200-206)),
#     relative to MAIN u (DataParams2.lean:232-234) at the SizesTextAt sizes (CORE.md s5). Log-scale (mpmath).
# P5: idxColLoc g k = |GL_k(Z/g)| / #{unitriangular} (GapAtLeast.lean:56-69) and LevelIndexAtPlace
#     |GL_k(F_q)| = |Stab(phi)| (q^k-1) (StandardLattice.lean:446-451) by brute force over all k x k matrices.
# P6: ClassNumberFormula2 (StandardLattice.lean:453-457): RHS = massI(k) idxPi(g0,g1) (T^k-1)^f must be a positive
#     integer (it is #cls); host E = Q(zeta12), F = Q(sqrt3), f = 2, masses from the GHY product (P2), GBig: 2fk < g_i - 1.
import itertools
from fractions import Fraction as Fr
from math import factorial, prod
import mpmath as mp
from sympy import isprime, bernoulli, Poly, symbols, factorint
mp.mp.dps = 30
print('P4: relative DIV errors vs MAIN ~ m|TT| a log a; single-lift C|TT|a La, per-lift C m (B0+Nc) x^(Cabs eps); Cabs=2, f=n')
for eps in (mp.mpf('0.01'), mp.mpf('0.02')):
    for lx in (12, 24, 48):
        x = mp.mpf(10) ** lx
        n = x ** (2 + 8 * eps); TT = x ** (-eps / 8) * n / 2; T2 = x ** (-eps) * n
        k = x ** (mp.mpf(9) / 2 + 10 * eps); f = n
        loga = mp.log(x) + x ** eps        # a = x b, b = exp(x^eps)
        best = None; row = []
        for mu in (0.3, 0.4, 0.45, 0.5, 0.55, 0.6, 0.7):
            m = x ** mu
            r1 = loga / (m * loga)          # |TT| a La / (m |TT| a log a), La ~ log a
            B0b = f * m * x ** eps          # B0 / b
            Ncb = T2 * k * x ** eps / (n * (n * x ** (-eps)))   # Nc / b
            r2 = m * (B0b + Ncb) * x ** (2 * eps) / (m * TT * x * loga)
            tot = r1 + r2
            row.append('mu=%.2f:%s' % (mu, mp.nstr(mp.log10(tot), 4)))
            if best is None or tot < best[1]: best = (mu, tot)
        print('  eps=%s x=1e%d log10(rel err): %s ; best mu=%.2f' % (eps, lx, ' '.join(row), best[0]))
print('P5: brute force over all k x k matrices mod q (q prime)')
def det_mod(A, q):
    k = len(A); A = [r[:] for r in A]; d = 1
    for c in range(k):
        p = next((r for r in range(c, k) if A[r][c] % q), None)
        if p is None: return 0
        if p != c: A[c], A[p] = A[p], A[c]; d = -d
        d = d * A[c][c] % q; inv = pow(A[c][c], q - 2, q)
        for r in range(c + 1, k):
            t = A[r][c] * inv % q
            A[r] = [(A[r][j] - t * A[c][j]) % q for j in range(k)]
    return d % q
for q in (2, 3, 5):
    for k in (1, 2, 3):
        if q == 5 and k == 3: continue
        GL = 0; U = 0; stab = 0
        for ent in itertools.product(range(q), repeat=k * k):
            A = [list(ent[i * k:(i + 1) * k]) for i in range(k)]
            if det_mod(A, q) == 0: continue
            GL += 1
            if all(A[i][j] == 0 for i in range(k) for j in range(k) if j < i) and all(A[i][i] == 1 for i in range(k)): U += 1
            # phi = e_1 (row vector); vecMul phi A = first row of A
            if A[0] == [1] + [0] * (k - 1): stab += 1
        print('  q=%d k=%d |GL|=%d idxColLoc=%s prod(q^i-1)=%d  |GL|==|Stab|(q^k-1): %s' % (q, k, GL, Fr(GL, U), prod(q ** i - 1 for i in range(1, k + 1)), GL == stab * (q ** k - 1)))
print('P6: ClassNumberFormula2 integrality, E = Q(zeta12), F = Q(sqrt3), f = 2')
X = symbols('X')
def kron(d, nn):
    r = 1
    for p, e in factorint(nn).items():
        if p == 2:
            if d % 2 == 0: return 0
            v = 1 if d % 8 in (1, 7) else -1
        else:
            dm = d % p
            if dm == 0: return 0
            v = 1 if pow(dm, (p - 1) // 2, p) == 1 else -1
        r *= v ** e
    return r
def genB(nn, d):
    f_ = abs(d); P = Poly(bernoulli(nn, X), X); s = Fr(0)
    for a in range(1, f_ + 1):
        c = kron(d, a)
        if c: s += c * Fr(str(P.eval(Fr(a, f_))))
    return Fr(f_) ** (nn - 1) * s
def LF(M):
    return (genB(M, -4) / M) * (genB(M, -3) / M) if M % 2 else (Fr(str(bernoulli(M))) / M) * (genB(M, 12) / M)
mass = {1: Fr(1, 12)}
for M in range(2, 7): mass[M] = mass[M - 1] * abs(LF(M)) / 4
splitF = [g for g in range(5, 400) if isprime(g) and g % 12 in (1, 11)]
bad = 0; tot = 0; ex = []
for k in (3, 4, 5, 6):
    for g0 in splitF[:12]:
        for g1 in splitF[:12]:
            if g1 <= g0: continue
            if not (2 * 2 * k < g0 - 1 or 2 * 2 * k < g1 - 1): continue
            for T in splitF[:6]:
                rhs = mass[k] * (prod(g0 ** i - 1 for i in range(1, k + 1)) * prod(g1 ** i - 1 for i in range(1, k + 1))) ** 2 * (T ** k - 1) ** 2
                tot += 1
                if rhs.denominator != 1:
                    bad += 1
                    if len(ex) < 4: ex.append((k, g0, g1, T, str(rhs)[:40]))
print('  masses', {M: str(v) for M, v in mass.items()})
print('  cases', tot, 'non-integral RHS', bad, 'examples', ex)

