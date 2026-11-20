"""gl2gapaudit: independent checks of the GL(2) order-6 certificate data (no import of cert5/exact5/exact6).
(a) rationalisation of the float form (SRC) vs the certified JSON; a1 = 1, a6 > 0.
(b) H exact (Fraction) and H by a different formula: H = -(1/Lam) int_0^Lam p(v^2) dv (mpmath Gauss-Legendre, 60 digits).
(c) uniform-measure identity sum_{k,l} C(j,k)C(j,l)(-1)^{k+l}2^{k+l+1}/(k+l+1) = 2/(2j+1), j <= 6 (exact).
(d) falsification of (V6) on h(u) = u(1-u)u^m, m < M, on [0,1] scaled to [0,g]: exact rational matrices, exact LDL^T signs,
    min generalized eigenvalue of (Q_g, g*G0) at 60 digits.
usage: indep.py SRCjson CERTjson M"""
import sys, json
from fractions import Fraction as Fr
from math import comb
import mpmath as mp
mp.mp.dps = 60
src = json.load(open(sys.argv[1]))[0]; cert = json.load(open(sys.argv[2])); M = int(sys.argv[3])
af = src['a6']; phf = {k.strip('()').replace(', ', ''): v for k, v in src['phi6'].items()}
a = [Fr(x).limit_denominator(10**12) for x in af]; ph = {k: Fr(v).limit_denominator(10**12) for k, v in phf.items()}
ac = [Fr(x) for x in cert['a']]; phc = {k: Fr(v) for k, v in cert['phi'].items()}
print('(a) a == cert a:', a == ac, '; phi == cert phi:', ph == phc, '; a1 =', ac[1], '; a6 = %.6g > 0:' % float(ac[6]), ac[6] > 0)
print('    max |float - rational| = %.2e' % max([abs(float(x) - y) for x, y in zip(ac, af)] + [abs(float(phc[k]) - phf[k]) for k in phf]))
lam = Fr(218, 125)
H = -sum(ac[j] * lam ** (2 * j) / (2 * j + 1) for j in range(7))
p = lambda x: sum(mp.mpf(ac[j].numerator) / ac[j].denominator * x ** j for j in range(7))
L = mp.mpf(218) / 125
Hq = -mp.quad(lambda v: p(v * v), [0, L], method='gauss-legendre') / L
print('(b) H exact = %s ~ %s; H by quadrature = %s; diff = %.2e' % ('%d/%d digits' % (len(str(H.numerator)), len(str(H.denominator))), mp.nstr(mp.mpf(H.numerator) / H.denominator, 25), mp.nstr(Hq, 25), float(abs(Hq - mp.mpf(H.numerator) / H.denominator))))
# H(A) for A slightly above 218/125: continuity margin
for A in (Fr(1744, 1000), Fr(17445, 10000), Fr(1745, 1000), Fr(1746, 1000)):
    print('    H(%s) = %.6e' % (A, float(-sum(ac[j] * A ** (2 * j) / (2 * j + 1) for j in range(7)))))
ok = all(sum(comb(j, k) * comb(j, l) * Fr((-1) ** (k + l) * 2 ** (k + l + 1), k + l + 1) for k in range(j + 1) for l in range(j + 1)) == Fr(2, 2 * j + 1) for j in range(7))
print('(c) binomial recombination = 2/(2j+1) for j <= 6:', ok)
# (d) polynomial falsification
def pmul(P, Q):
    R = [Fr(0)] * (len(P) + len(Q) - 1)
    for i, x in enumerate(P):
        if x:
            for j, y in enumerate(Q):
                R[i + j] += x * y
    return R
def pder(P, k):
    for _ in range(k):
        P = [i * P[i] for i in range(1, len(P))] or [Fr(0)]
    return P
def pint01(P):
    return sum(c / (i + 1) for i, c in enumerate(P))
def pev(P, u):
    return sum(c * u ** i for i, c in enumerate(P))
basis = [[Fr(0)] * (m + 1) + [Fr(1), Fr(-1)] for m in range(M)]  # u^{m+1} - u^{m+2}
D = [[pder(b, j) for b in basis] for j in range(7)]
G = [[[pint01(pmul(D[j][m], D[j][n])) for n in range(M)] for m in range(M)] for j in range(7)]
pairs = {'12': (1, 2), '14': (1, 4), '32': (3, 2), '34': (3, 4), '52': (5, 2), '54': (5, 4)}
B = {}
for k, (pp, qq) in pairs.items():
    B[k] = [[(pev(D[pp][m], 1) * pev(D[qq][n], 1) + pev(D[pp][n], 1) * pev(D[qq][m], 1) - pev(D[pp][m], 0) * pev(D[qq][n], 0) - pev(D[pp][n], 0) * pev(D[qq][m], 0)) / 2 for n in range(M)] for m in range(M)]
def ldl_signs(Q):
    A = [row[:] for row in Q]; n = len(A); piv = []
    for i in range(n):
        d = A[i][i]; piv.append(d)
        if d == 0:
            return piv
        for r in range(i + 1, n):
            f = A[r][i] / d
            if f:
                for c in range(i, n):
                    A[r][c] -= f * A[i][c]
    return piv
gs = [Fr(1, 100), Fr(1, 10), Fr(1, 2), Fr(1), Fr(3, 2), Fr(2), Fr(2108718, 10**6), Fr(5, 2), Fr(3), Fr(31, 10), Fr(314, 100), Fr(3141592, 10**6), Fr(31415926, 10**7)]
for g in gs:
    Q = [[sum(ac[j] * g ** (1 - 2 * j) * G[j][m][n] for j in range(7)) + sum(phc[k] * g ** (-pairs[k][0] - pairs[k][1]) * B[k][m][n] for k in pairs) for n in range(M)] for m in range(M)]
    piv = ldl_signs(Q)
    Qm = mp.matrix([[mp.mpf(x.numerator) / x.denominator for x in row] for row in Q])
    Gm = mp.matrix([[mp.mpf((g * x).numerator) / (g * x).denominator for x in row] for row in G[0]])
    Lc = mp.cholesky(Gm); Li = mp.inverse(Lc)
    S = Li * Qm * Li.T; S = (S + S.T) / 2
    ev = mp.eigsy(S, eigvals_only=True)
    print('(d) g = %-12s exact LDL: all pivots > 0: %s (min pivot %.3e); min eig of form/||f||^2 = %s' % (g, all(x > 0 for x in piv), float(min(piv)), mp.nstr(min(ev), 8)), flush=True)
