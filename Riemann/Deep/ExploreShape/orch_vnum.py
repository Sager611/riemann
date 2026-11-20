# orchestrator 49: independent numeric test of explore53's per-gap inequality (V) on polynomial bases (exact rationals, then mpmath)
import json, sys
from fractions import Fraction as Fr
import mpmath as mp
mp.mp.dps = 60
J = json.load(open('Riemann/Deep/ExploreShape/logs/cert4_k3.34.json'))
a = [Fr(s) for s in J['a']]; p12 = Fr(J['phi12']); p32 = Fr(J['phi32'])
M = int(sys.argv[1]) if len(sys.argv) > 1 else 20
def pmul(p, q):
    r = [Fr(0)] * (len(p) + len(q) - 1)
    for i, x in enumerate(p):
        if x:
            for j, y in enumerate(q):
                r[i + j] += x * y
    return r
def pder(p): return [i * p[i] for i in range(1, len(p))] or [Fr(0)]
def pint01(p): return sum(c / (i + 1) for i, c in enumerate(p))
def pval(p, x): return sum(c * x**i for i, c in enumerate(p))
# shifted Legendre P_m(2u-1) by recurrence
P = [[Fr(1)], [Fr(-1), Fr(2)]]
while len(P) < M:
    m = len(P) - 1
    t = pmul([Fr(-1), Fr(2)], P[m])
    nxt = [((2 * m + 1) * (t[i] if i < len(t) else 0) - m * (P[m - 1][i] if i < len(P[m - 1]) else 0)) / (m + 1) for i in range(len(t))]
    P.append(nxt)
H = [pmul([Fr(0), Fr(1), Fr(-1)], P[m]) for m in range(M)]          # u(1-u) P_m(2u-1)
D = [[h] for h in H]
for h in D:
    for j in range(4): h.append(pder(h[-1]))
I = [[[pint01(pmul(D[m][j], D[n][j])) for n in range(M)] for m in range(M)] for j in range(5)]
B12 = [[(pval(D[m][1], 1) * pval(D[n][2], 1) + pval(D[n][1], 1) * pval(D[m][2], 1) - pval(D[m][1], 0) * pval(D[n][2], 0) - pval(D[n][1], 0) * pval(D[m][2], 0)) / 2 for n in range(M)] for m in range(M)]
B32 = [[(pval(D[m][3], 1) * pval(D[n][2], 1) + pval(D[n][3], 1) * pval(D[m][2], 1) - pval(D[m][3], 0) * pval(D[n][2], 0) - pval(D[n][3], 0) * pval(D[m][2], 0)) / 2 for n in range(M)] for m in range(M)]
G = mp.matrix([[mp.mpf(I[0][m][n].numerator) / I[0][m][n].denominator for n in range(M)] for m in range(M)])
L = mp.cholesky(G); Li = mp.inverse(L)
out = []
for g in [mp.pi, mp.mpf('3.1'), mp.mpf('3.0'), mp.mpf('2.7707'), mp.mpf('2.5'), mp.mpf('2.0'), mp.mpf('1.0'), mp.mpf('0.3')]:
    A = mp.matrix(M, M)
    for m in range(M):
        for n in range(M):
            s = sum(mp.mpf(a[j].numerator) / a[j].denominator * g**(1 - 2 * j) * mp.mpf(I[j][m][n].numerator) / I[j][m][n].denominator for j in range(5))
            s += mp.mpf(p12.numerator) / p12.denominator * g**(-3) * mp.mpf(B12[m][n].numerator) / B12[m][n].denominator
            s += mp.mpf(p32.numerator) / p32.denominator * g**(-5) * mp.mpf(B32[m][n].numerator) / B32[m][n].denominator
            A[m, n] = s
    # form relative to ||f||^2 = g * ||h||^2_{L2(0,1)}
    C = Li * A * Li.T
    C = (C + C.T) / 2
    ev = mp.eigsy(C, eigvals_only=True)
    mn = min(ev) / g
    out.append((float(g), mp.nstr(mn, 12)))
    print('g = %.6f  min eigenvalue of the form / ||f||^2 = %s' % (float(g), mp.nstr(mn, 12)), flush=True)
print('M =', M, 'NEGATIVE' if any(mp.mpf(v) < 0 for _, v in out) else 'all >= 0')
