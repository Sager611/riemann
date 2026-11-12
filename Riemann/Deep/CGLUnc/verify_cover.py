"""Arb verification of the covering certificate cover_S_P.json: M(lambda) is positive definite (LDL^T pivots > 0),
hence for every window c' in R^16 some library configuration X_i (density nu0 = P/S) has
nu0 - H(c) - e(X_i; c) > 0, so no TEN-type certificate on any 16-harmonic window proves kappa > nu0.
All quantities are recomputed from closed forms in ball arithmetic (positions and lambda rationalised)."""
import sys, json
from fractions import Fraction
from flint import arb, acb, ctx

ctx.prec = 160
d = json.load(open(sys.argv[1]))
J = int(sys.argv[2]) if len(sys.argv) > 2 else 16
NJ = J + 1
S, P = int(d["S"]), int(d["P"])
pi = arb.pi(); s2 = arb(2).sqrt(); t0 = 1 / s2
om = [s2] + [2 * pi * j for j in range(1, J + 1)]
sin2 = t0.sin() ** 2
HMT = arb(3) / 2 - (t0.cos() / t0.sin()) / s2
W = [arb(0)] + [(arb(1) / 2 - 1 / (4 * pi ** 2 * j ** 2)) / (2 * sin2) for j in range(1, J + 1)]


def A(q):
    return arb(q.numerator) / arb(q.denominator)


def I(a, b, Om, ph):
    if Om is None:
        return (b - a) * ph.cos()
    return ((Om * b + ph).sin() - (Om * a + ph).sin()) / Om


def Tsym(al):
    lo, hi = arb(-1) / 2, arb(1) / 2 - al
    T = [[None] * NJ for _ in range(NJ)]
    for j in range(NJ):
        for k in range(NJ):
            dlt = None if j == k else om[j] - om[k]
            T[j][k] = (I(lo, hi, dlt, -om[k] * al) + I(lo, hi, om[j] + om[k], om[k] * al)) / 2
    return [[(T[j][k] + T[k][j]) / 2 for k in range(NJ)] for j in range(NJ)]


Ts = [Tsym(arb(m) / S) for m in range(S)]
lam = [Fraction(v).limit_denominator(10 ** 12) for v in d["lam"]]
lam = [max(l, Fraction(0)) for l in lam]
tot = sum(lam); lam = [l / tot for l in lam]
Qb = [[arb(0)] * NJ for _ in range(NJ)]
for li, xs in zip(lam, d["x"]):
    if li == 0:
        continue
    xr = [A(Fraction(round(v * 10 ** 9), 10 ** 9)) for v in xs]
    assert len(xr) == P
    wts = [arb(P * P)] + [arb(0)] * (S - 1)
    for m in range(1, S):
        mu = acb(0)
        for x in xr:
            mu += (2 * pi * acb(0, 1) * m * x / S).exp()
        wts[m] = 2 * abs(mu) ** 2
    L = A(li)
    for j in range(NJ):
        for k in range(NJ):
            acc = arb(0)
            for m in range(S):
                acc += wts[m] * Ts[m][j][k]
            Qb[j][k] += L * acc / (2 * sin2 * S * S)
nu0 = arb(P) / S
M = [[None] * NJ for _ in range(NJ)]
M[0][0] = 2 * nu0 - HMT - Qb[0][0]
for j in range(1, NJ):
    M[0][j] = -Qb[0][j]; M[j][0] = -Qb[j][0]
    for k in range(1, NJ):
        M[j][k] = (W[j] if j == k else arb(0)) - Qb[j][k]
# LDL^T with ball arithmetic
n = NJ; Lm = [[arb(0)] * n for _ in range(n)]; D = [arb(0)] * n
ok = True
for i in range(n):
    s = M[i][i]
    for k in range(i):
        s -= Lm[i][k] ** 2 * D[k]
    D[i] = s
    if not (D[i] > 0):
        ok = False
        print("pivot", i, "not certified positive:", D[i].str(10)); break
    for j in range(i + 1, n):
        s = M[j][i]
        for k in range(i):
            s -= Lm[j][k] * Lm[i][k] * D[k]
        Lm[j][i] = s / D[i]
print("J =", J, " nu0 =", Fraction(P, S), "=", nu0.str(10), " configs used:", sum(1 for l in lam if l > 0))
print("min pivot:", min((x for x in D if x != 0), key=lambda z: float(z.mid())).str(10) if ok else "-")
print("M positive definite (certified):", ok)
