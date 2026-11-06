#!/usr/bin/env python3
"""degree2mol certificate.  Rigorous (ball-arithmetic) evaluation of Conrey's mollified-moment constant
    c(P,Q,R) = 1 + (1/xi) * int_0^1 int_0^1 e^{2Ry} [ Q(y) P'(x) + xi (R Q(y) + Q'(y)) P(x) ]^2 dx dy
(Conrey 1989 Thm 2; for GL(2): Bernard 2015 Thm 5 / Andersen-Thorner 2021 eq. (2.3)-(2.4), with xi = nu/2)
and of kappa = 1 - log(c)/R, for explicit polynomials P, Q with exact rational coefficients.
All integrals are done exactly: the x-integral is polynomial algebra over Q, the y-integral uses
int_0^1 e^{ay} G(y) dy = sum_k (-1)^k (e^a G^(k)(1) - G^(k)(0)) / a^(k+1)   (finite sum, G polynomial),
so the only inexact step is one exp and one log in Arb (4000 bits), with outward rounding.
Run: uv run --no-project --with python-flint --with mpmath python3 cert.py
"""
import json, os, sys
from fractions import Fraction
from flint import fmpq, fmpq_poly, arb, ctx
import mpmath as mp
ctx.prec = 4000   # the exact sums reach ~1e300 before cancelling to O(1) when R is small (family_PGL2: R = 1.48)
HERE = os.path.dirname(os.path.abspath(__file__))
def F(s):
    fr = Fraction(s); return fmpq(fr.numerator, fr.denominator)
X = fmpq_poly([0, 1]); ONE = fmpq_poly([1])
def cheb(m):
    Ts = [ONE, X]
    for k in range(2, m + 1): Ts.append(2 * X * Ts[-1] - Ts[-2])
    return Ts[m]
def Q_cheb(h):            # Q(y) = 1 + sum h_n (T_{2n-1}(1-2y) - 1)
    z = 1 - 2 * X; Q = ONE
    for n, hn in enumerate(h, start=1): Q += hn * (cheb(2 * n - 1)(z) - 1)
    return Q
def Q_mono(h):            # AT21 basis: Q(y) = 1 + sum h_n ((1-2y)^{2n-1} - 1)
    z = 1 - 2 * X; Q = ONE
    for n, hn in enumerate(h, start=1): Q += hn * (z ** (2 * n - 1) - 1)
    return Q
def P_sinh(lam, K=12):    # P = p/p(1), p = Taylor polynomial of sinh(lam x), degree 2K+1
    p = fmpq_poly([0]); fact = 1
    for j in range(1, 2 * K + 2):
        fact *= j
        if j % 2 == 1: p += fmpq_poly([0] * j + [lam ** j / fact])
    return p / p(fmpq(1))
def int01(p):
    I = p.integral(); return I(fmpq(1)) - I(fmpq(0))
def int_exp(G, a):        # exact pair (Xc, Yc): int_0^1 e^{a y} G(y) dy = e^a Xc - Yc
    Xc = fmpq(0); Yc = fmpq(0); D = G; k = 0; ap = a
    while not D.is_zero():
        s = 1 if k % 2 == 0 else -1
        Xc += s * D(fmpq(1)) / ap; Yc += s * D(fmpq(0)) / ap
        D = D.derivative(); ap *= a; k += 1
    return Xc, Yc
def check_admissible(P, Q):
    dQ = Q.derivative()
    ok = P(fmpq(0)) == 0 and P(fmpq(1)) == 1 and Q(fmpq(0)) == 1 and (dQ - dQ(1 - X)).is_zero()
    assert ok, "P(0)=0, P(1)=1, Q(0)=1, Q'(y)=Q'(1-y) violated"
    return ok
def c_kappa(P, Q, R, xi):
    check_admissible(P, Q)
    dQ = Q.derivative(); U = R * Q + dQ; dP = P.derivative()
    a1, a2, a3 = int01(dP * dP), int01(P * dP), int01(P * P)
    G = a1 * Q * Q + 2 * xi * a2 * Q * U + xi * xi * a3 * U * U
    Xc, Yc = int_exp(G, 2 * R)
    c = 1 + (arb(Xc) * arb(2 * R).exp() - arb(Yc)) / arb(xi)
    return c, 1 - c.log() / arb(R)
def W01(Q, R):            # W0 = int e^{2Ry} Q^2, W1 = int e^{2Ry} (RQ+Q')^2 (Arb)
    U = R * Q + Q.derivative(); e = arb(2 * R).exp(); out = []
    for G in (Q * Q, U * U):
        Xc, Yc = int_exp(G, 2 * R); out.append(arb(Xc) * e - arb(Yc))
    return out
def lam_opt(Q, R, xi):    # Conrey's optimal P = sinh(lam x)/sinh(lam), lam = xi sqrt(W1/W0); rounded to 1e-8
    W0, W1 = W01(Q, R); l = float((arb(xi) * (W1 / W0).sqrt()).mid())
    return F(str(round(l, 8)))
def mpf(a): return mp.mpf(a.mid().str(60, radius=False))
def cross_checks(P, Q, R, xi, c):
    """(i) Andersen-Thorner closed form inf_P c; (ii) the general-P formula (kappa: general) of
    arXiv 2508.11108 l.879-881 with K(S) by mpmath quadrature. Both are consistency checks only."""
    mp.mp.dps = 140   # Q has monomial coefficients up to ~1e45 (T_59(1-2y)), so evaluate S with many digits
    W0, W1 = W01(Q, R); W0, W1 = mpf(W0), mpf(W1); Rm = mp.mpf(Fraction(str(R))); xim = mp.mpf(Fraction(str(xi)))
    Q1 = mp.mpf(Fraction(str(Q(fmpq(1)))))
    cAT = (1 + mp.e ** (2 * Rm) * Q1 ** 2) / 2 + mp.sqrt(W0 * W1) / mp.tanh(xim * mp.sqrt(W1 / W0))
    B = mp.mpf(Fraction(str(int01(P * P)))); Cc = mp.mpf(Fraction(str(int01(P.derivative() ** 2))))
    beta = 1 + Q1; c0 = Cc / xim - xim * B * Rm ** 2; c1 = 4 * xim * B * Rm ** 2
    qc = [mp.mpf(Fraction(str(x))) for x in Q.coeffs()]; dqc = [mp.mpf(Fraction(str(x))) for x in Q.derivative().coeffs()]
    S = lambda t: mp.polyval(qc[::-1], mp.mpf(1) / 2 + t / (2 * Rm))
    dS = lambda t: mp.polyval(dqc[::-1], mp.mpf(1) / 2 + t / (2 * Rm)) / (2 * Rm)
    K = mp.quad(lambda t: mp.e ** t * (c0 * S(t) ** 2 + c1 * dS(t) ** 2) + mp.e ** (-t) * (c0 * (beta - S(t)) ** 2 + c1 * dS(t) ** 2),
                mp.linspace(0, Rm, 9))
    cG = (1 + mp.e ** (2 * Rm) * (beta - 1) ** 2) / 2 + c1 * mp.e ** Rm / (2 * Rm) * (mp.e ** Rm * (1 - beta) ** 2 - mp.e ** (-Rm)) / 2 + mp.e ** Rm / (2 * Rm) * K
    cm = mpf(c)
    return float((cAT - cm) / cm), float((cG - cm) / cm)

def run():
    prm = json.load(open(os.path.join(HERE, 'params.json')))
    print("== A. Reproductions ==")
    at21 = json.load(open(os.path.join(HERE, 'at21_tables.json')))
    for key, nu, claim in (('KS', '25/128', '0.0696'), ('RSC', '1/4', '0.0896')):
        R = F(at21[key]['R']); h = [F(s) for s in at21[key]['h']]; xi = F(nu) / 2
        Q = Q_mono(h); lam = lam_opt(Q, R, xi); P = P_sinh(lam)
        c, k = c_kappa(P, Q, R, xi)
        print(f"AT21 table {key:3s} (nu={nu}, xi={xi}): R={float(R):.6f} lam={float(lam):.8f} kappa in {k.str(12)}  [AT21 claim kappa_f > {claim}]")
    for key, claim in (('paper_Bernard_uncond', '0.0632'), ('paper_RSC', '0.0854')):
        d = prm[key]; R = F(d['R']); xi = F(d['xi']); h = [F(s) for s in d['h']]
        h[-1] = fmpq(1, 2) - sum(h[:-1], fmpq(0))          # enforce Q(1)=0 exactly (beta=1)
        Q = Q_cheb(h); assert Q(fmpq(1)) == 0
        c, k = c_kappa(X, Q, R, xi)
        print(f"2508.11108 P=x, beta=1, R=sqrt(3/5)/xi ~ {float(R):.6f}, xi={xi}: kappa in {k.str(12)}  [paper table {claim}]")
    print("\n== B. New certified values (P, Q jointly optimised; Q(1) free) ==")
    res = {}
    for key, label in (('AT21_uncond', 'unconditional, AT21 length nu < 25/128'), ('RSC', 'Ramanujan (RSC), nu < 1/4'),
                       ('Bernard_uncond', 'unconditional, Bernard length nu < 5/27'),
                       ('family_PGL3', 'arXiv 2607.00282 Cor 6.1 r=3: kappa(theta), theta -> 1/6-'),
                       ('family_PGL2', 'arXiv 2607.00282 Cor 6.1 r=2: kappa(theta), theta -> 1/2-'),
                       ('family_PGL2_simple', 'arXiv 2607.00282 Cor simplpergl2: simple zeros, linear Q, theta -> 1/2-')):
        d = prm[key]; R = F(d['R']); xi = F(d['xi']); lam = F(d['lam']); h = [F(s) for s in d['h']]
        bound = {'AT21_uncond': fmpq(25, 256), 'RSC': fmpq(1, 8), 'Bernard_uncond': fmpq(5, 54),   # nu/2 for nu < 25/128, 1/4, 5/27
                 'family_PGL3': fmpq(1, 6), 'family_PGL2': fmpq(1, 2), 'family_PGL2_simple': fmpq(1, 2)}[key]
        assert xi < bound, 'xi must lie strictly inside the admissible range'
        Q = Q_cheb(h); P = P_sinh(lam)
        if key.endswith('_simple'): assert Q.degree() == 1, 'simple zeros need linear Q'
        c, k = c_kappa(P, Q, R, xi)
        e1, e2 = cross_checks(P, Q, R, xi, c)
        res[key] = k
        print(f"[{label}] xi = {xi} (nu = {2*xi} admissible), R = {R}, lam = {lam}")
        print(f"   deg P = {P.degree()}, deg Q = {Q.degree()}, Q(1) = {float(Q(fmpq(1))):.6e}")
        print(f"   c(P,Q,R) in {c.str(20)}")
        print(f"   kappa    in {k.str(20)}   lower endpoint {k.lower().str(15, radius=False)}")
        print(f"   cross-checks (rel. diff to exact c): AT21 inf_P closed form {e1:.2e} (<= 0 expected: inf over P <= our polynomial P); 2508.11108 (kappa: general) via quad {e2:.2e}")
    print("\n== C. Exact coefficients ==")
    for key in ('AT21_uncond', 'RSC', 'Bernard_uncond', 'family_PGL3', 'family_PGL2', 'family_PGL2_simple'):
        d = prm[key]
        print(f"{key}: xi={d['xi']} R={d['R']} lam={d['lam']}  P(x) = p(x)/p(1), p = sum_(j odd<=25) lam^j x^j / j!")
        print("   Q(y) = 1 + sum_n h_n (T_{2n-1}(1-2y) - 1), S(t) = Q(1/2 + t/(2R)); h_n (exact decimals):")
        print("   " + " ".join(d['h']))
    print("\n== Verdict ==")
    k1, k2, k3 = res['AT21_uncond'], res['RSC'], res['Bernard_uncond']
    print("kappa_f > 0.07005 (uncond.):", bool(k1 > arb(F('0.07005'))), "| > AT21 0.0696:", bool(k1 > arb(F('0.0696'))))
    print("kappa_f > 0.08971 (RSC):", bool(k2 > arb(F('0.08971'))), "| > AT21 0.0896:", bool(k2 > arb(F('0.0896'))))
    print("Bernard length: kappa > 0.06642:", bool(k3 > arb(F('0.06642'))), "| > 2508.11108 0.0632:", bool(k3 > arb(F('0.0632'))))
    print("2607.00282 families: PGL3 kappa(1/6-) > 0.1199:", bool(res['family_PGL3'] > arb(F('0.1199'))), "(stated 1/9);",
          "PGL2 kappa(1/2-) > 0.3658:", bool(res['family_PGL2'] > arb(F('0.3658'))), "(stated 1/3);",
          "PGL2 simple > 0.3562:", bool(res['family_PGL2_simple'] > arb(F('0.3562'))), "(stated 1/3)")
if __name__ == '__main__':
    run()
