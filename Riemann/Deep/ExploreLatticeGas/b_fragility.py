#!/usr/bin/env python3
"""b_fragility.py -- explore15, transfer A (N6): why no contraction/local-inequality proof can work.

(1) Fragility.  For J = J^{d,0} of Xi (positive coefficients a_j, zeros r_1 > ... > r_d < 0) let
      eps(d) = min over the local extrema x* between consecutive zeros of |J(x*)| / sum_j a_j |x*|^j.
    A relative perturbation of the coefficients of size eps(d) can merge two zeros into a complex pair, so any
    argument that propagates coefficient inequalities with a fixed relative slack stops at d ~ log(1/slack)/c.
    Same measure for Hermite He_d (= matching polynomial of K_d, real-rooted by Heilmann-Lieb) and for the
    Laguerre-Jensen polynomials of e^y.
(2) Favard.  Heilmann-Lieb works through an exact three-term vertex recursion.  Test whether J^{d+1,0} =
    (A X + B) J^{d,0} - C J^{d-1,0} holds for Xi (it does exactly for e^y: Laguerre).
gamma(n) as in Riemann/Deep/Jensen9/verify.py: 8 xi(1/2+z) = sum gamma(n) z^{2n}/n!  (Arb power series).
Run: uv run --no-project --with python-flint python3 -u b_fragility.py
"""
import time
from math import comb, factorial
from flint import arb, acb, arb_series, acb_poly, arb_poly, ctx

def gammas(N, prec):
    ctx.prec = prec; ctx.cap = 2 * N + 2
    s = arb_series([arb(1) / 2, 1], prec=2 * N + 2)
    xi = (s * (s - 1) / 2) * (-(s / 2) * arb.pi().log()).exp() * arb_series.gamma(s / 2) * arb_series.zeta(s)
    c = xi.coeffs(); out = []; fact = arb(1)
    for n in range(N + 1):
        if n: fact *= n
        out.append(8 * fact * c[2 * n])
    return out

def certified_real_count(coef):
    """lower bound for the number of real zeros: sign changes of P at exact points between the sorted real
    parts of approximate roots (ball arithmetic; a ball containing 0 counts as no sign change)"""
    P = arb_poly([arb(c) for c in coef])
    rs = sorted(float(r.real.mid()) for r in acb_poly([acb(c) for c in coef]).roots())
    xs = sorted(set(arb(r.real.mid()) for r in acb_poly([acb(c) for c in coef]).roots()), key=lambda a: float(a))
    d = len(coef) - 1
    pts = [xs[0] - 1 - abs(xs[0])] + [arb(((xs[i] + xs[i + 1]) / 2).mid()) for i in range(len(xs) - 1)] + [xs[-1] + 1 + abs(xs[-1])]
    sg = []
    for p in pts:
        v = P(p)
        sg.append(1 if v > 0 else (-1 if v < 0 else 0))
    cnt = sum(1 for i in range(len(sg) - 1) if sg[i] * sg[i + 1] < 0)
    return cnt, d

def fragility(coef):
    """coef: list of arb/ints, low degree first; returns (log10 eps, all zeros real?)"""
    P = acb_poly([acb(c) for c in coef])
    cnt, d = certified_real_count(coef)
    if cnt != d:
        return float('nan'), False
    D = P.derivative()
    cs = D.roots()
    A = arb_poly([abs(arb(c)) for c in coef])
    best = None
    for c in cs:
        x = c.real.mid()
        num = abs(arb_poly(coef)(arb(x)))
        den = A(abs(arb(x)))
        v = (num / den).log() / arb(10).log()
        v = float(v.mid())
        best = v if best is None or v < best else best
    return best, True

def main():
    t0 = time.time()
    N = 160
    g = gammas(N, 4400)
    print(f'gamma(0..{N}) computed ({time.time() - t0:.1f} s); gamma(0) = {g[0].str(15)}')
    ctx.prec = 4400
    print('(1) fragility log10 eps(d):  d | Xi J^{d,0} | Xi J^{d,40} | Hermite He_d | Laguerre L_d(-X) [gamma_j = 1/j!]')
    rows = []
    for d in [4, 8, 12, 16, 24, 32, 48, 64, 80, 100, 120]:
        cx = [comb(d, j) * g[j] / g[0] for j in range(d + 1)]
        cx40 = [comb(d, j) * g[40 + j] / g[40] for j in range(d + 1)]
        he = [0] * (d + 1)
        for k in range(d // 2 + 1):
            he[d - 2 * k] = (-1) ** k * factorial(d) // (factorial(k) * 2 ** k * factorial(d - 2 * k))
        la = [arb(comb(d, j)) / factorial(j) for j in range(d + 1)]
        fx, okx = fragility(cx); fx40, ok40 = fragility(cx40); fh, okh = fragility(he); fl, okl = fragility(la)
        rows.append((d, fx, fx40, fh, fl))
        print(f'   {d:4d} | {fx:9.2f} | {fx40:9.2f} | {fh:9.2f} | {fl:9.2f}   (all real: {okx},{ok40},{okh},{okl})  [{time.time() - t0:.0f} s]')
    # slopes
    import math
    d1, d2 = rows[-4][0], rows[-1][0]
    for name, i in [('Xi n=0', 1), ('Xi n=40', 2), ('Hermite', 3), ('Laguerre L_d(-X)', 4)]:
        s = (rows[-1][i] - rows[-4][i]) / (d2 - d1)
        print(f'   slope d log10 eps / dd on [{d1},{d2}]: {name}: {s:.4f}  (eps ~ 10^({s:.3f} d))')
    print('(2) Favard three-term test: J^{d+1} - (A X + B) J^d + C J^{d-1}, A,B,C fitted on X^{d+1}, X^0, X^1;'
          ' max relative residual of the other coefficients')
    for name, seq in [('Xi', [g[j] / g[0] for j in range(42)]), ('(1+X)^d [gamma_j = 1]', [arb(1)] * 42),
                      ('Laguerre L_d(-X) [gamma_j = 1/j!]', [arb(1) / arb(factorial(j)) for j in range(42)])]:
        out = []
        for d in [2, 4, 8, 16, 30]:
            Jm, J0, Jp = ([comb(e, j) * seq[j] for j in range(e + 1)] for e in (d - 1, d, d + 1))
            A = Jp[-1] / J0[-1]
            # X^0: Jp[0] = B J0[0] - C Jm[0];  X^1: Jp[1] = A J0[0] + B J0[1] - C Jm[1]
            a11, a12, b1 = J0[0], -Jm[0], Jp[0]
            a21, a22, b2 = J0[1], -Jm[1], Jp[1] - A * J0[0]
            det = a11 * a22 - a12 * a21
            B = (b1 * a22 - a12 * b2) / det; C = (a11 * b2 - b1 * a21) / det
            res = 0.0
            for k in range(2, d + 1):
                lhs = Jp[k]
                rhs = A * J0[k - 1] + B * (J0[k] if k <= d else 0) - C * (Jm[k] if k <= d - 1 else 0)
                r = abs(float(((lhs - rhs) / lhs).mid()))
                res = max(res, r)
            out.append(f'd={d}: {res:.2e} (C={float(C.mid()):.3g})')
        print(f'   {name}: ' + '; '.join(out))
    print(f'done ({time.time() - t0:.0f} s)')

if __name__ == '__main__':
    main()
