#!/usr/bin/env python3
"""q1_margins.py -- explore16, transfers T1/T2 (N6, N2): q-deformed Jensen polynomials of Xi and of Davenport-Heilbronn.

Families, for the GORZ sequence gamma(n) (8 xi(1/2+z) = sum gamma(n) z^{2n}/n!) and J^{d,n}(X) = sum_j C(d,j) gamma(n+j) X^j:
  L_s  (Laguerre / theta family):  sum_j C(d,j) e^{s j^2} gamma(n+j) X^j.   s <= 0: q^{j^2}-damping (q = e^{s}), a multiplier
       sequence; in v = log(-X) it is the backward heat flow on the cylinder (de Bruijn-Newman flow, periodic version;
       Piotrowski-Sorensen-White 1802.02641: zero-sector reducing, cos(angle') = cos(angle)/q).  s > 0: anti-damping.
       margin s*(d,n) = sup{s: L_s hyperbolic}; RH <=> s*(d,n) >= 0 for all d, n (hyperbolic set in s is (-inf, s*]).
  RS_q (Rogers-Szego family):      sum_j [d j]_q gamma(n+j) X^j; gamma == 1 gives the Rogers-Szego polynomial (zeros on |X|=1).
       Prediction: [d j]_q = C(d,j) e^{-(eps/2) j(d-j) + O(eps^2)}, q = e^{-eps}: anti-damping e^{+eps j^2/2}, so RS_q stays
       real-rooted only for eps <~ 2 s*(d,n).
  SS_q (Schur-Szego family):       sum_j [d j]_q q^{j(j-1)/2} gamma(n+j) X^j = J^{d,n} composed with (-X;q)_d; real-rooted for
       every q in [0,1] whenever J^{d,n} is (Schur-Szego composition) -- propagation from q = 1 (RH end) down to q = 0 (linear).
GORZ: gamma(n+j)/gamma(n) = exp(A j - delta^2 j^2 + ...), delta(n)^2 = -(1/2) second difference of log gamma.
Run: uv run --no-project --with python-flint python3 -u q1_margins.py   (single thread; a few minutes)
"""
import sys, time, math
from math import comb
from flint import arb, acb, arb_series, arb_poly, acb_poly, ctx

def gam_xi(N):
    s = arb_series([arb(1) / 2, 1], prec=2 * N + 2)
    xi = (s * (s - 1) / 2) * (-(s / 2) * arb.pi().log()).exp() * arb_series.gamma(s / 2) * arb_series.zeta(s)
    c = xi.coeffs(); out = []; fact = arb(1)
    for n in range(N + 1):
        if n: fact *= n
        out.append(8 * fact * c[2 * n])
    return out

def gam_dh(N):
    s = arb_series([arb(1) / 2, 1], prec=2 * N + 2)
    r5 = arb(5).sqrt(); kap = ((10 - 2 * r5).sqrt() - 2) / (r5 - 1)
    Z = lambda a: arb_series.zeta(s, a)
    f = (-s * arb(5).log()).exp() * (Z(arb(1) / 5) - Z(arb(4) / 5) + kap * (Z(arb(2) / 5) - Z(arb(3) / 5)))
    xi = (((s + 1) / 2) * (arb(5) / arb.pi()).log()).exp() * arb_series.gamma((s + 1) / 2) * f
    c = xi.coeffs(); out = []; fact = arb(1)
    for n in range(N + 1):
        if n: fact *= n
        out.append(fact * c[2 * n])
    return out

def real_count(c):
    """certified lower bound for the number of real zeros of a polynomial with positive coefficients (all real
    zeros negative): sign changes of the ball polynomial at exact points between approximate root real parts"""
    P = arb_poly(c)
    try:
        rr = acb_poly([acb(a.mid()) for a in c]).roots(tol=arb(2) ** (-(ctx.prec // 4)), maxprec=4 * ctx.prec)
    except ValueError:
        return -1                        # near-double root not resolved: counted as not certified
    xs = sorted([arb(r.real.mid()) for r in rr], key=lambda a: float(a))
    pts = [xs[0] * 2 - 1] + [arb(((xs[i] + xs[i + 1]) / 2).mid()) for i in range(len(xs) - 1)] + [arb(0)]
    sg = []
    for p in pts:
        v = P(p); sg.append(1 if v > 0 else (-1 if v < 0 else 0))
    return sum(1 for i in range(len(sg) - 1) if sg[i] * sg[i + 1] < 0)

def hyp(c):
    return real_count(c) == len(c) - 1

def lag(g, d, n, s):
    S = arb(s)
    return [comb(d, j) * (S * j * j).exp() * g[n + j] for j in range(d + 1)]

def qbin(d, q):
    # Gaussian binomials [d j]_q as arb, q an arb in (0, 1]
    row = [arb(1)]
    for m in range(1, d + 1):
        new = [arb(1)] * (m + 1)
        for j in range(1, m):
            new[j] = row[j - 1] + q ** j * row[j]     # [m j] = [m-1 j-1] + q^j [m-1 j]
        row = new
    return row

def rs(g, d, n, eps):
    q = (-arb(eps)).exp(); b = qbin(d, q)
    return [b[j] * g[n + j] for j in range(d + 1)]

def ss(g, d, n, q):
    Q = arb(q); b = qbin(d, Q)
    return [b[j] * Q ** (j * (j - 1) // 2) * g[n + j] for j in range(d + 1)]

def margin(g, d, n, iters=26):
    lo, hi = 0.0, 0.5
    if not hyp(lag(g, d, n, 0.0)):
        return None                      # J^{d,n} itself not hyperbolic (DH beyond its reach)
    while hyp(lag(g, d, n, hi)):
        lo, hi = hi, hi * 2
        if hi > 64: return float('inf')
    for _ in range(iters):
        mid = (lo + hi) / 2
        if hyp(lag(g, d, n, mid)): lo = mid
        else: hi = mid
    return lo

def eps_rs(g, d, n, smax, K=60):
    """first eps on a grid in (0, 6 smax] where RS_{e^-eps} stops being certified hyperbolic; also whether it recovers"""
    first = None; back = False
    for k in range(1, K + 1):
        e = 6 * smax * k / K
        h = hyp(rs(g, d, n, e))
        if not h and first is None: first = e
        if h and first is not None: back = True
    return first, back

def delta2(g, n):
    L = lambda m: float((g[m].log()).mid())
    return -(L(n + 1) - 2 * L(n) + L(n - 1)) / 2

def main():
    t0 = time.time()
    N = 90
    ctx.prec = 4000; ctx.cap = 2 * N + 2
    gx = gam_xi(N); gd = gam_dh(N)
    bad = [(k, i) for k, g in (('Xi', gx), ('DH', gd)) for i, x in enumerate(g) if not x > 0]
    assert not bad, bad[:5]
    ctx.prec = 1200
    print(f'gamma_Xi(0..{N}), gamma_DH(0..{N}) positive ({time.time() - t0:.0f} s)')
    for name, g in (('Xi', gx), ('DH', gd)):
        print(f'\n== {name}: margins s*(d,n) (anti-damping tolerance), GORZ delta(n+d/2)^2, RS threshold eps_RS vs 2 s*')
        print(' d   n     s*(d,n)      delta^2(c)   s*/delta^2   2d*s*     eps_RS     eps_RS/(2s*)  recovers')
        for n in (0, 4, 10, 20, 40):
            for d in (2, 3, 4, 6, 8, 12, 16, 24):
                if n + d + 2 > N: continue
                s = margin(g, d, n)
                c = n + d / 2
                dl = delta2(g, max(1, int(round(c))))
                if s is None:
                    print(f'{d:2d} {n:3d}   J^(d,n) not hyperbolic'); continue
                e, back = eps_rs(g, d, n, s)
                er = (e / (2 * s)) if (e is not None and s > 0) else float('nan')
                print(f'{d:2d} {n:3d}  {s:11.4e}  {dl:11.4e}  {s / dl:9.3f}  {2 * d * s:8.3f}  {e if e is not None else float("nan"):10.4e}  {er:8.3f}  {back}')
            print(f'   ({time.time() - t0:.0f} s)')
    # Schur-Szego family: real-rooted for all q (theorem); numeric confirmation
    bad = 0; tot = 0
    for g in (gx, gd):
        for (d, n) in ((8, 0), (16, 0), (24, 0), (12, 10), (24, 20)):
            for k in range(1, 20):
                tot += 1; bad += 0 if hyp(ss(g, d, n, k / 20)) else 1
    print(f'\nSchur-Szego family SS_q, q = 0.05..0.95, 5 (d,n) pairs, Xi and DH: {tot - bad}/{tot} certified hyperbolic')
    # Rogers-Szego-Jensen zeros below threshold: on a circle?
    print('\nRS_q zeros below the threshold (Xi, d = 12, n = 0): spread of log|X| vs angle spread')
    for eps in (0.5, 1.0, 2.0, 4.0):
        c = rs(gx, 12, 0, eps)
        rr = acb_poly([acb(a.mid()) for a in c]).roots(tol=arb(2) ** (-300), maxprec=4800)
        mods = [math.log(abs(complex(float(r.real.mid()), float(r.imag.mid())))) for r in rr]
        angs = [abs(math.atan2(float(r.imag.mid()), -float(r.real.mid()))) for r in rr]
        print(f'  eps = {eps}: log|X| in [{min(mods):.3f}, {max(mods):.3f}], |arg(-X)| in [{min(angs):.3f}, {max(angs):.3f}], real: {real_count(c)}/12')
    print(f'\ndone ({time.time() - t0:.0f} s)')

if __name__ == '__main__':
    main()
