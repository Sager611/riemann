#!/usr/bin/env python3
"""q5audit: check that the paper's closed form (eqn: SRT), TeX l.440-473 of arXiv 2508.11108, is the BVP solution (FULLPROOF C4).
(a) g1, g2 solve S'' + tanh S' - S = 0; (b) W(u) of (spec: Wronsk) l.473 equals g1 g2' - g1' g2; (c) the paper's constants (l.~750-760);
(d) C1(R) of (const: C1R) l.468 gives S_R(R) = 0 and S_R'(0) = C1 f'(0) + g0'(0) + g1(0)v1(0) + g2(0)v2(0) (finderSP) equals the
independent BVP value (indep_bvp.log). mpmath only. Run: uv run --no-project --with mpmath python3 indep_closed.py
"""
from mpmath import mp, mpf, sqrt, exp, tanh, hyp2f1, quad, diff, nstr
mp.dps = 30
phi = (1 + sqrt(5)) / 2; ip = 1 / phi
Fp = lambda t: hyp2f1(mpf(1)/2, phi, mpf(1)/2 + phi, -exp(2*t))
Fm = lambda t: hyp2f1(mpf(1)/2, -ip, mpf(1)/2 - ip, -exp(2*t))
F1p = lambda u: hyp2f1(mpf(1)/2, 1 + phi, mpf(1)/2 + phi, -exp(2*u))
F1m = lambda u: hyp2f1(mpf(3)/2, -ip, mpf(1)/2 - ip, -exp(2*u))
g1 = lambda t: exp(-ip*t) * Fm(t); g2 = lambda t: exp(phi*t) * Fp(t)
G10, G20 = g1(0), g2(0)
f = lambda t: g1(t) - G10 / G20 * g2(t); g0 = lambda t: g2(t) / (2 * G20)
W = lambda u: exp(u) * (2*phi*Fm(u)*F1p(u) - Fp(u)*F1m(u))
v1 = lambda u: g2(u) / (W(u) * (1 + exp(2*u))); v2 = lambda u: -g1(u) / (W(u) * (1 + exp(2*u)))
print("(a) homogeneous residuals g'' + tanh g' - g:")
for t in [mpf('0.3'), mpf(1), mpf('2.5'), mpf(6)]:
    for name, g in [("g1", g1), ("g2", g2)]:
        res = diff(g, t, 2) + tanh(t) * diff(g, t, 1) - g(t)
        print(f"   {name}({nstr(t,3)}): residual {nstr(res, 3)} (|g| = {nstr(abs(g(t)), 5)})")
print("(b) Wronskian: W(u) vs g1 g2' - g1' g2:")
for u in [mpf(0), mpf('0.7'), mpf(3)]:
    wr = g1(u) * diff(g2, u) - diff(g1, u) * g2(u)
    print(f"   u={nstr(u,3)}: W = {nstr(W(u), 15)}, g1g2'-g1'g2 = {nstr(wr, 15)}, cosh-law W*cosh(u) = {nstr(W(u)*mp.cosh(u), 12)}")
fp0, g0p0 = diff(f, 0), diff(g0, 0)
print(f"(c) f'(0)={nstr(fp0,8)} [paper -1.47277]  g0'(0)={nstr(g0p0,8)} [0.602775]  g1(0)={nstr(G10,8)} [-1.07479]  g2(0)={nstr(G20,8)} [0.759136]"
      f"  v1(0)={nstr(v1(0),8)} [0.339496]  v2(0)={nstr(v2(0),8)} [0.480664]")
bvp = {mpf('0.5'): '-0.4415648493', mpf('0.25'): '-0.3919336928', mpf('0.1'): '-0.3900533226'}
for th, ref in bvp.items():
    R = sqrt(mpf(3)/5) / th
    w1 = quad(v1, [0, R/2, R]); w2 = quad(v2, [0, R/2, R])
    C1 = (-2*Fm(R)*Fp(0)*w1 - exp(sqrt(5)*R)*Fp(R)*(1 + 2*Fp(0)*w2)) / (2*Fm(R)*Fp(0) - 2*exp(sqrt(5)*R)*Fm(0)*Fp(R))
    SR = C1*f(R) + g0(R) + g1(R)*w1 + g2(R)*w2
    Sp0 = C1*fp0 + g0p0 + G10*v1(0) + G20*v2(0)
    print(f"(d) theta={nstr(th,3)} R={nstr(R,10)}: C1={nstr(C1,12)}  S_R(R)={nstr(SR,3)}  S_R'(0) closed={nstr(Sp0,12)}  BVP={ref}")
