#!/usr/bin/env python3
"""Exact analysis of Liu's refined stage map (arXiv 2610.26835, Sec. v:sec:refined).
Run: python3 bflat_exact.py   (sympy; ~10 s, <150 MB)
Notation: delta = 2a-1 zero class, x = g/delta amplitude, b, ell lengths, kappa moment parameter.
Constants (defaults = published): lam = 6 (M2 prime capacity n1+n2+lam*kappa*z<=M),
mu = 2 (Moebius moment r+mu*z<=1), al = 5/6 (sixth-power amplification exponent).
"""
import sympy as sp

d, x, b, l, k, r, t = sp.symbols('delta x b ell kappa r t', real=True)
Q = sp.Rational

def count(lam=6, mu=2, al=Q(5, 6)):
    """R_kappa(delta,x): crossing of short (Moebius vs M2) and long (amplified) counts."""
    g = d*x
    AI = 1 - d*r - 2*g*(1 - r)/mu                      # Moebius moment, z_M=(1-r)/mu
    AS = 1 - 2*d*(t - r) - 2*g*(1 - 2*(t - r))/(lam*k)  # M2 with S_m twice, z_S=(1-2m)/(lam k), m=t-r
    rs = sp.solve(sp.Eq(AI, AS), r)[0]
    Rs = sp.simplify(AI.subs(r, rs))
    Rl = 1 - al + (al - d)*t                            # amplified long Moebius, r=t>=1
    ts = sp.solve(sp.Eq(Rs, Rl), t)[0]
    return sp.simplify(Rl.subs(t, ts)), sp.simplify(ts), sp.simplify(rs.subs(t, ts))

def endpoint(R, lam=6):
    """Liu (eq:v-refined-endpoint) with theta = g/(lam kappa)."""
    h = (1 + b + 3*l)/2; ly = (1 + b - l)/2; c = Q(2)/(lam*k)
    C0 = -Q(1, 4) + b/6 + 5*l/4
    return C0 + d/2 + d*x*(l + c*(h/4 - ly/2)) - Q(3, 4)*h*(1 - R)

R, tstar, rstar = count()
# published formula
c1 = 1/(3*k); D = 3 - (1 + 2*c1)*x; P = (2 - 2*c1*x)*(1 - x); J = (Q(5, 6) - d)*D + d*P
R_liu = 1 - d + (Q(5, 6) - d)*d*P/(2*J)
print("count matches Liu (eq:v-refined-count):", sp.simplify(R - R_liu) == 0)
E = endpoint(R)
E2 = E.subs(x, Q(1, 2))
cr_liu = (3 - 324*k**2)*d**2 + (756*k**2 - 54*k - 5)*d + 40*k - 300*k**2
num_b = sp.factor(sp.together(sp.diff(E2, b)))
print("dE/db at x=1/2 =", num_b)
print("  proportional to Liu cr:", sp.simplify(sp.numer(sp.together(sp.diff(E2, b)))/cr_liu))
lstar = sp.solve(sp.Eq(E2, 0), l)[0]
lstar_liu = (324*d**2*k**2 - 72*d**2*k + 3*d**2 - 684*d*k**2 + 126*d*k - 5*d + 450*k**2 - 60*k) / \
            (324*d**2*k**2 - 144*d**2*k + 15*d**2 - 2340*d*k**2 + 450*d*k - 25*d + 2250*k**2 - 300*k)
crit = sp.numer(sp.together(sp.diff(E2, b)))
# on cr=0 the ell-root is independent of b; compare with Liu there (reduce modulo cr)
diffl = sp.numer(sp.together(lstar - lstar_liu))
_k0 = Q(3, 4); _d0 = sp.nsolve(cr_liu.subs(k, _k0), d, 0.469)
print("  ell*(mine,b=0,1/10) - ell*(Liu) at kappa=3/4, delta_c:",
      [sp.N((lstar - lstar_liu).subs({k: _k0, d: _d0, b: bb}), 6) for bb in (0, Q(1, 10))])
# fixed point: F(kappa)=11/12-ell*/4=(1+kappa)/2  <=> ell* = 5/3-2kappa, on cr=0
lst_c = sp.simplify(lstar_liu)
fp = sp.numer(sp.together(lst_c - (Q(5, 3) - 2*k)))
res = sp.factor(sp.resultant(sp.Poly(fp, d), sp.Poly(cr_liu, d)))
print("resultant(delta) factors:", res)
quart = 16686*k**4 - 27405*k**3 + 9936*k**2 + 1290*k - 275
Bs = sp.symbols('B')
qB = sp.expand(quart.subs(k, 2*Bs - 1))
print("quartic in B:", sp.factor(qB), " (Liu: 133488B^4-376596B^3+384534B^2-167541B+26231)")
# exact rational bracket of kappa_flat and B_flat
lo, hi = Q(7497732940799, 10**13), Q(7497732940800, 10**13)
print("sign change of Liu quartic on [lo,hi]:", sp.sign(quart.subs(k, lo)), sp.sign(quart.subs(k, hi)))
import mpmath; mpmath.mp.dps = 40
kf = sp.nsolve(quart, k, 0.7497733, prec=40)
Bf = (1 + kf)/2
dc = sp.nsolve(cr_liu.subs(k, kf), d, 0.469, prec=40)
lf = lstar_liu.subs({k: kf, d: dc})
bf = sp.nsolve(sp.diff(E2, d).subs({k: kf, l: lf}).subs(d, dc), b, 0.0865, prec=40)
print(f"kappa_flat = {sp.N(kf, 25)}\nB_flat     = {sp.N(Bf, 25)}\ndelta_c    = {sp.N(dc, 20)}\nell*       = {sp.N(lf, 20)}\nb*         = {sp.N(bf, 20)}")
print("check 11/12-ell*/4 - B_flat =", sp.N(Q(11, 12) - lf/4 - Bf, 5))
vals = {k: kf, d: dc, x: Q(1, 2), b: bf, l: lf}
print("E(delta_c,1/2) =", sp.N(E.subs(vals), 5), " dE/ddelta =", sp.N(sp.diff(E, d).subs(vals), 5),
      " dE/dx =", sp.N(sp.diff(E, x).subs(vals), 6), " dE/db =", sp.N(sp.diff(E, b).subs(vals), 5))
Rv = R.subs(vals); tv = tstar.subs(vals); rv = rstar.subs(vals); mv = tv - rv
print(f"count R = {sp.N(Rv, 12)}  (density hyp. 1-delta = {sp.N(1 - dc, 8)});  t* = {sp.N(tv, 10)}  r* = {sp.N(rv, 10)}  m* = {sp.N(mv, 10)}")
