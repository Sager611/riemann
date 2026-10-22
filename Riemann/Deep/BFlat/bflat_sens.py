#!/usr/bin/env python3
"""Part 2: stage map, side constraints and sensitivities of B_flat (Liu arXiv 2610.26835).
Run: python3 bflat_sens.py  (sympy+mpmath+numpy; ~30 s, <200 MB)"""
import sympy as sp, mpmath as mm, numpy as np
mm.mp.dps = 30
d, x, b, l, k, r, t = sp.symbols('delta x b ell kappa r t', real=True)
lam, mu, al, p, xm = sp.symbols('lam mu al p xm', positive=True)
Q = sp.Rational

def count_sym():
    g = d*x
    AI = 1 - d*r - 2*g*(1 - r)/mu
    AS = 1 - 2*d*(t - r) - 2*g*(1 - 2*(t - r))/(lam*k)
    rs = sp.solve(sp.Eq(AI, AS), r)[0]
    Rs = AI.subs(r, rs)
    Rl = 1 - al + (al - d)*t
    ts = sp.solve(sp.Eq(Rs, Rl), t)[0]
    return sp.simplify(Rl.subs(t, ts)), sp.simplify(ts), sp.simplify(rs.subs(t, ts))

R, TS, RS = count_sym()
h = (1 + b + 3*l)/2; ly = (1 + b - l)/2; lx = (1 - b - l)/2
C0 = -Q(1, 4) + b/6 + 5*l/4
g = d*x
Kp = ((p - 1)*R + 1)/p
theta = g/(lam*k)                       # used when its coefficient is negative (checked below)
coef_th = h*(1 - 2/p) - ly
E = C0 + d/2 + l*g + h*(Kp - 1) + theta*coef_th
pub = {lam: 6, mu: 2, al: Q(5, 6), p: 4}
E_pub = E.subs(pub)
c1 = 1/(3*k); Dk = 3 - (1 + 2*c1)*x; Pk = (2 - 2*c1*x)*(1 - x); Jk = (Q(5, 6) - d)*Dk + d*Pk
R_liu = 1 - d + (Q(5, 6) - d)*d*Pk/(2*Jk)
E_liu = C0 + d/2 + d*x*(l + c1*(h/4 - ly/2)) - Q(3, 4)*h*(1 - R_liu)
print("generic E reduces to Liu's endpoint:", sp.simplify(E_pub - E_liu) == 0)

# ---- generic fixed point: unknowns (delta, kappa, b); ell = 5/3 - 2 kappa  (B = (1+kappa)/2 = 11/12 - ell/4)
def solve_fp(params, xcap=Q(1, 2), guess=(0.4693, 0.74977, 0.0865)):
    Ex = E.subs(params).subs(x, xcap).subs(l, Q(5, 3) - 2*k)
    f1 = sp.lambdify((d, k, b), sp.diff(Ex, b), 'mpmath')
    f2 = sp.lambdify((d, k, b), Ex, 'mpmath')
    f3 = sp.lambdify((d, k, b), sp.diff(Ex, d), 'mpmath')
    sol = mm.findroot(lambda u, v, w: [f1(u, v, w), f2(u, v, w), f3(u, v, w)], guess)
    dd, kk, bb = [sol[i] for i in range(3)]
    return dd, kk, bb, (1 + kk)/2

def rect_max(params, dd, kk, bb, xcap=Q(1, 2)):
    """max of E over delta in [1/8,5/6], x in [0,xcap] at the fixed point (tangent check)."""
    lv = Q(5, 3) - 2*kk
    fE = sp.lambdify((d, x), E.subs(params).subs({k: kk, b: bb, l: lv}), 'numpy')
    D_, X_ = np.meshgrid(np.linspace(1/8, 5/6, 1401), np.linspace(0, float(xcap), 101), indexing='ij')
    V = fE(D_, X_)
    i = np.unravel_index(np.nanargmax(V), V.shape)
    return float(V[i]), float(D_[i]), float(X_[i])

def report(name, params, xcap=Q(1, 2), guess=(0.4693, 0.74977, 0.0865)):
    dd, kk, bb, BB = solve_fp(params, xcap, guess)
    vals = dict(params); vals.update({d: dd, k: kk, b: bb, x: xcap})
    tv = float(TS.subs(vals)); rv = float(RS.subs(vals)); mv = tv - rv
    ct = float(coef_th.subs(vals).subs(l, Q(5, 3) - 2*kk))
    mx = rect_max(params, dd, kk, bb, xcap)
    ok = (1 <= tv <= 1.5) and (tv - 0.5 <= rv <= min(tv, 1)) and (0 <= mv <= 0.5) and ct < 0
    print(f"{name:42s} B={mm.nstr(BB, 13):>16s} dB={mm.nstr(BB - BF, 4):>11s} delta_c={mm.nstr(dd, 6)} t={tv:.4f} r={rv:.4f} m={mv:.4f} "
          f"rectmax={mx[0]:+.1e}@({mx[1]:.3f},{mx[2]:.2f}) {'ok' if ok else 'CHECK'}")
    return BB

BF = mm.mpf('0.8748866470399641309446140')
print("\n== fixed point under changed constants (same pipeline; only the named constant changes) ==")
B0 = report("published (lam=6, mu=2, al=5/6, p=4, x<=1/2)", pub)
for v in (Q(11, 2), 5, 4, 3):
    report(f"M2 prime capacity 6kappa -> {v}kappa", {**pub, lam: v})
for v in (Q(9, 5), Q(7, 4), Q(3, 2)):
    report(f"Moebius capacity r+2z -> r+{v}z", {**pub, mu: v})
for v in (Q(4, 5), Q(3, 4), Q(2, 3)):
    report(f"amplification 5/6 -> {v}", {**pub, al: v})
for v in (Q(9, 20), Q(2, 5)):
    report(f"amplitude cap x<=1/2 -> x<={v}", pub, xcap=v)

print("\n== derivatives at the published point (central differences) ==")
for name, sym, base in (("lam", lam, 6), ("mu", mu, 2), ("al", al, Q(5, 6)), ("p (numerator moment order)", p, 4)):
    hstep = Q(1, 10**6)
    Bp = solve_fp({**pub, sym: base + hstep})[3]; Bm = solve_fp({**pub, sym: base - hstep})[3]
    print(f"dB_flat/d{name} = {mm.nstr((Bp - Bm)/(2*hstep), 8)}")
Bp = solve_fp(pub, xcap=Q(1, 2) + Q(1, 10**6))[3]; Bm = solve_fp(pub, xcap=Q(1, 2) - Q(1, 10**6))[3]
print(f"dB_flat/dxcap = {mm.nstr((Bp - Bm)/(2*Q(1, 10**6)), 8)}")

# ---- stage map Phi(B) = F_flat(2B-1)
print("\n== stage map B -> B' = F_flat(2B-1) (published constants) ==")
E2 = E_pub.subs(x, Q(1, 2))
crf = sp.lambdify((d, k), sp.numer(sp.together(sp.diff(E2, b))), 'mpmath')
lsol = sp.solve(sp.Eq(E2, 0), l)[0]
lsf = sp.lambdify((d, k, b), lsol, 'mpmath')
def Fflat(kk):
    dd = mm.findroot(lambda u: crf(u, kk), 0.469)
    return mm.mpf(11)/12 - lsf(dd, kk, 0)/4      # independent of b on cr=0
Bj = mm.mpf(7)/8
for j in range(6):
    Bn = Fflat(2*Bj - 1)
    print(f"  B_{j} = {mm.nstr(Bj, 16)} -> B_{j+1} = {mm.nstr(Bn, 16)}   B_{j+1}-B_flat = {mm.nstr(Bn - BF, 3)}")
    Bj = Bn
eps = mm.mpf('1e-8')
print("  Phi'(B_flat) =", mm.nstr((Fflat(2*(BF + eps) - 1) - Fflat(2*(BF - eps) - 1))/(2*eps), 6),
      "  F_flat(1) =", mm.nstr(Fflat(1), 12), "  F_flat(1/2) [formal] =", mm.nstr(Fflat(mm.mpf(1)/2), 12))
