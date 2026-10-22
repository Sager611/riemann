#!/usr/bin/env python3
"""Part 3: side constraints at B_flat, prime supply over the rectangle, and three improvement attempts.
Run: python3 bflat_checks.py  (numpy+mpmath; ~20 s, <150 MB)"""
import numpy as np, mpmath as mm
mm.mp.dps = 30
kf = mm.mpf('0.7497732940799282618892280'); Bf = (1 + kf)/2
dc = mm.mpf('0.46934485869814306749'); lf = mm.mpf('0.16712007850681014289'); bf = mm.mpf('0.086534959255841828205')
h = (1 + bf + 3*lf)/2; ly = (1 + bf - lf)/2; lx = (1 - bf - lf)/2; L = lf/h
def cnt(dl, xx, kk, al=5/6):
    c = 2/(3*kk); D = 3 - (1 + c)*xx; P = (2 - c*xx)*(1 - xx); J = (al - dl)*D + dl*P
    s = dl*P/(2*J); t = 1 + s; R = 1 - dl + (al - dl)*s
    r = ((2 - c*xx)*t - xx + c*xx/2)/D
    return R, t, r, t - r
R, t, r, m = cnt(dc, 0.5, kf)
g = dc/2
print("== side constraints at the B_flat tangent point (value, bound, slack) ==")
rows = [
 ("ell >= 1/6 (Liu adm. 1)", lf, mm.mpf(1)/6, lf - mm.mpf(1)/6),
 ("ell <= 1/5 (13/15 floor inactive)", lf, mm.mpf(1)/5, mm.mpf(1)/5 - lf),
 ("B(ell) > 437/500", Bf, mm.mpf(437)/500, Bf - mm.mpf(437)/500),
 ("l_x > ell", lx, lf, lx - lf),
 ("l_y - ell > 11b/6", ly - lf, 11*bf/6, ly - lf - 11*bf/6),
 ("6 kappa ell > 3/4", 6*kf*lf, mm.mpf(3)/4, 6*kf*lf - mm.mpf(3)/4),
 ("ell/h < 1/(6 kappa)", L, 1/(6*kf), 1/(6*kf) - L),
 ("ell/h > 21/100", L, mm.mpf(21)/100, L - mm.mpf(21)/100),
 ("h/2 < l_y < 1/2", ly, h/2, min(ly - h/2, mm.mpf(1)/2 - ly)),
 ("count: 1 <= t* <= 3/2", t, 1, min(t - 1, 1.5 - t)),
 ("count: t*-1/2 <= r* <= 1", r, t - 0.5, min(r - (t - 0.5), 1 - r)),
 ("count: 0 <= m* <= 1/2", m, 0.5, 0.5 - m),
 ("Moebius 2r+8z <= 3 (z=(1-r)/2)", 2*r + 4*(1 - r), 3, 3 - 2*r - 4*(1 - r)),
 ("supply: z_M=(1-r*)/2 <= ell/h", (1 - r)/2, L, L - (1 - r)/2),
 ("supply: z_S=(1-2m*)/(6k) <= ell/h", (1 - 2*m)/(6*kf), L, L - (1 - 2*m)/(6*kf)),
 ("selected plain length z_S <= c/12=1/(18k)", (1 - 2*m)/(6*kf), 1/(18*kf), 1/(18*kf) - (1 - 2*m)/(6*kf)),
 ("Xu-Sang 1/(18k) <= 1/13 (k>=13/18)", 1/(18*kf), mm.mpf(1)/13, mm.mpf(1)/13 - 1/(18*kf)),
 ("OpenAI plain capacity 1/(18k) <= 2/27 (k>=3/4)", 1/(18*kf), mm.mpf(2)/27, mm.mpf(2)/27 - 1/(18*kf)),
 ("numerator gate R <= (6k-2g) ell/h", R, (6*kf - 2*g)*L, (6*kf - 2*g)*L - R),
 ("theta=g/(6k) <= min(a-1/2,1/10)", g/(6*kf), min(dc/2, mm.mpf(1)/10), min(dc/2, mm.mpf(1)/10) - g/(6*kf)),
 ("M2 centred range used: 2m*+z_S > 5/6", 2*m + (1 - 2*m)/(6*kf), mm.mpf(5)/6, 2*m + (1 - 2*m)/(6*kf) - mm.mpf(5)/6),
 ("kappa_flat >= 13/18 (M2 range)", kf, mm.mpf(13)/18, kf - mm.mpf(13)/18),
 ("kappa_flat >= 11/15 (13/15 floor)", kf, mm.mpf(11)/15, kf - mm.mpf(11)/15),
]
C0v = -mm.mpf(1)/4 + bf/6 + 5*lf/4; Asm = mm.mpf(1)/2 + (lf/2)*(mm.mpf(1)/2 + 2*lf/h) - h/2
rows += [
 ("adm.4a max(C0, C0+A_sm/8) < -1/250", max(C0v, C0v + Asm/8), -mm.mpf(1)/250, -mm.mpf(1)/250 - max(C0v, C0v + Asm/8)),
 ("adm.4b C0+(1/2+3ell/2)/50 < -1/100", C0v + (mm.mpf(1)/2 + 3*lf/2)/50, -mm.mpf(1)/100, -mm.mpf(1)/100 - C0v - (mm.mpf(1)/2 + 3*lf/2)/50),
 ("adm.4c 13h/75-l_y/2+1/50 < -1/20", 13*h/75 - ly/2 + mm.mpf(1)/50, -mm.mpf(1)/20, -mm.mpf(1)/20 - (13*h/75 - ly/2 + mm.mpf(1)/50)),
]
for nm, v, bd, s in rows:
    print(f"  {nm:48s} {mm.nstr(v, 10):>14s} vs {mm.nstr(bd, 10):>14s}   slack {mm.nstr(s, 4):>10s}")
print(f"  derived: h={mm.nstr(h,10)} l_x={mm.nstr(lx,10)} l_y={mm.nstr(ly,10)} R={mm.nstr(R,12)} K=(3R+1)/4={mm.nstr((3*R+1)/4,10)} theta={mm.nstr(g/(6*kf),10)}")

# true small-delta margin with the exact count instead of the coarse bound R <= 1-2delta/3
Ds, Xs = np.meshgrid(np.linspace(1e-4, 1/8, 400), np.linspace(0, 0.5, 51), indexing='ij')
_h = float(h); _l = float(lf); _C0 = float(C0v)
Es = _C0 + Ds/2 + _l*Ds*Xs*(0.5 + 2*_l/_h) - 0.75*_h*(1 - cnt(Ds, Xs, float(kf))[0])
print(f"  small-delta endpoint (delta<1/8) with the exact count: max = {Es.max():+.5f} (coarse bound used by Liu: {float(max(C0v, C0v + Asm/8)):+.5f})")
# ---- prime supply over the whole rectangle at the B_flat point
kk = float(kf); bb = float(bf); ll = float(lf)
hh = (1 + bb + 3*ll)/2; LL = ll/hh
Dg, Xg = np.meshgrid(np.linspace(1/8, 5/6, 711), np.linspace(0, 0.5, 101), indexing='ij')
Rg, tg, rg, mg = cnt(Dg, Xg, kk)
zM = (1 - rg)/2; zS = (1 - 2*mg)/(6*kk)
def Eav(D_, X_, b_=bb, l_=ll, k_=kk):
    h_ = (1 + b_ + 3*l_)/2; ly_ = (1 + b_ - l_)/2; c_ = 1/(3*k_)
    R_ = cnt(D_, X_, k_)[0]
    return -0.25 + b_/6 + 5*l_/4 + D_/2 + D_*X_*(l_ + c_*(h_/4 - ly_/2)) - 0.75*h_*(1 - R_)
Eg = Eav(Dg, Xg)
need = np.maximum(zM*(Xg > 0), zS*(Xg > 0))
i = np.unravel_index(np.argmax(need), need.shape)
print(f"\n== prime supply at the endpoint d=h: supply ell/h={LL:.5f}; max needed over classes = {need[i]:.5f} at delta={Dg[i]:.3f}, x={Xg[i]:.2f} (E there {Eg[i]:+.4f})")
print(f"   range checks over rectangle: min t={tg.min():.4f} max t={tg.max():.4f}; min(r-(t-1/2))={np.min(rg-(tg-0.5)):.4f}; max r={rg.max():.4f}; max m={mg.max():.4f}; max E={Eg.max():+.2e}")

# ---- attempt 1: class-by-class choice of averaged / pointwise numerator at kappa_flat
def Ept(D_, X_, b_, l_, k_=kk):
    h_ = (1 + b_ + 3*l_)/2; R_ = cnt(D_, X_, k_)[0]
    return -0.25 + b_/6 + 5*l_/4 + (0.5 + l_)*D_ + l_*D_*X_ - h_*(1 - R_)
Df, Xf = np.meshgrid(np.linspace(0.40, 0.55, 3001), np.array([0.46, 0.48, 0.49, 0.5]), indexing='ij')
def V(l_, b_):      # worst class of the mixed system: full rectangle (coarse) + critical window (fine)
    return max(np.max(np.minimum(Eav(Dg, Xg, b_, l_), Ept(Dg, Xg, b_, l_))),
               np.max(np.minimum(Eav(Df, Xf, b_, l_), Ept(Df, Xf, b_, l_))))
def minb(l_):
    bs = np.linspace(0.0, 0.25, 251); vs = [V(l_, b_) for b_ in bs]; j = int(np.argmin(vs))
    lo, hi = bs[max(j-1, 0)], bs[min(j+1, 250)]; gr = (5**0.5 - 1)/2
    for _ in range(40):
        m1, m2 = hi - gr*(hi - lo), lo + gr*(hi - lo)
        if V(l_, m1) < V(l_, m2): hi = m2
        else: lo = m1
    return V(l_, (lo + hi)/2), (lo + hi)/2
for dl_ in (0.0, 1e-6, 1e-5, 1e-4):
    v_, b_ = minb(ll + dl_)
    print(f"== attempt: min(averaged, pointwise) numerator per class, kappa_flat, ell = ell*+{dl_:.0e}: min_b max_class = {v_:+.3e} (b={b_:.5f})")
print(f"   pointwise endpoint at the critical class (tangent point) = {float(Ept(np.array(float(dc)), np.array(0.5), bb, ll)):+.2e}")

# ---- attempt 2: Heath-Brown identity split of the long Moebius witness (t*,m=0): M_a M_a S_b, a<=t/2, b<=1/2
c = 2/(3*kk); x = 0.5; dl = float(dc); tt = float(t)
A_I = lambda a: 1 - dl*(x + (1 - x)*a)
A_S = lambda bb_: 1 - dl*(c*x/2 + (2 - c*x)*bb_)
aa = np.linspace(max(0, (tt - 0.5)/2), tt/2, 20001)
split = np.max(np.minimum(A_I(aa), A_S(tt - 2*aa)))
print(f"== attempt: long witness M_t (t={tt:.4f}) split as M_a M_a S_b: best count exponent {split:.5f} vs amplified e_L(t)={float(R):.5f}")

# ---- attempt 3: amplification with prime factors in the enlarged family (H=U^(t+2z))
zz = np.linspace(0, 0.3, 7)
eL = 1/6 + 5*(tt + 2*zz)/6 - dl*tt - 2*(dl/2)*zz
print("== attempt: amplified count with selected primes, exponent vs z:", " ".join(f"{z_:.2f}:{e_:.4f}" for z_, e_ in zip(zz, eL)))

# ---- exact rational re-evaluation of every side constraint at a rational point within 1e-20 of the tangent point
from fractions import Fraction as F
qk, qd, ql, qb = (F(mm.nstr(v, 25)) for v in (kf, dc, lf, bf))
qh = (1 + qb + 3*ql)/2; qly = (1 + qb - ql)/2; qlx = (1 - qb - ql)/2; qL = ql/qh
qc = F(2)/(3*qk); qx = F(1, 2); qD = 3 - (1 + qc)*qx; qP = (2 - qc*qx)*(1 - qx); qJ = (F(5, 6) - qd)*qD + qd*qP
qs = qd*qP/(2*qJ); qt = 1 + qs; qR = 1 - qd + (F(5, 6) - qd)*qs; qr = ((2 - qc*qx)*qt - qx + qc*qx/2)/qD; qm = qt - qr
qg = qd/2; qC0 = -F(1, 4) + qb/6 + 5*ql/4; qA = F(1, 2) + (ql/2)*(F(1, 2) + 2*ql/qh) - qh/2
exact = {
 "ell-1/6": ql - F(1, 6), "1/5-ell": F(1, 5) - ql, "B-437/500": F(11, 12) - ql/4 - F(437, 500),
 "lx-ell": qlx - ql, "ly-ell-11b/6": qly - ql - 11*qb/6, "6k ell-3/4": 6*qk*ql - F(3, 4),
 "1/(6k)-ell/h": 1/(6*qk) - qL, "ell/h-21/100": qL - F(21, 100), "ly-h/2": qly - qh/2, "1/2-ly": F(1, 2) - qly,
 "-max(C0,C0+A/8)-1/250": -max(qC0, qC0 + qA/8) - F(1, 250), "-C0-(1/2+3ell/2)/50-1/100": -qC0 - (F(1, 2) + 3*ql/2)/50 - F(1, 100),
 "-(13h/75-ly/2+1/50)-1/20": -(13*qh/75 - qly/2 + F(1, 50)) - F(1, 20),
 "t-1": qt - 1, "3/2-t": F(3, 2) - qt, "r-(t-1/2)": qr - (qt - F(1, 2)), "1-r": 1 - qr, "1/2-m": F(1, 2) - qm,
 "3-2r-8zM": 3 - 2*qr - 4*(1 - qr), "L-zM": qL - (1 - qr)/2, "L-zS": qL - (1 - 2*qm)/(6*qk),
 "1/(18k)-zS": 1/(18*qk) - (1 - 2*qm)/(6*qk), "1/13-1/(18k)": F(1, 13) - 1/(18*qk), "2/27-1/(18k) (OpenAI, unused)": F(2, 27) - 1/(18*qk),
 "gate (6k-2g)L-R": (6*qk - 2*qg)*qL - qR, "2m+zS-5/6 (centred M2 range)": 2*qm + (1 - 2*qm)/(6*qk) - F(5, 6),
 "k-13/18": qk - F(13, 18), "k-11/15": qk - F(11, 15),
}
neg = [kk_ for kk_, v in exact.items() if v <= 0]
print("\n== exact rational slacks at the rational point (all > 0 except the unused OpenAI 2/27 bound) ==")
print("  min slack:", min((float(v), kk_) for kk_, v in exact.items() if v > 0), "  nonpositive:", neg)
