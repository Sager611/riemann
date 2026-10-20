#!/usr/bin/env python3
"""Refined comparison system of the quasi-RH cubic-theta argument (numerical check).

Sources: OpenAI 30 Sep 2026 (7/8); Liu arXiv 2610.12234; Xu-Sang arXiv 2610.15291 (eq. 27-28);
Liu arXiv 2610.26835 (Sec. 14, eq. 14.3-14.8; B_flat root of 133488B^4-376596B^3+384534B^2-167541B+26231).
Exponents are in base Z; rows u have norm U=Z^d, zero class delta=2a-1, prime amplitude x=q/delta.
Run: python3 refined_system.py   (numpy only, ~20 s, <100 MB)
"""
import math
import numpy as np

def R_count(d, x, k):            # row-count exponent R_kappa(delta,x)  (Xu-Sang (21), Liu-NF (14.3))
    c = 1/(3*k); D = 3-(1+2*c)*x; P = (2-2*c*x)*(1-x); J = (5/6-d)*D+d*P
    return 1-d+(5/6-d)*d*P/(2*J)

def E_end(d, x, lx, ly, l, k, B, numerator="av"):
    """Endpoint exponent at the critical row scale h=1-lx+l, relative to the direct bound Z^{C(B)}."""
    h = 1-lx+l; a = (1+d)/2; q = d*x; R = R_count(d, x, k); c = 1/(3*k)
    Eav = a-B-ly/2-h/6-l/2+l*q+h*(3*R+1)/4+(q*c/2)*(h/2-ly)      # averaged numerator, theta=q/(6k)
    if numerator == "av":
        return Eav
    Ept = a*(1-ly)-B-h/6-l/2+l*q+h*(R+d/2)                         # pointwise numerator U^{delta/2}
    return np.minimum(Eav, Ept)

def B_low(lx, ly, l):            # direct (low) bound matched to C(B)=s-5/6+lx/3+l/6  (Liu Prop. 22.2)
    M, b = lx+ly, ly-lx
    k = max(0, (1+3*l-2*M)/4, M+l-1, (3*l-1)/2)
    return 5/6+lx/6-l/6+max(0, b/12, b-ly/2)+k/2

dg = np.linspace(1/8, 5/6, 1201); xg = np.linspace(0, 0.5, 241)
DD, XX = np.meshgrid(dg, xg, indexing="ij")

def max_E(lx, ly, l, k, B, numerator="av"):
    E = E_end(DD, XX, lx, ly, l, k, B, numerator)
    i = np.unravel_index(E.argmax(), E.shape)
    lo, hi = max(1/8, dg[i[0]]-0.002), min(5/6, dg[i[0]]+0.002)     # refine along delta at that x
    f = lambda d: float(E_end(d, xg[i[1]], lx, ly, l, k, B, numerator))
    g = (math.sqrt(5)-1)/2
    for _ in range(60):
        m1, m2 = hi-g*(hi-lo), lo+g*(hi-lo)
        if f(m1) > f(m2): hi = m2
        else: lo = m1
    return max(E.max(), f((lo+hi)/2)), dg[i[0]], xg[i[1]]

def golden_min(f, a, b, tol=1e-9):
    g = (math.sqrt(5)-1)/2; c = b-g*(b-a); d = a+g*(b-a); fc, fd = f(c), f(d)
    while b-a > tol:
        if fc < fd: b, d, fd = d, c, fc; c = b-g*(b-a); fc = f(c)
        else: a, c, fc = c, d, fd; d = a+g*(b-a); fd = f(d)
    return (a+b)/2, min(fc, fd)

def balanced(b, l):
    return (1-b-l)/2, (1+b-l)/2

def best_balanced(k):
    lo, hi = 0.165, 0.169
    for _ in range(34):
        l = (lo+hi)/2
        _, v = golden_min(lambda b: max_E(*balanced(b, l), l, k, 11/12-l/4)[0], 0.05, 0.12, 1e-8)
        lo, hi = (l, hi) if v <= 0 else (lo, l)
    b, _ = golden_min(lambda b: max_E(*balanced(b, lo), lo, k, 11/12-lo/4)[0], 0.05, 0.12, 1e-8)
    return lo, b

if __name__ == "__main__":
    k = 0.75
    l, b = best_balanced(k)
    B = 11/12-l/4
    for _ in range(3):                       # licensing fixed point kappa = 2B-1
        k = 2*B-1
        l, b = best_balanced(k)
        B = 11/12-l/4
    k = 2*B-1
    Bflat = [r.real for r in np.roots([133488, -376596, 384534, -167541, 26231]) if 0.87 < r.real < 0.88][0]
    v, dc, xc = max_E(*balanced(b, l), l, k, B)
    print(f"1) refined system optimum: B={B:.12f}  (published B_flat={Bflat:.12f})  l={l:.6f} b={b:.6f} kappa={k:.6f}")
    print(f"   critical class delta={dc:.4f} x={xc:.2f}; max endpoint exponent there = {v:.2e}")
    print("2) unbalanced total length M=lx+ly at the same B (largest M with B_low<=B), best b:")
    for dl in (-2e-3, -5e-4, 0.0, 5e-4, 2e-3):
        ll = l+dl
        Mlo, Mhi = (1+3*ll)/2, 1.1
        for _ in range(60):                   # B_low is independent of b on this branch
            Mm = (Mlo+Mhi)/2
            Mlo, Mhi = (Mm, Mhi) if B_low((Mm-0.08)/2, (Mm+0.08)/2, ll) <= B else (Mlo, Mm)

        def worst(bv, M=Mlo, lv=ll):
            return max_E((M-bv)/2, (M+bv)/2, lv, k, B)[0]
        bb, vv = golden_min(worst, 0.0, 0.25, 1e-8)
        print(f"   l-l*={dl:+.4f}  M={Mlo:.6f} (1-l={1-ll:.6f})  min_b max E = {vv:+.2e}")
    vmix, _, _ = max_E(*balanced(b, l), l, k, B, numerator="min")
    print(f"3) class-by-class min(averaged, pointwise) numerator: max E = {vmix:.2e} (unchanged)")
    print("4) low-side ceiling of this probe: min over M,l of B_low with l<=1/5 is 13/15 =", 13/15)

