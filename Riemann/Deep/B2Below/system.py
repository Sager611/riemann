# Scalar contact system (Hexagon2/BeyondB2 audited model of PR103) with PR137's improved count window.
# usage: python3 system.py
import numpy as np, math, sys
from mpmath import mp, mpf, sqrt, findroot
mp.dps = 40
S1081 = sqrt(1081)
B2 = (777 - S1081) / 852; D2 = (61 - S1081) / 66; R2 = (281 - S1081) / 396
b2 = -mpf(11) / 71 + 595 * S1081 / 76751
G0 = mpf(552429) / 10**7; G1 = mpf(7227847) / 10**8
WIN = (189 / 512, 15 / 32, 15 / 32, 0.5)   # delta range, x range of the improved window

def Rd_native(d, x):
    md = (5 / 6 - d * x) / (5 / 2 - d - 2 * d * x)
    return 1 - d * (x + 2 * (1 - x) * md)

def R_line(d, x, g0, g1):
    # ImprovedWindowQ.lean: k=1-x(1+2g1), V=5/6-d, J=3V+2dk, U=1-x(1+6g0-4g1), N=(1-d)J+V d U ; value N/J
    k = 1 - x * (1 + 2 * g1); V = 5 / 6 - d; J = 3 * V + 2 * d * k; U = 1 - x * (1 + 6 * g0 - 4 * g1)
    N = (1 - d) * J + V * d * U
    return N / J, 2 * (J + d * U) / (3 * J)   # value, crossing

def R_curve(d, x, G, lo=0.70, hi=0.80):
    # general gain function G(w): balance 1-d*F(c) = 1-d+V(3c/2-1), F(c)=c+2x((1-c)/2+G(c)); bisection in c
    V = 5 / 6 - d
    f = lambda c: (1 - d * (c + 2 * x * ((1 - c) / 2 + G(c)))) - (1 - d + V * (1.5 * c - 1))
    a, b = lo, hi
    fa = f(a)
    for _ in range(80):
        m = (a + b) / 2; fm = f(m)
        if (fm > 0) == (fa > 0): a, fa = m, fm
        else: b = m
    c = (a + b) / 2
    return 1 - d + V * (1.5 * c - 1), c

def margin(B, b, d, x, R):
    l = 11 / 3 - 4 * B; h = (1 + b + 3 * l) / 2; T = (1 + b - l) - h
    return 0.25 - b / 6 - 5 * l / 4 - d / 2 - d * x * l + 3 * h * (1 - R) / 4 + d * x * T / 4

def Broot(b, d, x, R):
    # margin affine in B: solve margin=0
    m0 = margin(0.0, b, d, x, R); m1 = margin(1.0, b, d, x, R)
    return -m0 / (m1 - m0)

if __name__ == '__main__':
    print("B2 =", mp.nstr(B2, 20), " delta2 =", mp.nstr(D2, 15), " R2 =", mp.nstr(R2, 15), " b2 =", mp.nstr(b2, 15))
    d2 = float(D2)
    print("native Rd(d2,1/2) =", repr(Rd_native(d2, 0.5)))
    v, c = R_line(d2, 0.5, float(G0), float(G1)); print("improved line value at contact =", repr(v), "crossing", repr(c), "dR", v - Rd_native(d2, 0.5))
    v2, c2 = R_curve(d2, 0.5, lambda w: float(G0) - float(G1) * w); print("curve solver check:", repr(v2), repr(c2))
