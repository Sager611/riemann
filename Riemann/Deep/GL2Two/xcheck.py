"""Independent check of the jet algebra in tp.py: evaluate BCY's c12 (eq:c12) and c2 (eq:c2) literally,
integrand as written in the TeX with complex x, y, and take d^2/dxdy resp. d^4/dx^2dy^2 at 0 by Cauchy's formula
on circles |x| = |y| = rho (trapezoid rule, exact for entire integrands up to aliasing).  Same quadrature nodes as tp.py."""
import numpy as np, sys
from tp import *

def F12(x, y, P1, Q, R, th1, th2, P2dd, n=40):
    p, wp = gl01(n); q, wq = gl01(n); u, wu = gl01(n)
    PP, QQ, UU = np.meshgrid(p, q, u, indexing='ij')
    W = (wp[:, None, None] * wq[None, :, None] * wu[None, None, :]) * (1 - PP)
    a = PP; b = (1 - PP) * QQ
    integrand = (UU ** 2 * (1 - UU) * np.exp(R * (th1 * (y - x) + UU * th2 * (a - b)))
                 * Q(-x * th1 + a * UU * th2) * Q(1 + y * th1 - b * UU * th2)
                 * P1(x + y + 1 - (1 - UU) * th2 / th1) * P2dd((1 - a - b) * UU))
    return np.sum(W * integrand)

def F2(x, y, Q, R, th2, P2dd, nt=60, nr=20, nu=20):
    t, wt = gl01(nt); r, wr = gl01(nr); u, wu = gl01(nu)
    tot = 0
    for ir in range(nr):
        rr = r[ir]
        T, U, V = np.meshgrid(t, u, u, indexing='ij')
        S = x + y - V * (y + rr) - U * (x + rr)
        integrand = ((1 - rr) ** 4 * (1 / th2 + S) * np.exp(-th2 * R * S)
                     * Q(th2 * (-y + U * (x + rr)) + T * (1 + th2 * S))
                     * np.exp(2 * R * T * (1 + th2 * S))
                     * Q(th2 * (-x + V * (y + rr)) + T * (1 + th2 * S))
                     * (x + rr) * (y + rr) * P2dd((1 - U) * (x + rr)) * P2dd((1 - V) * (y + rr)))
        tot = tot + wr[ir] * np.sum(wt[:, None, None] * wu[None, :, None] * wu[None, None, :] * integrand)
    return tot

def cauchy(F, i, j, rho=0.3, N=16):
    """i! j! [x^i y^j] F."""
    from math import factorial
    ph = 2 * np.pi * np.arange(N) / N
    acc = 0
    for a in ph:
        for b in ph:
            x = rho * np.exp(1j * a); y = rho * np.exp(1j * b)
            acc += F(x, y) * x ** (-i) * y ** (-j)
    return (acc / N ** 2).real * factorial(i) * factorial(j)

class CQ:   # complex-capable polynomial in (1-2y) via Horner (low degree)
    def __init__(self, az): self.az = np.asarray(az, float)
    def __call__(self, y):
        z = 1 - 2 * y; s = 0
        for c in self.az[::-1]: s = s * z + c
        return s
class CP:
    def __init__(self, a): self.a = np.asarray(a, float)
    def __call__(self, x):
        s = 0
        for c in self.a[::-1]: s = s * x + c
        return s

cases = [("BCY zeta", [0.492, 0.604, 0, -0.08, 0, -0.06, 0, 0.046], [0, .842706, .00845721, .093117, .118788, -.0630687],
          [0, 0, 0, .0245412, -.00635566, .00603128], 1.28, 4 / 7, 1 / 2),
         ("KRZ uncond (R_eff=6.42)", [.499386, 1.58992, 0, -2.99061, 0, 3.01825, 0, -1.11694], [0, .93271, .147723, -.35572, .444208, -.168921],
          [0, 0, 0, 1.0, -2.0, 1.5], 6.42, 5 / 54, 25 / 298)]
for name, Qz, P1c, P2c, R, th1, th2 in cases:
    Qz = np.asarray(Qz) / np.sum(Qz)
    Qc = CQ(Qz); P1c_ = CP(P1c)
    P2dd = CP(np.polynomial.polynomial.polyder(np.asarray(P2c, float), 2))
    # literal evaluation
    v12 = 4 * (th2 / th1) ** 2 * np.exp(R) * cauchy(lambda x, y: F12(x, y, P1c_, Qc, R, th1, th2, P2dd), 1, 1)
    v2 = 2 / 3 * cauchy(lambda x, y: F2(x, y, Qc, R, th2, P2dd), 2, 2)
    # jet evaluation (tp.py)
    K = 6; B = GBasis(K); g = gvec_from_P2(P2c, K)
    Q = QCheb.from_z_mono(Qz); P1 = PMono(P1c)
    j12 = g @ c12_vec(P1, Q, R, th1, th2, B); M, _ = c2_mat(Q, R, th2, B, 60, 20, 20); j2 = g @ M @ g
    print(f"{name}: c12 literal {v12:.12e} jets {j12:.12e} rel {abs(v12-j12)/abs(j12):.1e};  c2 literal {v2:.12e} jets {j2:.12e} rel {abs(v2-j2)/abs(j2):.1e}", flush=True)
