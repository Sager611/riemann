"""explore53: the amplifier for the certificate. Float: top eigenvector of (C0, R) for the Arb-certified rational J = 4 form
(logs/cert4_k3.34.json) at kappa, degree d, twist v optimised; converted to GapBar's convention
f(U) = sum_n A_n P_n(2U/theta - 1) (optmid.py: a = R^{-1} y, A = a * norm), A_0 = 1, 16 digits.
usage: amp.py d kappa D_cert out.json [form.json]"""
import sys, os, json
from math import comb, factorial
from fractions import Fraction as Fr
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import moments as Mo
import gcore as G
import joint as Jt
from scipy.optimize import minimize_scalar


def forms_with_basis(r, th, d, kap, K, v0, D=40):
    Nb = d + 2 * r + 6
    H = [h.real for h in Mo.arrays_K(kap * np.pi, r, K, D)]
    N0 = Nb + D + 2
    B = G.Basis(th, N0, nq=N0 + 2 * r * r + 12 + K)
    q0 = r * r - 1
    Wc = [(B.V * (B.w * B.X ** q0 * (B.X - v0) ** p)[:, None]).T @ B.V for p in range(K + 1)]
    Ipow = [np.eye(N0)]
    for k in range(D):
        Ipow.append(B.I @ Ipow[-1])
    Ip = [P[:, :Nb] for P in Ipow]
    s0 = 2 * r - 1
    fac = [float(factorial(k)) for k in range(D + s0 + 2)]
    def Gm(p, h):
        out = np.zeros((Nb, Nb)); L = [Wc[p] @ Ip[k] for k in range(D)]
        for i in range(D):
            for j in range(D):
                if h[i, j] != 0:
                    out += fac[i + s0] * fac[j + s0] * h[i, j] * (Ip[i].T @ L[j])
        return 0.5 * (out + out.T)
    Ir = np.linalg.matrix_power(B.I, 2 * r)[:Nb, :d + 1]
    Q, R = np.linalg.qr(Ir)
    C = [(Q.T @ sum(comb(k, i) * Gm(i, H[k - i]) for i in range(k + 1)) @ Q) * kap ** k for k in range(K + 1)]
    return C, R, B


if __name__ == '__main__':
    d = int(sys.argv[1]); kap = float(sys.argv[2]); Dc = int(sys.argv[3]); out = sys.argv[4]
    cert = json.load(open(os.path.join(HERE, sys.argv[5] if len(sys.argv) > 5 else 'logs/cert4_k3.34.json')))
    a = [float(Fr(x)) for x in cert['a']]
    J = len(a) - 1
    r, th, v0 = 2, 0.2499, 1.3617
    C0, R, B = forms_with_basis(r, th, d, kap, 2 * J, v0)
    def h(dv):
        Cs = Jt.shift(C0, dv, kap)
        return -a[0] * Jt.lmax(Cs[0], sum(a[j] * Cs[2 * j] for j in range(1, J + 1)))[0]
    rr = minimize_scalar(lambda x: -h(x), bounds=(-0.005, 0.005), method='bounded', options={'xatol': 1e-10})
    dv = rr.x
    v = round(v0 + dv, 6)
    Cs = Jt.shift(C0, v - v0, kap)
    lam, y = Jt.lmax(Cs[0], sum(a[j] * Cs[2 * j] for j in range(1, J + 1)))
    co = np.linalg.solve(R, y)
    An = co * B.norm[:d + 1]; An = An / An[0]
    A = [('%.15f' % x) if abs(x) >= 1 else ('%.16f' % x) for x in An]
    # float check of the rounded amplifier
    Ar = np.array([float(x) for x in A]) / B.norm[:d + 1]
    yr = R @ Ar
    Hf = -sum(a[j] * (yr @ Cs[2 * j] @ yr) for j in range(J + 1)) / (yr @ Cs[0] @ yr)
    print('d=%d kappa=%.4f: float h_J = %.8f at v = %.6f; rounded amplifier: H/c0 (float) = %.8f' % (d, kap, -a[0] * lam, v, Hf))
    json.dump({'r': r, 'theta': '0.2499', 'd': d, 'kappa': ('%.6f' % kap).rstrip('0'), 'v': '%.6f' % v, 'A': A, 'a': cert['a'], 'D': Dc,
               'out': 'logs/certJ_d%d_k%s_J%d.json' % (d, sys.argv[2], J), 'H_over_c0_float': Hf}, open(os.path.join(HERE, out), 'w'), indent=1)
