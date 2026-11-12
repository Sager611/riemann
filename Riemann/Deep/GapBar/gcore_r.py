"""general r (float): weight x^{r^2-1}(x1x2x3x4)^{r-1} (BM Thm smoothed). phi = I^{2r} f;
Q = sum_{i,j} (i+2r-1)!(j+2r-1)! rhobar^(r)_ij <I^i phi, X^{r^2-1+p} I^j phi>, rhobar^(r) = fiber average with weight (t(1-t))^{r-1} per pair."""
import numpy as np
from numpy.polynomial.legendre import leggauss
from math import factorial
import gcore as G


def gbar_r(r1, r2, c, r, nth=28):
    t, w = G.gl01(nth)
    wr = w * (t * (1 - t)) ** (r - 1)
    T1, T2 = np.meshgrid(t, t, indexing='ij')
    W = np.outer(wr, wr)
    sh = r1.shape
    R1a = r1.reshape(-1, 1, 1); R2a = r2.reshape(-1, 1, 1)
    out = [np.zeros(R1a.shape[0], complex) for _ in range(3)]
    ch = 128
    for s in range(0, R1a.shape[0], ch):
        R1 = R1a[s:s + ch]; R2 = R2a[s:s + ch]
        p = G.rho_parts(R1 * T1, R1 * (1 - T1), R2 * T2, R2 * (1 - T2), c)
        for k in range(3):
            out[k][s:s + ch] = (p[k] * W).sum(axis=(1, 2))
    return [o.reshape(sh) for o in out]


def kap_arrays_r(c, r, D=36, R=0.4, K=64):
    ang = np.exp(2j * np.pi * np.arange(K) / K)
    Z1, Z2 = np.meshgrid(R * ang, R * ang, indexing='ij')
    gs = []
    for g in gbar_r(Z1, Z2, c, r):
        C = np.fft.fft2(g) / K ** 2
        sc = R ** -np.arange(K)
        gs.append((C * np.outer(sc, sc))[:D, :D])
    g0, g1, g2 = gs
    a = g0; b = G.shift_sum(g0) + g1
    cc = G.shift_sum(G.shift_sum(g0)) + 2 * G.shift_sum(g1) + g2
    sym = lambda h: 0.5 * (h + h.T).real
    return sym(a), sym(b), sym(cc)


def mats_r(c, th, r, N, D=36):
    a, b, cc = kap_arrays_r(c, r, D)
    N0 = N + D + 2
    B = G.Basis(th, N0, nq=N0 + 2 * r * r + 12)
    # weights X^q, q up to r^2+1
    B.W = [(B.V * (B.w * B.X ** q)[:, None]).T @ B.V for q in range(r * r + 2)]
    Ipow = [np.eye(N0)]
    for k in range(D):
        Ipow.append(B.I @ Ipow[-1])
    Ip = [P[:, :N] for P in Ipow]
    s0 = 2 * r - 1
    fac = [float(factorial(k)) for k in range(D + s0 + 2)]
    q0 = r * r - 1

    def Gm(p, h):
        out = np.zeros((N, N))
        WP = B.W[q0 + p]
        L = [WP @ Ip[k] for k in range(D)]
        for i in range(D):
            for j in range(D):
                hv = h[i, j]
                if hv == 0:
                    continue
                out += fac[i + s0] * fac[j + s0] * hv * (Ip[i].T @ L[j])
        return 0.5 * (out + out.T)
    M0 = Gm(0, a)
    MS = Gm(1, a) + Gm(0, b)
    MSS = Gm(2, a) + 2 * Gm(1, b) + Gm(0, cc)
    return M0, MS, MSS, B
