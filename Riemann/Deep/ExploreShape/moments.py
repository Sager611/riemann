"""explore53: higher spectral moments <S^k> of BM's main terms (midpoint mollifier, general r), k <= K.
Spectral picture: c_j(v) = <(v - S)^{2j}> is the mean of |F^{(j)}|^2 (u = tL units), F = e^{ivtL} zeta(s) zeta(s+i kappa pi/L) M(s+i kappa pi/2L).
Reuses GapBar's gcore/gcore_r (read-only import); rho_parts extended to tau-moments (tau1+tau2)^m, m <= K.
usage: moments.py r theta d kappa K  -> logs/mom_r{r}_d{d}_k{kappa}.json (raw moments s_k of the optimal amplifier, v, h1)"""
import sys, os, json
from math import comb, factorial
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'GapBar'))
import gcore as G


def mkn(A, z, K):
    if z == 0:
        return [A ** (n + 1) / (n + 1) for n in range(K + 1)]
    e = np.exp(z * A)
    m = [(e - 1) / z]
    for n in range(1, K + 1):
        m.append((A ** n * e - n * m[-1]) / z)
    return m


def rho_K(X1, X2, X3, X4, c, K):
    A1 = 1 - X1 - X3; A2 = 1 - X2 - X4
    out = [0] * (K + 1)
    terms = [(1, X1 - X3, 0, 0), (1, X2 - X4, 0, 0), (-1, X1 - X4, 1j * c, -1j * c), (-1, X2 - X3, -1j * c, 1j * c)]
    for s, ph, z1, z2 in terms:
        a = mkn(A1, z1, K); b = mkn(A2, z2, K)
        e = s * np.exp(1j * c * ph)
        for m in range(K + 1):
            out[m] = out[m] + e * sum(comb(m, i) * a[i] * b[m - i] for i in range(m + 1))
    tw = np.exp(1j * c * (X3 + X4 - X1 - X2) / 2)   # midpoint mollifier twist (GapBar optmid.py)
    return [o * tw for o in out]


def gbar_K(r1, r2, c, r, K, nth=28):
    t, w = G.gl01(nth)
    wr = w * (t * (1 - t)) ** (r - 1)
    T1, T2 = np.meshgrid(t, t, indexing='ij')
    W = np.outer(wr, wr)
    sh = r1.shape
    R1a = r1.reshape(-1, 1, 1); R2a = r2.reshape(-1, 1, 1)
    out = [np.zeros(R1a.shape[0], complex) for _ in range(K + 1)]
    ch = 64
    for s in range(0, R1a.shape[0], ch):
        R1 = R1a[s:s + ch]; R2 = R2a[s:s + ch]
        p = rho_K(R1 * T1, R1 * (1 - T1), R2 * T2, R2 * (1 - T2), c, K)
        for k in range(K + 1):
            out[k][s:s + ch] = (p[k] * W).sum(axis=(1, 2))
    return [o.reshape(sh) for o in out]


def arrays_K(c, r, K, D=40, R=0.4, NF=64):
    ang = np.exp(2j * np.pi * np.arange(NF) / NF)
    Z1, Z2 = np.meshgrid(R * ang, R * ang, indexing='ij')
    gs = []
    for g in gbar_K(Z1, Z2, c, r, K):
        C = np.fft.fft2(g) / NF ** 2
        sc = R ** -np.arange(NF)
        gs.append((C * np.outer(sc, sc))[:D, :D])
    H = []
    for n in range(K + 1):
        h = np.zeros_like(gs[0])
        for j in range(n + 1):
            x = gs[n - j]
            for _ in range(j):
                x = G.shift_sum(x)
            h = h + comb(n, j) * x
        H.append(h)
    sym = lambda h: 0.5 * (h + h.T)
    return [sym(h) for h in H]


def mats_K(c, th, r, N, K, D=40):
    H = arrays_K(c, r, K, D)
    imax = max(np.abs(h.imag).max() / max(np.abs(h.real).max(), 1e-300) for h in H)
    H = [h.real for h in H]
    N0 = N + D + 2
    B = G.Basis(th, N0, nq=N0 + 2 * r * r + 12 + K)
    B.W = [(B.V * (B.w * B.X ** q)[:, None]).T @ B.V for q in range(r * r + K + 1)]
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
    Ms = []
    for k in range(K + 1):
        Ms.append(sum(comb(k, i) * Gm(i, H[k - i]) for i in range(k + 1)))
    return Ms, B, imax


def setup(r, th, d, kap, K):
    Nb = d + 2 * r + 6
    Ms, B, imax = mats_K(kap * np.pi, th, r, Nb, K)
    Ir = np.linalg.matrix_power(B.I, 2 * r)[:Nb, :d + 1]
    Q, R = np.linalg.qr(Ir)
    return [Q.T @ A @ Q for A in Ms], imax


if __name__ == '__main__':
    r = int(sys.argv[1]); th = float(sys.argv[2]); d = int(sys.argv[3]); kap = float(sys.argv[4]); K = int(sys.argv[5])
    Ms, imax = setup(r, th, d, kap, K)
    lam, v = G.best_v(Ms[0], Ms[1], Ms[2], lo=0.8, hi=2.0)
    ev, y = G.top(Ms[0], Ms[1], Ms[2], v)
    s = [float(y @ M @ y) for M in Ms]
    s = [x / s[0] for x in s]
    mu = s[1]
    cen = [sum(comb(k, i) * s[i] * (-mu) ** (k - i) for i in range(k + 1)) for k in range(K + 1)]
    out = {'r': r, 'theta': th, 'd': d, 'kappa': kap, 'K': K, 'v': v, 'h1': lam / kap ** 2, 'raw': s, 'central_about_mean': cen,
           'mean': mu, 'imag_rel_max': imax}
    # Hankel check of the raw moments (must be PSD: they are Gram data of F, F', ...)
    J = K // 2
    Hk = np.array([[s[i + j] for j in range(J + 1)] for i in range(J + 1)])
    out['hankel_eigs'] = [float(x) for x in np.linalg.eigvalsh(Hk)]
    sd = cen[2] ** 0.5
    out['standardized'] = [cen[k] / sd ** k for k in range(K + 1)]
    fn = os.path.join(HERE, 'logs', 'mom_r%d_d%d_k%s.json' % (r, d, sys.argv[4]))
    json.dump(out, open(fn, 'w'), indent=1)
    print('r=%d th=%.4f d=%d kappa=%.4f: h1=%.8f v=%.6f mean S=%.6f var=%.6g (1/kappa^2=%.6g) imag_rel=%.1e' % (r, th, d, kap, lam / kap ** 2, v, mu, cen[2], 1 / kap ** 2, imax))
    print('standardized central moments:', ' '.join('%.5f' % x for x in out['standardized']))
    print('hankel eigs:', ' '.join('%.3e' % x for x in out['hankel_eigs']))
