"""gapceil float core: Bui-Milinovich criterion (r=1) as an operator pencil in phi = I^2 f.
X-variables: c_j ∝ ∫ f(X+X1+X2) f(X+X3+X4) rho_j(X,X1..X4), rho_j = ∫_0^{A1}∫_0^{A2} (v-S)^{2j} U conj(U') dτ1 dτ2,
A1 = 1-X1-X3, A2 = 1-X2-X4, U conj U' = e^{ic(X1-X3)} + e^{ic(X2-X4)} - e^{ic(X1-X4)} e^{ic(τ1-τ2)} - e^{ic(X2-X3)} e^{-ic(τ1-τ2)}, c = kappa*pi,
S = X + X1+X2+X3+X4 + τ1+τ2.  f supported on [0, th]; f(U) = P(1-U/th).
Q_j = sum_{m,n>=1} m! n! ∫_0^th kap_mn(X) (I^{m-1}phi)(I^{n-1}phi) dX, phi(X) = ∫_X^th (u-X) f(u) du, I g(X) = ∫_X^th g.
"""
import numpy as np
from numpy.polynomial.legendre import leggauss
from math import factorial


def gl01(n):
    x, w = leggauss(n)
    return (x + 1) / 2, w / 2


def mk(A, z):
    if z == 0:
        return A, A ** 2 / 2, A ** 3 / 3
    e = np.exp(z * A)
    m0 = (e - 1) / z
    m1 = (A * e - m0) / z
    m2 = (A * A * e - 2 * m1) / z
    return m0, m1, m2


def rho_parts(X1, X2, X3, X4, c):
    """(rho0, rho1, rho2): ∫∫ U conj U' * (1, τ1+τ2 [J1], J2) ; rho_(j=1) = V^2 rho0 - 2V rho1 + rho2, V = v - S0."""
    A1 = 1 - X1 - X3
    A2 = 1 - X2 - X4
    out = [0, 0, 0]
    terms = [(1, X1 - X3, 0, 0), (1, X2 - X4, 0, 0), (-1, X1 - X4, 1j * c, -1j * c), (-1, X2 - X3, -1j * c, 1j * c)]
    for s, ph, z1, z2 in terms:
        a0, a1, a2 = mk(A1, z1)
        b0, b1, b2 = mk(A2, z2)
        e = s * np.exp(1j * c * ph)
        out[0] = out[0] + e * (a0 * b0)
        out[1] = out[1] + e * (a1 * b0 + a0 * b1)
        out[2] = out[2] + e * (a2 * b0 + 2 * a1 * b1 + a0 * b2)
    return out


def gbar(r1, r2, c, nth=24):
    """theta-averages of rho_parts with X1=r1 t1, X2=r1(1-t1), X3=r2 t2, X4=r2(1-t2); r1, r2 arrays (same shape), complex ok."""
    t, w = gl01(nth)
    T1, T2 = np.meshgrid(t, t, indexing='ij')
    W = np.outer(w, w)
    sh = r1.shape
    R1a = r1.reshape(-1, 1, 1); R2a = r2.reshape(-1, 1, 1)
    out = [np.zeros(R1a.shape[0], complex) for _ in range(3)]
    ch = 128
    for s in range(0, R1a.shape[0], ch):
        R1 = R1a[s:s + ch]; R2 = R2a[s:s + ch]
        p = rho_parts(R1 * T1, R1 * (1 - T1), R2 * T2, R2 * (1 - T2), c)
        for k in range(3):
            out[k][s:s + ch] = (p[k] * W).sum(axis=(1, 2))
    return [o.reshape(sh) for o in out]


def taylor2(c, D=36, R=0.4, K=64, nth=24):
    """2D Taylor coefficients [i,j] of r1^i r2^j for g0,g1,g2 (i,j < D) via FFT on |r|=R."""
    ang = np.exp(2j * np.pi * np.arange(K) / K)
    Z1, Z2 = np.meshgrid(R * ang, R * ang, indexing='ij')
    gs = []
    for g in gbar(Z1, Z2, c, nth):
        C = np.fft.fft2(g) / K ** 2
        sc = R ** -np.arange(K)
        C = C * np.outer(sc, sc)
        gs.append(C[:D, :D])
    return gs


def shift_sum(h):
    """coefficients of (r1+r2) h"""
    out = np.zeros_like(h)
    out[1:, :] += h[:-1, :]
    out[:, 1:] += h[:, :-1]
    return out


def kap_arrays(c, D=36, **kw):
    g0, g1, g2 = taylor2(c, D, **kw)
    a = g0
    b = shift_sum(g0) + g1
    cc = shift_sum(shift_sum(g0)) + 2 * shift_sum(g1) + g2
    sym = lambda h: 0.5 * (h + h.T).real
    return sym(a), sym(b), sym(cc)


class Basis:
    """L2-orthonormal shifted Legendre on [0, th]; coefficient space dimension N0."""
    def __init__(self, th, N0, nq=None):
        self.th, self.N0 = th, N0
        nq = nq or N0 + 8
        x, w = leggauss(nq)
        self.X = th * (x + 1) / 2
        self.w = w * th / 2
        from numpy.polynomial.legendre import legvander
        V = legvander(x, N0 - 1)
        self.norm = np.sqrt((2 * np.arange(N0) + 1) / th)
        self.V = V * self.norm      # values of orthonormal basis at nodes
        # tail-integration I in coefficient space of the *unnormalized* P_n(2X/th-1): ∫_X^th P_n = th/2 (P_{n-1}-P_{n+1})/(2n+1), n>=1; n=0: th/2 (P0 - P1)
        I = np.zeros((N0, N0))
        I[0, 0] += th / 2; I[1, 0] -= th / 2
        for n in range(1, N0):
            I[n - 1, n] += th / 2 / (2 * n + 1)
            if n + 1 < N0:
                I[n + 1, n] -= th / 2 / (2 * n + 1)
        # convert to orthonormal coefficients: phi = sum a_n hat_n, hat_n = norm_n P_n ; I hat_n = norm_n I P_n
        self.I = (I * self.norm[None, :]) / self.norm[:, None]
        self.W = [ (self.V * (self.w * self.X ** p)[:, None]).T @ self.V for p in range(3)]


def mats(c, th, N, D=36, N0=None, **kw):
    """M0, MS, MSS (N x N) in the orthonormal Legendre basis of phi on [0, th]."""
    a, b, cc = kap_arrays(c, D, **kw)
    N0 = N0 or N + D + 2
    B = Basis(th, N0)
    Ipow = [np.eye(N0)]
    for k in range(D):
        Ipow.append(B.I @ Ipow[-1])
    Ip = [P[:, :N] for P in Ipow]
    fac = [float(factorial(k)) for k in range(D + 2)]

    def G(p, h):
        out = np.zeros((N, N))
        WP = B.W[p]
        L = [WP @ Ip[k] for k in range(D)]
        for m in range(1, D + 1):
            for n in range(1, D + 1):
                hv = h[m - 1, n - 1]
                if hv == 0:
                    continue
                out += fac[m] * fac[n] * hv * (Ip[m - 1].T @ L[n - 1])
        return 0.5 * (out + out.T)
    M0 = G(0, a)
    MS = G(1, a) + G(0, b)
    MSS = G(2, a) + 2 * G(1, b) + G(0, cc)
    return M0, MS, MSS, B


def top(M0, MS, MSS, v):
    M1 = v * v * M0 - 2 * v * MS + MSS
    L = np.linalg.cholesky(M1)
    Li = np.linalg.inv(L)
    ev, U = np.linalg.eigh(Li @ M0 @ Li.T)
    return ev[-1], Li.T @ U[:, -1]


def best_v(M0, MS, MSS, v0=1.3, lo=0.8, hi=2.0):
    from scipy.optimize import minimize_scalar
    r = minimize_scalar(lambda v: -top(M0, MS, MSS, v)[0], bounds=(lo, hi), method='bounded', options={'xatol': 1e-9})
    return -r.fun, r.x
