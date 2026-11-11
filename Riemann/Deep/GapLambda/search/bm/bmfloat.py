"""Float estimate for Bui-Milinovich (arXiv 1410.3635v3) Thm smoothed (l.300-319): h1 = c0/(kappa^2 c1) (l.327-330).
c_j(v,kappa) = int_{[0,1]^9, x+x1+x2<=1, x+x3+x4<=1} e^{i Phi} a b (th(x1-x2)+a t1-b t2)(th(x3-x4)+a t1-b t2)
               x^{r^2-1}(x1x2x3x4)^{r-1} (v - S)^{2j} P(1-x-x1-x2) P(1-x-x3-x4),
a = 1-th(x1+x3), b = 1-th(x2+x4), S = th(x+x1+x2+x3+x4) + a t1 + b t2,
Phi = th k pi (x2-x4-(x3-x4)t3+(x1-x2)t4) - k pi (a t1 - b t2)(t3-t4)   (l.307).
t3, t4 enter only through e^{i Phi}: int_0^1 int_0^1 = e^{i th k pi (x2-x4)} E(w3) E(w4), E(w) = (e^{iw}-1)/(iw),
w3 = -th k pi (x3-x4) - k pi (a t1 - b t2), w4 = th k pi (x1-x2) + k pi (a t1 - b t2).
Remaining 7-dim integral: Gauss-Legendre, collapsed coordinates x1 = R s1, x2 = R(1-s1)e1 (R = 1-x), same for x3, x4.
P in the shifted-Legendre basis of w on [0,1]. Returns real symmetric matrices M0, M1, M2 with
c0 = P'M0P, c1(v) = P'(v^2 M0 - 2v M1 + M2)P."""
import numpy as np
from numpy.polynomial.legendre import leggauss, Legendre


def gl01(n):
    x, w = leggauss(n); return (x + 1) / 2, w / 2


def leg01(D, w):
    """(len(w), D) values of shifted Legendre polys on [0,1]."""
    out = np.empty((len(w), D))
    for m in range(D):
        out[:, m] = Legendre.basis(m, domain=[0, 1])(w)
    return out


def E(w):
    return np.exp(0.5j * w) * np.sinc(w / (2 * np.pi))


def mats(th, r, kap, D, n):
    z, wz = gl01(n)
    S1, E1, S2, E2 = np.meshgrid(z, z, z, z, indexing='ij')
    W4 = np.einsum('i,j,k,l->ijkl', wz, wz, wz, wz)
    S1, E1, S2, E2, W4 = S1.ravel(), E1.ravel(), S2.ravel(), E2.ravel(), W4.ravel()
    kp = kap * np.pi
    M = [np.zeros((D, D), complex) for _ in range(3)]
    for ix, x in enumerate(z):
        R = 1 - x
        x1 = R * S1; x2 = R * (1 - S1) * E1; x3 = R * S2; x4 = R * (1 - S2) * E2
        jac = (R ** 2 * (1 - S1)) * (R ** 2 * (1 - S2)) * W4 * wz[ix]
        a = 1 - th * (x1 + x3); b = 1 - th * (x2 + x4)
        wt = jac * a * b * x ** (r * r - 1) * (x1 * x2 * x3 * x4) ** (r - 1) * np.exp(1j * th * kp * (x2 - x4))
        F1 = leg01(D, 1 - x - x1 - x2); F2 = leg01(D, 1 - x - x3 - x4)
        base0 = th * (x1 + x2 + x3 + x4 + x)
        for it1, t1 in enumerate(z):
            for it2, t2 in enumerate(z):
                s = a * t1 - b * t2
                g = wt * wz[it1] * wz[it2] * (th * (x1 - x2) + s) * (th * (x3 - x4) + s)
                g = g * E(-th * kp * (x3 - x4) - kp * s) * E(th * kp * (x1 - x2) + kp * s)
                Sv = base0 + a * t1 + b * t2
                for k in range(3):
                    M[k] += (F1 * g[:, None]).T @ F2
                    g = g * Sv
    out = []
    for A in M:
        Ar = A.real
        out.append(0.5 * (Ar + Ar.T))
    imag = max(np.abs(0.5 * (A.imag + A.imag.T)).max() / np.abs(A.real).max() for A in M)
    return out, imag


def h1max(M, kap, v):
    """max over P of c0/(kappa^2 c1) for given v: largest generalised eigenvalue."""
    M0, M1, M2 = M
    A1 = v * v * M0 - 2 * v * M1 + M2
    L = np.linalg.cholesky(A1)
    Li = np.linalg.inv(L)
    ev, evec = np.linalg.eigh(Li @ M0 @ Li.T)
    P = Li.T @ evec[:, -1]
    return ev[-1] / kap ** 2, P / P[0]


def h1P(M, kap, v, Pleg):
    M0, M1, M2 = M
    A1 = v * v * M0 - 2 * v * M1 + M2
    return (Pleg @ M0 @ Pleg) / (kap ** 2 * (Pleg @ A1 @ Pleg))
