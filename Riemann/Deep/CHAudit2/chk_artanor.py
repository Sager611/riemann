# chaudit2: end-to-end numerical test of Lemma artanor (l.5633) of 2511.14736v1:
# ||e^{-a u} f||_{L1[u0,oo)} <= c0(a) + e^{-a u0}/(11 pi a u0^2),  ||e^{-a u} f||_{TV[u0,oo)} <= c1(a) + c2(a) e^{-a u0}/(a u0^2),
# f = hat(phi_lambda) - I_lambda = K_nu(v) - E_nu(v), v = 2 sgn(lambda) u, nu = |lambda|/2, via eq. monthmayending.
import numpy as np, mpmath as mp
from numpy.polynomial.laguerre import laggauss
from numpy.polynomial.legendre import leggauss

xl, wl = laggauss(120)
xg, wg = leggauss(12)
# composite Gauss-Legendre on [0, 120] with 600 panels
edges = np.linspace(0, 120, 601)
W_nodes = ((edges[:-1, None] + edges[1:, None]) / 2 + (edges[1:, None] - edges[:-1, None]) / 2 * xg[None, :]).ravel()
W_wts = (((edges[1:] - edges[:-1]) / 2)[:, None] * wg[None, :]).ravel()

def mbp(y):  # -b'(y) = 1/(4 cosh^2(y/2))
    return 0.25 / np.cosh(y / 2) ** 2

def G(v, nu):
    v = np.asarray(v, float); out = np.empty_like(v)
    big = np.abs(v) >= 0.5
    vb = v[big]; av = np.abs(vb)
    w = xl[None, :] / av[:, None]
    neg = vb < 0
    h = np.where(neg[:, None], mbp(nu + w), mbp(nu - w))
    val = (h * wl[None, :]).sum(1) / av
    out[big] = np.where(neg, val, -val)
    vs = v[~big]
    if vs.size:
        e = np.exp(-np.abs(vs)[:, None] * W_nodes[None, :])
        neg = vs < 0
        h = np.where(neg[:, None], mbp(nu + W_nodes[None, :]), mbp(nu - W_nodes[None, :]))
        val = (h * e * W_wts[None, :]).sum(1)
        out[~big] = np.where(neg, val, -val)
    return out

def D(v, nu):
    return np.sinc(v) * G(v, nu)   # np.sinc(v) = sin(pi v)/(pi v)

# cross-check D against the series definition K_nu(z) = sin(pi z)/pi * sum_{n>=1} (-1)^n e^{-nu n}(1/(z-n) - 1/z), E_nu
mp.mp.dps = 20
def D_series(v, nu):
    q = -mp.e**(-nu)
    K = mp.sin(mp.pi*v)/mp.pi * (q*mp.lerchphi(q, 1, 1 - v)*(-1) + 1/(v*(mp.e**nu + 1)))
    E = mp.e**(-nu*v) if v > 0 else 0
    return K - E
mx = 0
for nu in [0.0005, 0.25, 1.5]:
    for v in [-37.3, -3.7, -1.2, -0.45, -0.05, 0.03, 0.4, 0.77, 2.6, 15.1, 80.4]:
        a = D(np.array([v]), nu)[0]; b = float(D_series(mp.mpf(v), mp.mpf(nu)))
        mx = max(mx, abs(a - b))
print("max |D_quadrature - D_series| on 33 test points =", mx)

def c012(a):
    a = mp.mpf(a); k = 2*mp.e**(a/2) - a*mp.ei(a/2) + mp.mpf('2.525')*a; p = mp.pi
    c0 = mp.e**(a/2)/2 + mp.e**(-a/2)/(8*p) + k/(16*p)
    c1 = 2*mp.e**(a/2) + (a/(8*p) + mp.mpf(3)/11)*mp.e**(-a/2) + (a/(16*p) + mp.mpf(3)/22)*k
    c2 = mp.mpf(16)/11*(a/(16*p) + mp.mpf(3)/22)
    return float(c0), float(c1), float(c2)

CH = 10000
def grid(u0, U):
    parts = []
    a1 = min(max(u0, -5.0), 5.0)
    if u0 < -5.0: parts.append(np.arange(u0, -5.0, 0.004))
    parts.append(np.arange(a1, 5.0, 0.001))
    parts.append(np.arange(5.0, U + 1e-9, 0.004))
    g = np.concatenate(parts)
    return np.where(np.abs(g) < 1e-12, 1e-12, g)
worst_L1 = (0,); worst_TV = (0,)
rows = []
for lam in [-6.0, -0.5, -1e-4, 1e-4, 0.5, 6.0]:
    nu = abs(lam) / 2; sg = np.sign(lam)
    for a in [0.02, 1/6, 0.5, 2.0]:
        U = 1500.0 if a < 0.1 else 300.0
        for u0 in [-0.5, -1.0, -3.0, -10.0, -40.0, -0.3, 0.2]:
            if a * abs(u0) > 60: continue
            us = grid(u0, U)
            F = np.empty_like(us)
            for i in range(0, us.size, CH):
                F[i:i+CH] = np.exp(-a * us[i:i+CH]) * D(2 * sg * us[i:i+CH], nu)
            L1 = np.sum(0.5 * (np.abs(F[1:]) + np.abs(F[:-1])) * np.diff(us)) + np.exp(-a*U)/(16*np.pi*U)
            TV = np.sum(np.abs(np.diff(F))) + (a/(16*np.pi) + 3/22) * np.exp(-a*U) / U
            c0, c1, c2 = c012(a)
            extra = np.exp(-a*u0)/(a*u0**2)
            rL = L1 / (c0 + extra/(11*np.pi)); rT = TV / (c1 + c2*extra)
            rows.append((lam, a, u0, L1, rL, TV, rT))
            if rL > worst_L1[0]: worst_L1 = (rL, lam, a, u0)
            if rT > worst_TV[0]: worst_TV = (rT, lam, a, u0)
            print("lam=%-8g a=%-6.4g u0=%-6g L1=%.6f (ratio %.4f)  TV=%.6f (ratio %.4f)" % rows[-1], flush=True)
print("cases:", len(rows), " max L1 ratio:", worst_L1, " max TV ratio:", worst_TV)
