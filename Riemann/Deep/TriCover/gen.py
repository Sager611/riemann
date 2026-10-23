#!/usr/bin/env python3
"""Generate an optimized cyclic chain (arXiv 2609.06558, Sections 5-6) and its exact rational certificate.
config.txt is the output of 'python3 gen.py 2000' (renamed).

Needs numpy, scipy, mpmath.  Usage: python3 gen.py M [EXP]   (writes config_M.txt; about 15 s for M = 2000)
  M   = bridges per chain (family has 3(M+1) members)
  EXP = coordinates are written as integers over 10^EXP (default 15)

Steps
 1. Parameters c, d, tau, h_1..h_M of the paper's recurrence (Section 5: support (c, d, 0),
    last vertex on the conic parametrised by tau, steps h_k).  Minimise
    S = 3(1 - c - d + sum_k h_k E_k) subject to the two product equations of Theorem 5.3, by
    SLSQP at M = 100 (start: the paper's continuum candidate, uniform steps), then an
    augmented Lagrangian with L-BFGS-B along M = 200, 300, 500, 1000, 2000, ... up to the
    target (step profile interpolated from the previous M).
 2. With mpmath (60 digits) re-solve the two product equations exactly for (tau, a common
    step factor), run the recurrence, round every nonzero coordinate down to the grid
    10^-EXP and subtract one more unit.  Rounding down enlarges every triangle, so each
    vertex contact of the recurrence becomes a strict segment relation; verify.py checks
    this exactly.
"""
import sys
import numpy as np
from scipy.optimize import minimize
import mpmath as mp

def params(c, d, tau):
    den = 1 + tau + tau * tau
    v = (1 - d + d * tau) / den; u = tau * v; kap = 1 - u - v; D = v - d
    E0 = d * (kap + D) - kap * c; R0 = (v - c) / E0; Y0 = d / E0
    Ef = kap * (v - u) + D * (v - kap); Rf = (v - u) / Ef; Yf = v / Ef
    return u, v, kap, D, E0, R0, Y0, Ef, Rf, Yf

def chain(c, d, tau, h):
    u, v, kap, D, E0, R0, Y0, Ef, Rf, Yf = params(c, d, tau)
    R = R0 * np.cumprod(1 + D * h); Y = Y0 * np.cumprod(1 - kap * h)
    E = D * kap / (kap * R + D * Y - 1)
    return R, Y, E

def S_of(z):
    R, Y, E = chain(z[0], z[1], z[2], z[3:])
    return 3 * (1 - z[0] - z[1] + np.sum(z[3:] * E))

def eq(z):
    c, d, tau = z[:3]; h = z[3:]
    u, v, kap, D, E0, R0, Y0, Ef, Rf, Yf = params(c, d, tau)
    return np.array([np.sum(np.log1p(D * h)) - np.log(Rf / R0), np.sum(np.log1p(-kap * h)) - np.log(Yf / Y0)])

def ineq(z):
    c, d, tau = z[:3]
    u, v, kap, D, E0, R0, Y0, Ef, Rf, Yf = params(c, d, tau)
    return np.array([1 - u - c - d, E0, v - c, d - v, u + v - d])

def fd3(f, z):
    cols = []
    for i in range(3):
        zp = z.copy(); zm = z.copy(); zp[i] += 1e-7; zm[i] -= 1e-7
        cols.append((f(zp) - f(zm)) / 2e-7)
    return cols

def gradS(z):
    c, d, tau = z[:3]; h = z[3:]
    u, v, kap, D = params(c, d, tau)[:4]
    R, Y, E = chain(c, d, tau, h)
    sa = np.cumsum((h * E * E * R)[::-1])[::-1]; sb = np.cumsum((h * E * E * Y)[::-1])[::-1]
    g = np.empty_like(z)
    g[3:] = 3 * (E - sa / (1 + D * h) + sb / (1 - kap * h))
    g[:3] = fd3(S_of, z)
    return g

def jac_eq(z):
    c, d, tau = z[:3]; h = z[3:]
    u, v, kap, D = params(c, d, tau)[:4]
    J = np.zeros((2, len(z)))
    J[0, 3:] = D / (1 + D * h); J[1, 3:] = -kap / (1 - kap * h)
    J[:, :3] = np.array(fd3(eq, z)).T
    return J

def slsqp(m):
    c, d, tau, s = 0.3104341279739202, 0.41586782506519804, 0.6737630444333738, 1.1213014301809736
    z0 = np.concatenate([[c, d, tau], np.full(m, s / m)])
    res = minimize(S_of, z0, jac=gradS, method='SLSQP',
                   bounds=[(0.2, 0.45), (0.3, 0.5), (0.3, 0.95)] + [(1e-9, 1.0)] * m,
                   constraints=[{'type': 'eq', 'fun': eq, 'jac': jac_eq}, {'type': 'ineq', 'fun': ineq}],
                   options={'maxiter': 500, 'ftol': 1e-16})
    return res.x

def auglag(z, outer=30, mu=1e4):
    lam = np.zeros(2); m = len(z) - 3
    bounds = [(0.2, 0.45), (0.3, 0.5), (0.3, 0.95)] + [(1e-12, 1.0)] * m
    for it in range(outer):
        def f(w):
            g = eq(w)
            return S_of(w) + lam @ g + 0.5 * mu * g @ g, gradS(w) + jac_eq(w).T @ (lam + mu * g)
        z = minimize(f, z, jac=True, method='L-BFGS-B', bounds=bounds,
                     options={'maxiter': 5000, 'maxcor': 30, 'ftol': 1e-18, 'gtol': 1e-14}).x
        g = eq(z); lam = lam + mu * g
        if np.max(np.abs(g)) < 1e-15 and it > 3:
            break
        mu = min(mu * 3, 1e9)
    return z

def interp(z, m):
    h = z[3:]; mo = len(h)
    xo = (np.arange(1, mo + 1) - 0.5) / mo; x = (np.arange(1, m + 1) - 0.5) / m
    return np.concatenate([z[:3], np.interp(x, xo, h) * mo / m])

LADDER = [100, 200, 300, 500, 1000, 2000, 3000, 5000, 10000, 20000, 50000]

def optimise(M):
    z = slsqp(min(M, 100))
    for m in [k for k in LADDER if 100 < k < M] + ([M] if M > 100 else []):
        z = auglag(interp(z, m))
    return z

def certificate(z, EXP):
    mp.mp.dps = 60
    c, d, tau0 = [mp.mpf(float(v)) for v in z[:3]]
    h0 = [mp.mpf(float(v)) for v in z[3:]]
    def resid(tau, al):
        u, v, kap, D, E0, R0, Y0, Ef, Rf, Yf = params(c, d, tau)
        f = mp.e ** al
        return [mp.fsum(mp.log(1 + D * h * f) for h in h0) - mp.log(Rf / R0),
                mp.fsum(mp.log(1 - kap * h * f) for h in h0) - mp.log(Yf / Y0)]
    x = [tau0, mp.mpf(0)]
    for _ in range(40):
        r = resid(*x)
        if max(abs(r[0]), abs(r[1])) < mp.mpf(10) ** -50:
            break
        J = mp.matrix(2, 2); eps = mp.mpf(10) ** -25
        for j in range(2):
            xp = list(x); xp[j] += eps; rp = resid(*xp)
            J[0, j] = (rp[0] - r[0]) / eps; J[1, j] = (rp[1] - r[1]) / eps
        dx = mp.lu_solve(J, mp.matrix(r)); x = [x[0] - dx[0], x[1] - dx[1]]
    tau, al = x
    u, v, kap, D, E0, R0, Y0, Ef, Rf, Yf = params(c, d, tau)
    R, Y = R0, Y0
    pts = [(c, d, mp.mpf(0))]
    for hk in h0:
        hk = hk * mp.e ** al
        R *= 1 + D * hk; Y *= 1 - kap * hk
        E = D * kap / (kap * R + D * Y - 1)
        xk = v - E * R; yk = E * Y; tk = hk * E
        pts.append((xk, yk, 1 - xk - yk - tk))
    N = 10 ** EXP
    return N, [[0 if q == 0 else int(mp.floor(q * N)) - 1 for q in p] for p in pts]

if __name__ == '__main__':
    M = int(sys.argv[1]); EXP = int(sys.argv[2]) if len(sys.argv) > 2 else 15
    z = optimise(M)
    print('M = %d: float S = %.13f, lambda = %.13f' % (M, S_of(z), 1 / S_of(z)))
    N, rows = certificate(z, EXP)
    with open('config_%d.txt' % M, 'w') as f:
        f.write('# cyclic NS family of 3(M+1) = %d positive homothetic triangles, M = %d\n' % (3 * (M + 1), M))
        f.write('# p^(k) = (X_k, Y_k, Z_k)/DEN, k = 0..M (k = 0: support); members T(rho^i p^(k)), i = 0,1,2\n')
        f.write('DEN %d\n' % N)
        for k, r in enumerate(rows):
            f.write('%d %d %d %d\n' % (k, r[0], r[1], r[2]))
    print('wrote config_%d.txt' % M)

