"""explore53: third-order local Wirtinger criterion (J = 3) for Bui-Milinovich's F, joint with the amplifier.
Per-gap inequality (normalised units: gaps <= pi, frequency omega' = kappa*omega): for real or complex f with f(0)=f(g)=0,
g <= pi:   |f|^2 <= gam |f'|^2 + a3 |f''' + mu f'|^2                                              (*)
Exact validity (sine coefficients b_k of f, A = f''(0), B = f''(g) free; REPORT.md Lemma 1):
  (i) p(w) = -1 + gam w^2 + a3 w^2 (w^2 - mu)^2 > 0 for w >= 1;
  (ii) S_odd(g) = sum_{k odd} r_k >= 0, (iii) S_even(g) + 1/2 >= 0, r_k = (gam w_k^2 - 1)/p(w_k), w_k = k pi/g, all g <= pi.
Globally (*) summed over the gaps of F gives c0 <= gam c1 + a3 (c3 - 2 mu c2 + mu^2 c1) (Re int F''' conj F' = -int |F''|^2),
c_j = <(v-S)^{2j}> kappa^{2j}.  h3 = lambda_max(M0, gam M2 + a3(M6 - 2mu M4 + mu^2 M2)) > 1 => Lambda > kappa.
usage: crit3.py r theta d kappa  (uses moments.setup; prints h1 (Wirtinger) and the optimised h3)"""
import sys, os, json
from math import comb
import numpy as np
from scipy.optimize import minimize, minimize_scalar
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)


def pval(w2, gam, a3, mu):
    return -1 + gam * w2 + a3 * w2 * (w2 - mu) ** 2


def cond_i(gam, a3, mu):
    xs = [1.0]
    disc = mu * mu - 3 * gam / a3 if a3 > 0 else -1
    if disc >= 0:
        for sg in (1, -1):
            x = (2 * mu + sg * np.sqrt(disc)) / 3
            if x > 1:
                xs.append(x)
    return min(pval(x, gam, a3, mu) for x in xs)


KM = 4000
KK = np.arange(1, KM + 1)


def sums(g, gam, a3, mu):
    w2 = (KK * np.pi / g) ** 2
    p = pval(w2, gam, a3, mu)
    r = (gam * w2 - 1) / p
    # tail k > KM: r_k ~ gam/(a3 w^4) > 0 (dropping it is conservative for (ii),(iii))
    return r[0::2].sum(), r[1::2].sum() + 0.5, p.min()


GGRID = np.pi * np.concatenate([np.linspace(0.05, 1, 400)])


def valid(gam, a3, mu):
    if a3 < 0:
        return False
    if cond_i(gam, a3, mu) <= 0:
        return False
    for g in GGRID:
        so, se, pm = sums(g, gam, a3, mu)
        if so < 0 or se < 0 or pm <= 0:
            return False
    return True


def gam_min(a3, mu, lo=0.5, hi=1.0):
    if not valid(hi, a3, mu):
        return None
    for _ in range(40):
        m = 0.5 * (lo + hi)
        if valid(m, a3, mu):
            hi = m
        else:
            lo = m
    return hi


def central(Ms, v, kap, jmax=3):
    out = []
    for j in range(jmax + 1):
        n = 2 * j
        A = sum(comb(n, i) * (-v) ** (n - i) * Ms[i] for i in range(n + 1))
        out.append(A * kap ** n)
    return out


def lmax(A, B):
    L = np.linalg.cholesky(B)
    Li = np.linalg.inv(L)
    return np.linalg.eigvalsh(Li @ A @ Li.T)[-1]


def h1(Ms, v, kap):
    C = central(Ms, v, kap, 1)
    return lmax(C[0], C[1])


def h3(Ms, v, kap, a3, mu, gam=None):
    gam = gam if gam is not None else gam_min(a3, mu)
    if gam is None:
        return 0.0, None
    C = central(Ms, v, kap, 3)
    R = gam * C[1] + a3 * (C[3] - 2 * mu * C[2] + mu * mu * C[1])
    return lmax(C[0], R), gam


if __name__ == '__main__':
    import moments as Mo
    import gcore as G
    r = int(sys.argv[1]); th = float(sys.argv[2]); d = int(sys.argv[3]); kap = float(sys.argv[4])
    Ms, _ = Mo.setup(r, th, d, kap, 6)
    lam, v1 = G.best_v(Ms[0], Ms[1], Ms[2], lo=0.8, hi=2.0)
    print('Wirtinger: h1 = %.8f at v = %.6f' % (h1(Ms, v1, kap), v1), flush=True)
    best = [0, None]
    def f(x):
        a3, mu, dv = np.exp(x[0]), x[1], x[2]
        val, gam = h3(Ms, v1 + dv, kap, a3, mu)
        if val > best[0]:
            best[0] = val; best[1] = (a3, mu, v1 + dv, gam)
        return -val
    for a30, mu0 in [(0.00137, 4.84), (0.003, 4.0), (0.0007, 6.0)]:
        res = minimize(f, [np.log(a30), mu0, 0.0], method='Nelder-Mead', options={'xatol': 1e-5, 'fatol': 1e-9, 'maxiter': 300})
        print('  start (%.4g, %.3g): h3 = %.8f' % (a30, mu0, -res.fun), flush=True)
    a3, mu, v, gam = best[1]
    print('J=3: h3 = %.8f  a3=%.6g mu=%.6g v=%.6f gam=%.8f' % (best[0], a3, mu, v, gam))
    json.dump({'r': r, 'theta': th, 'd': d, 'kappa': kap, 'h1': h1(Ms, v1, kap), 'v1': v1, 'h3': best[0], 'a3': a3, 'mu': mu,
               'v3': v, 'gam': gam}, open(os.path.join(HERE, 'logs', 'crit3_r%d_d%d_k%s.json' % (r, d, sys.argv[4])), 'w'), indent=1)
