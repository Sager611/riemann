"""Float64 evaluator of c, c1(u), c2(u) for the optimiser (exact tables from core.py, then float)."""
import sys
sys.path.insert(0, '.')
import numpy as np
from math import comb
from flint import fmpq
import core


def G_float_table(r, n0, nl):
    """Float table G[n1,n2,n3,n4] (offset n0) of core.G via log-factorials; all summands positive (no cancellation)."""
    from math import lgamma
    e = r * r - 1
    nmax = n0 + nl - 1
    M = 4 * nmax + 4 * e + 16
    lf = np.array([lgamma(m + 1) for m in range(M)])          # lf[m] = log m!
    lB = lambda a, b: lf[a - 1] + lf[b - 1] - lf[a + b - 1]    # log Beta(a,b), a,b >= 1 integers (arrays ok)
    nv = np.arange(n0, nmax + 1)
    H = np.zeros((nl, nl, nl, nl))
    for i2, n2 in enumerate(nv):
        for i3, n3 in enumerate(nv):
            j = np.arange(n3 + 1)[None, None, :, None]
            k = np.arange(n2 + 1)[None, None, None, :]
            n1 = nv[:, None, None, None]
            n4 = nv[None, :, None, None]
            Q = n3 - j + n2 - k
            L = (lf[n3] - lf[j] - lf[n3 - j] + lf[n2] - lf[k] - lf[n2 - k] + lB(Q + 1, e + 1)
                 + lB(n1 + j + 1, e + 1) + lB(n4 + k + 1, e + 1) + lB(n1 + j + n4 + k + 2 * e + 3, 2 * e + Q + 2))
            H[:, i2, i3, :] = np.exp(L).sum(axis=(2, 3))
    G = H + H.transpose(2, 3, 0, 1)
    return G.reshape(nl * nl, nl * nl)


class Fast:
    def __init__(self, theta, r, dmax, floatG=True):
        self.r, self.dmax = r, dmax
        parts = core.poly_parts(fmpq(theta))
        self.terms = {k: core.split_terms(p) for k, p in parts.items()}
        T = self.terms
        self.gAs = sorted({k[0] for t in T.values() for k in t})
        self.gBs = sorted({k[1] for t in T.values() for k in t})
        rA = {g: i for i, g in enumerate(self.gAs)}
        rB = {g: j for j, g in enumerate(self.gBs)}
        self.coef = {}
        for name, terms in T.items():
            deg = max(k for (_, _, k) in terms)
            C = np.zeros((deg + 1, len(self.gAs), len(self.gBs)))
            for (gA, gB, k), cf in terms.items():
                C[k, rA[gA], rB[gB]] += float(cf)
            self.coef[name] = C
        amax = int(max(max(g[0] + g[1] for g in self.gAs), max(g[0] + g[1] for g in self.gBs)))
        bmax = 1
        self.nmax = amax + bmax + dmax + 2 * r
        n0 = 2 * r
        self.n0 = n0
        self.nl = self.nmax - n0 + 1
        e = r * r - 1
        nl = self.nl
        import os
        fn = f'Gcache_r{r}_n0{n0}_nl{nl}.npy'
        if os.path.exists(fn):
            Gm = np.load(fn)
        elif floatG:
            Gm = G_float_table(r, n0, nl)
        else:
            Gm = np.zeros((nl * nl, nl * nl))
            for n1 in range(nl):
                for n2 in range(nl):
                    for n3 in range(nl):
                        for n4 in range(nl):
                            Gm[n1 * nl + n2, n3 * nl + n4] = float(core.G(n1 + n0, n2 + n0, n3 + n0, n4 + n0, e))
            np.save(fn, Gm)
        self.Gm = Gm
        # Phi[(p,q,b)] : (dmax+1) x nl matrix
        self.Phi = {}
        for p in range(amax + 1):
            for q in range(amax + 1 - p):
                for b in range(bmax + 1):
                    M = np.zeros((dmax + 1, nl))
                    for c in range(dmax + 1):
                        n = p + q + b + c + 2 * r
                        M[c, n - n0] = float(core.dirichlet(p, q, b + c, r))
                    self.Phi[(p, q, b)] = M

    def V(self, g, P):
        aa, ab, b1, b2 = g
        out = np.zeros((self.nl, self.nl))
        for i in range(aa + 1):
            for k in range(ab + 1):
                f1 = P @ self.Phi[(i, k, b1)]
                f2 = P @ self.Phi[(aa - i, ab - k, b2)]
                out += comb(aa, i) * comb(ab, k) * np.outer(f1, f2)
        return out.ravel()

    def Vall(self, g):
        """array (D, D, nl*nl): V(g, P) = sum_{c,c'} P_c P_c' Vall[c, c']."""
        aa, ab, b1, b2 = g
        D = self.dmax + 1
        out = np.zeros((D, D, self.nl, self.nl))
        for i in range(aa + 1):
            for k in range(ab + 1):
                A = self.Phi[(i, k, b1)]; B = self.Phi[(aa - i, ab - k, b2)]
                out += comb(aa, i) * comb(ab, k) * np.einsum('cn,dm->cdnm', A, B)
        return out.reshape(D, D, self.nl * self.nl)

    def tensors(self):
        """{name: array (K, D, D, D, D)}: c_name(u, P) = sum_k u^k T[k] contracted with P four times
        (indices: P-power in C1, C2 (group A) and C3, C4 (group B))."""
        VA = np.array([self.Vall(g) for g in self.gAs])      # (nA, D, D, nl2)
        VB = np.array([self.Vall(g) for g in self.gBs])      # (nB, D, D, nl2)
        W = np.einsum('ij,bcdj->bcdi', self.Gm, VB)          # (nB, D, D, nl2): G . VB
        out = {}
        for name, C in self.coef.items():
            Ts = []
            for k in range(C.shape[0]):
                Y = np.einsum('ab,bcdi->acdi', C[k], W)       # (nA, D, D, nl2)
                Ts.append(np.einsum('aefi,acdi->efcd', VA, Y))
            out[name] = np.array(Ts)
        return out

    def cs(self, P):
        P = np.asarray(P, float)
        P = np.concatenate([P, np.zeros(self.dmax + 1 - len(P))])
        VA = np.array([self.V(g, P) for g in self.gAs])
        VB = np.array([self.V(g, P) for g in self.gBs])
        M = VA @ self.Gm @ VB.T
        return {name: np.einsum('kij,ij->k', C, M) for name, C in self.coef.items()}


def lam0(v):
    return (1 + 4 * v + np.sqrt(1 + 8 * v)) / 8


def kappa2(c, c1, c2, v):
    return (np.sqrt(144 * v ** 2 * c2 ** 2 + 48 * lam0(v) * c * c1) - 12 * v * c2) / (16 * c1)


def best_uv(cs, u0=0.6, v0=0.4):
    c = cs['c'][0]
    p1 = np.poly1d(cs['c1'][::-1]); p2 = np.poly1d(cs['c2'][::-1])
    best = (-1, None)
    for u in np.linspace(u0 - 0.6, u0 + 0.6, 241):
        C1, C2 = p1(u), p2(u)
        if C1 <= 0:
            continue
        vs = np.linspace(0, 2, 801)
        k2 = kappa2(c, C1, C2, vs)
        i = int(np.nanargmax(k2))
        if k2[i] > best[0]:
            best = (k2[i], (u, vs[i]))
    # refine
    (u, v) = best[1]
    hu, hv = 0.005, 0.0025
    for _ in range(40):
        improved = False
        for du, dv in ((hu, 0), (-hu, 0), (0, hv), (0, -hv)):
            uu, vv = u + du, max(v + dv, 0)
            C1 = p1(uu)
            if C1 <= 0: continue
            k2 = kappa2(c, C1, p2(uu), vv)
            if k2 > best[0]:
                best = (k2, (uu, vv)); u, v = uu, vv; improved = True
        if not improved:
            hu /= 2; hv /= 2
    return np.sqrt(best[0]), best[1]
