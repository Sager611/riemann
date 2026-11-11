"""Well-conditioned float evaluator: P in shifted-Legendre basis (w on [0,1]); functions of each R_k in
shifted-Legendre basis (R on [0,1]); all basis changes done exactly (fmpq) before converting to float."""
import sys
sys.path.insert(0, '.')
import numpy as np
from math import comb
from flint import fmpq, fmpq_mat
import core


def shleg_mono(n):
    """exact monomial coefficients of the shifted Legendre polynomial on [0,1] of degree n."""
    return [fmpq((-1) ** (n + k) * comb(n, k) * comb(n + k, k)) for k in range(n + 1)]


def basis_mat(J):
    """S[k, m] = coeff of x^k in psi_m(x), m,k < J (exact fmpq_mat)."""
    S = fmpq_mat(J, J)
    for m in range(J):
        for k, cf in enumerate(shleg_mono(m)):
            S[k, m] = cf
    return S


def tofloat(M):
    return np.array([[float(M[i, j]) for j in range(M.ncols())] for i in range(M.nrows())])


class FastL:
    def __init__(self, theta, r, dmax):
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
        D = dmax + 1
        n0 = 2 * r
        nmax = amax + bmax + dmax + 2 * r
        nl = nmax - n0 + 1
        J = nmax + 1
        self.D, self.J, self.n0, self.nl = D, J, n0, nl
        # P basis: P_mono = BP pL
        BP = basis_mat(D)
        self.BP = BP
        BPt = BP.transpose()
        # R basis: mono (powers n0..nmax) -> Legendre coefficients (0..nmax): f_L = Kemb f_mono
        SR = basis_mat(J)
        SRinv = SR.inv()
        Kemb = fmpq_mat(J, nl, [SRinv[i, n0 + j] for i in range(J) for j in range(nl)])
        Kt = fmpq_mat(J, nl, [SR[n0 + j, i] for i in range(J) for j in range(nl)])   # K^T, K = SR[n0:, :]
        self.Phi = {}
        for p in range(amax + 1):
            for q in range(amax + 1 - p):
                for b in range(bmax + 1):
                    M = fmpq_mat(D, nl)
                    for c in range(D):
                        n = p + q + b + c + 2 * r
                        M[c, n - n0] = core.dirichlet(p, q, b + c, r)
                    self.Phi[(p, q, b)] = tofloat(BPt * M * Kemb.transpose())      # (D, J)
        # exact G table, then G_L = G x1 K^T x2 K^T x3 K^T x4 K^T
        from flint import arb, arb_mat, ctx
        oldprec = ctx.prec
        ctx.prec = 384                                   # transform in ball arithmetic (G > 0 exact; K has entries ~1e12)
        Hx = core.G_exact_table(r, n0, nl, return_H=True)
        Ha = {k: arb(val.p) / arb(val.q) for k, val in Hx.items()}
        del Hx
        arr = np.empty((nl, nl, nl, nl), dtype=object)
        for (a1, a2, a3, a4), val in Ha.items():
            arr[a1 - n0, a2 - n0, a3 - n0, a4 - n0] = val + Ha[(a3, a4, a1, a2)]
        del Ha
        Kta = arb_mat(J, nl, [arb(Kt[i, j].p) / arb(Kt[i, j].q) for i in range(J) for j in range(nl)])
        for ax in range(4):                              # mode product on axis 0, then rotate axes
            s = arr.shape
            M = arb_mat(s[0], int(np.prod(s[1:])), list(arr.reshape(s[0], -1).ravel()))
            del arr
            R = Kta * M                                  # (J, prod rest)
            del M
            out = np.empty((J,) + s[1:], dtype=object)
            flat = out.reshape(J, -1)
            nc = flat.shape[1]
            for i in range(J):
                for j in range(nc):
                    flat[i, j] = R[i, j]
            del R, flat
            arr = np.moveaxis(out, 0, -1)                # processed axis goes last; after 4 steps order is restored
            if ax == 3:
                break
        rad = max(float(x.rad()) for x in arr.ravel())
        GL = np.vectorize(lambda x: float(x.mid()), otypes=[float])(arr)
        self.Grad = rad / max(np.abs(GL).max(), 1e-300)  # relative radius of the ball transform
        del arr, out
        ctx.prec = oldprec
        self.Gm = GL.reshape(J * J, J * J)

    def Vall(self, g):
        aa, ab, b1, b2 = g
        D, J = self.D, self.J
        out = np.zeros((D, D, J, J))
        for i in range(aa + 1):
            for k in range(ab + 1):
                A = self.Phi[(i, k, b1)]; B = self.Phi[(aa - i, ab - k, b2)]
                out += comb(aa, i) * comb(ab, k) * np.einsum('cn,dm->cdnm', A, B)
        return out.reshape(D, D, J * J)

    def tensors(self):
        """{name: (K, D, D, D, D)} in the Legendre basis of P."""
        D, J2 = self.D, self.J * self.J
        nA, nB = len(self.gAs), len(self.gBs)
        VA = np.empty((nA, D * D * J2))
        for i, g in enumerate(self.gAs):
            VA[i] = self.Vall(g).ravel()
        acc = {name: np.zeros((C.shape[0], D * D, D * D)) for name, C in self.coef.items()}
        for j, g in enumerate(self.gBs):
            Wb = self.Vall(g).reshape(D * D, J2) @ self.Gm.T          # (D^2, J2)
            for name, C in self.coef.items():
                for k in range(C.shape[0]):
                    col = C[k][:, j]
                    if not col.any():
                        continue
                    U = (col @ VA).reshape(D * D, J2)                   # sum_a C[a,j] VA[a]
                    acc[name][k] += U @ Wb.T
        return {name: a.reshape(a.shape[0], D, D, D, D) for name, a in acc.items()}

    def mono_to_leg(self, Pmono):
        """exact: list of fmpq monomial coefficients -> list of fmpq Legendre coefficients."""
        D = self.D
        v = fmpq_mat(D, 1, list(Pmono) + [fmpq(0)] * (D - len(Pmono)))
        x = self.BP.solve(v)
        return [x[i, 0] for i in range(D)]
