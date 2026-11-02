"""Adversary: hardest unit polygonal arcs for a convex K = {n_i.x <= h_i} of the
paper's type (h0 = h1 = 0, sides j >= 2 in the negative cone of sides 0, 1).

rho(Z) = min over orientations M in O(2) of max_j Psi_j(M)/h_j, where
Psi_j = a_j0 h_Z(M n0) + a_j1 h_Z(M n1) + h_Z(M n_j)   (translation criterion).
Z fits in K iff rho(Z) <= 1; K scaled by s covers Z iff rho(Z) <= s.
"""
import numpy as np
from fractions import Fraction


def unit(t):
    t = float(Fraction(t))
    return np.array([(1 - t * t) / (1 + t * t), 2 * t / (1 + t * t)])


class Shape:
    def __init__(self, T, H, ngrid=2880):
        self.N = np.array([unit(t) for t in T])
        self.H = np.array([float(Fraction(h)) for h in H])
        k = len(T)
        n0, n1 = self.N[0], self.N[1]
        A = np.column_stack([n0, n1])
        rows = []
        for j in range(2, k):
            c = np.linalg.solve(A, -self.N[j])
            assert c.min() > 0, "side not in negative cone"
            rows.append(c)
        self.a = np.array(rows)  # (k-2, 2)
        self.k = k
        self.set_grid(ngrid)

    def set_grid(self, ngrid):
        th = np.linspace(0, 2 * np.pi, ngrid, endpoint=False)
        self.th = th
        c, s = np.cos(th), np.sin(th)
        # D[g, i, :] = R_theta n_i
        nx, ny = self.N[:, 0], self.N[:, 1]
        self.D = np.stack([c[:, None] * nx - s[:, None] * ny, s[:, None] * nx + c[:, None] * ny], axis=-1)
        self.D2T = np.ascontiguousarray(self.D.reshape(-1, 2).T)

    def ratios(self, Z, D=None):
        """Z: (p, 2). Returns (2G,) array of max_j Psi_j/h_j over grid orientations, for Z and mirror(Z)."""
        if D is None:
            D2T = self.D2T
            G = self.D.shape[0]
        else:
            G = D.shape[0]
            D2T = np.ascontiguousarray(D.reshape(-1, 2).T)
        out = []
        for ZZ in (Z, Z * np.array([-1.0, 1.0])):
            h = (ZZ @ D2T).max(axis=0).reshape(G, self.k)
            psi = h[:, [0]] * self.a[:, 0][None, :] + h[:, [1]] * self.a[:, 1][None, :] + h[:, 2:]
            out.append((psi / self.H[None, 2:]).max(axis=1))
        return np.concatenate(out)

    def ratio_at(self, ZZ, th):
        c, s = np.cos(th), np.sin(th)
        N = self.N
        D = np.empty((self.k, 2))
        D[:, 0] = c * N[:, 0] - s * N[:, 1]
        D[:, 1] = s * N[:, 0] + c * N[:, 1]
        h = (ZZ @ D.T).max(axis=0)
        psi = self.a[:, 0] * h[0] + self.a[:, 1] * h[1] + h[2:]
        return float((psi / self.H[2:]).max())

    def rho(self, Z, refine=True):
        r = self.ratios(Z)
        g = int(r.argmin())
        best = float(r[g])
        if not refine:
            return best
        G = len(self.th)
        ZZ = Z * np.array([-1.0, 1.0]) if g >= G else Z
        t0 = self.th[g % G]
        step = 2 * np.pi / G
        a, b = t0 - step, t0 + step
        gr = 0.6180339887498949
        x1, x2 = b - gr * (b - a), a + gr * (b - a)
        f1, f2 = self.ratio_at(ZZ, x1), self.ratio_at(ZZ, x2)
        for _ in range(40):
            if f1 < f2:
                b, x2, f2 = x2, x1, f1
                x1 = b - gr * (b - a)
                f1 = self.ratio_at(ZZ, x1)
            else:
                a, x1, f1 = x1, x2, f2
                x2 = a + gr * (b - a)
                f2 = self.ratio_at(ZZ, x2)
        return min(best, f1, f2)


def polyline(params, m):
    """params: m-1 turning directions (phi_1 = 0) and m length logits."""
    phi = np.concatenate([[0.0], params[: m - 1]])
    w = np.exp(params[m - 1:] - params[m - 1:].max())
    lens = w / w.sum()
    steps = np.stack([np.cos(phi), np.sin(phi)], axis=1) * lens[:, None]
    return np.vstack([[0.0, 0.0], np.cumsum(steps, axis=0)])


def hardest(shape, m=6, starts=40, iters=4000, seed=0, init=None):
    from scipy.optimize import minimize  # pyright: ignore[reportMissingImports]
    rng = np.random.default_rng(seed)
    best = (-1, None)
    for s in range(starts):
        if init is not None and s == 0:
            x0 = init
        else:
            x0 = np.concatenate([rng.uniform(-np.pi, np.pi, m - 1), rng.normal(0, 0.5, m)])
        f = lambda x: -shape.rho(polyline(x, m), refine=False)
        res = minimize(f, x0, method="Nelder-Mead", options=dict(maxiter=iters, xatol=1e-7, fatol=1e-9, adaptive=True))
        val = shape.rho(polyline(res.x, m), refine=True)
        if val > best[0]:
            best = (val, res.x)
    return best

