"""explore53: validate a per-gap form sum a_j |f^{(j)}|^2 + sum phi_pq [f^(p) f^(q)]_0^g on a large stable polynomial basis
(float64, Legendre series), for all gaps g in (0, pi] (fine s = g^2 grid).  Reports min over s of the smallest eigenvalue of
the form relative to the H^J Gram matrix (sum_j G_j at the same scale).  Negative => invalid (counterexample function).
usage: from validate import minform; minform(a, phi, J, N)"""
import numpy as np
from numpy.polynomial import legendre as Lg


def basis(N, J):
    x, w = Lg.leggauss(N + J + 10)
    xi = (x + 1) / 2; w = w / 2
    bub = Lg.poly2leg([0.25, 0, -0.25])
    polys = []
    for n in range(N):
        c = np.zeros(n + 1); c[n] = 1
        polys.append(Lg.legmul(bub, c))
    y = 2 * xi - 1
    der = lambda cf, k: Lg.legder(cf, k) * 2.0 ** k if k else cf
    G = [(np.array([Lg.legval(y, der(p, j)) for p in polys]) * w) @ np.array([Lg.legval(y, der(p, j)) for p in polys]).T
         for j in range(J + 1)]
    jets = {k: np.array([[Lg.legval(-1.0, der(p, k)), Lg.legval(1.0, der(p, k))] for p in polys]) for k in range(1, J)}
    B = {}
    for p_ in range(1, J, 2):
        for q in range(2, J, 2):
            M = np.outer(jets[p_][:, 1], jets[q][:, 1]) - np.outer(jets[p_][:, 0], jets[q][:, 0])
            B[(p_, q)] = 0.5 * (M + M.T)
    return G, B


def minform(a, phi, J, N=40, ns=600, smin=1e-8):
    G, B = basis(N, J)
    worst = (np.inf, None)
    for s in np.concatenate([np.pi ** 2 * np.geomspace(smin, 1, ns), [np.pi ** 2]]):
        g = np.sqrt(s)
        # form on [0,g] for f(x) = h(x/g): |f^(j)|^2 = g^{1-2j}|h^(j)|^2, [f^(p) f^(q)] = g^{-p-q}[h^(p) h^(q)]
        M = sum(a[j] * g ** (1 - 2 * j) * G[j] for j in range(J + 1))
        for (p_, q), Bm in B.items():
            M = M + phi[(p_, q)] * g ** (-p_ - q) * Bm
        Nrm = sum(g ** (1 - 2 * j) * G[j] for j in range(J + 1))
        L = np.linalg.cholesky(Nrm)
        Li = np.linalg.inv(L)
        ev = np.linalg.eigvalsh(Li @ M @ Li.T)[0]
        if ev < worst[0]:
            worst = (ev, g)
    return worst
