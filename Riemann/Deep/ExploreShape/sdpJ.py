"""explore53: per-gap cone test for J-th order local Wirtinger inequalities (real symmetric case), any J >= 1.
Valid form on [0,g], g <= pi:  sum_{j<=J} a_j |f^{(j)}|^2 + sum_{p odd, q even, 1<=p,q<=J-1} phi_pq [f^{(p)} f^{(q)}]_0^g >= 0
for real f in H^J with f(0)=f(g)=0 (the bracket terms telescope over consecutive gaps; R-odd jet forms suffice by reflection).
Data: m_{2j} = central moments of nu times kappa^{2j}.  Minimise sum a_j m_{2j} with a_1 = 1.  RELAXATION (finite
polynomial basis, finite grid): value >= 0 proves no gain for these data; value < Wirtinger's m2 - 1 suggests gain.
LMI after dividing by g^{1-2J}, s = g^2 in [0, pi^2]: sum_j a_j s^{J-j} G_j + sum phi_pq s^{(2J-1-p-q)/2} B_pq >= 0.
usage: sdpJ.py mom.json J N [nsig]"""
import sys, json
import numpy as np
from numpy.polynomial import legendre as Lg
from numpy.polynomial import Polynomial as Pn
import cvxpy as cp


def basis(N, J):
    x, w = Lg.leggauss(N + J + 10)
    xi = (x + 1) / 2; w = w / 2
    # stable: Legendre series on [0,1] (variable y = 2xi-1), h_n = xi(1-xi) P_n(y) = (1-y^2)/4 P_n(y)
    bub = Lg.poly2leg([0.25, 0, -0.25])
    polys = []
    for n in range(N):
        c = np.zeros(n + 1); c[n] = 1
        polys.append(Lg.legmul(bub, c))
    y = 2 * xi - 1
    def der(cf, k):
        return Lg.legder(cf, k) * 2.0 ** k if k else cf
    G = []
    for j in range(J + 1):
        V = np.array([Lg.legval(y, der(p, j)) for p in polys])
        G.append((V * w) @ V.T)
    jets = {k: np.array([[Lg.legval(-1.0, der(p, k)), Lg.legval(1.0, der(p, k))] for p in polys]) for k in range(1, J)}
    B = {}
    for p_ in range(1, J, 2):
        for q in range(2, J, 2):
            M = np.outer(jets[p_][:, 1], jets[q][:, 1]) - np.outer(jets[p_][:, 0], jets[q][:, 0])
            B[(p_, q)] = 0.5 * (M + M.T)
    S = sum(G)
    ev, U = np.linalg.eigh(S)
    ev = np.maximum(ev, ev.max() * 1e-15)   # J >= 6: clip roundoff-negative eigenvalues (basis change only)
    T = U / np.sqrt(ev)
    return [T.T @ g @ T for g in G], {k: T.T @ b @ T for k, b in B.items()}


def main():
    mom = json.load(open(sys.argv[1])); J = int(sys.argv[2]); N = int(sys.argv[3])
    nsig = int(sys.argv[4]) if len(sys.argv) > 4 else 50
    kap = mom['kappa']; cen = mom['central_about_mean']
    m = [cen[2 * j] * kap ** (2 * j) / cen[0] for j in range(J + 1)]
    G, B = basis(N, J)
    a = cp.Variable(J + 1); ph = {k: cp.Variable() for k in B}
    cons = [a[1] == 1, a[J] >= 0] + [cp.abs(v) <= 1e3 for v in ph.values()] + [cp.abs(a) <= 1e3]
    sg = np.concatenate([[0.0], np.pi ** 2 * np.geomspace(1e-8, 1, nsig)])
    for s in sg:
        M = sum(a[j] * s ** (J - j) * G[j] for j in range(J + 1))
        for (p_, q), Bm in B.items():
            M = M + ph[(p_, q)] * s ** ((2 * J - 1 - p_ - q) // 2) * Bm
        cons.append(0.5 * (M + M.T) >> 0)
    pr = cp.Problem(cp.Minimize(sum(a[j] * m[j] for j in range(J + 1))), cons)
    pr.solve(solver=cp.CLARABEL)
    print('kappa=%.4f J=%d N=%d nsig=%d: Wirtinger m2-1 = %.7f ; J-value = %.7f (%s)' % (kap, J, N, nsig, m[1] - 1, pr.value, pr.status))
    print('  a =', ' '.join('%.6g' % x for x in a.value), ' phi =', {k: round(float(v.value), 6) for k, v in ph.items()})


if __name__ == '__main__':
    main()
