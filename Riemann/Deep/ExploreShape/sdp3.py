"""explore53: per-gap cone test for third-order local Wirtinger inequalities (J = 3), real symmetric case.
Valid form on a gap [0,g], g <= pi (units: gaps <= kappa*pi in u = tL become <= pi after omega -> kappa*omega):
  a0|f|^2 + a1|f'|^2 + a2|f''|^2 + a3|f'''|^2 + phi [f'' f']_0^g >= 0  for all real f in H^3, f(0)=f(g)=0.
(The phi term telescopes over consecutive gaps: f', f'' are continuous at the zeros.)  Even symbol about the twist v =
mean of nu; nu's odd central moments are ~1e-3 (moments.py), handled separately.
Contradiction (Lambda > kappa) if a0 m0 + a1 m2 + a2 m4 + a3 m6 < 0, m = central moments of nu scaled by kappa^(2j).
SDP: minimise that subject to a1 = 1 (normalisation) and the LMI on a sigma-grid (sigma = pi/g >= 1), polynomial basis.
This is a RELAXATION (finite basis, finite grid): a value >= 0 proves no J=3 gain for these data; < 0 needs verification.
usage: sdp3.py mom.json N [nsig] [sigmax]"""
import sys, json
import numpy as np
from numpy.polynomial import legendre as Lg
import cvxpy as cp


def basis_mats(N, nq=None):
    nq = nq or N + 12
    x, w = Lg.leggauss(nq)
    xi = (x + 1) / 2; w = w / 2
    # h_n = xi(1-xi) P_n(2xi-1): derivatives via polynomial arithmetic in xi
    from numpy.polynomial import Polynomial as Pn
    polys = []
    for n in range(N):
        c = np.zeros(n + 1); c[n] = 1
        Pleg = Lg.Legendre(c, domain=[0, 1])        # P_n(2xi-1)
        p = Pleg.convert(kind=Pn) * Pn([0, 1, -1])
        polys.append(p)
    G = []
    for j in range(4):
        V = np.array([(p.deriv(j) if j else p)(xi) for p in polys])
        G.append((V * w) @ V.T)
    d1 = np.array([[p.deriv(1)(0.0), p.deriv(1)(1.0)] for p in polys])
    d2 = np.array([[p.deriv(2)(0.0), p.deriv(2)(1.0)] for p in polys])
    # [f'' f']_0^1 = f''(1) f'(1) - f''(0) f'(0): symmetric matrix
    Bm = np.outer(d2[:, 1], d1[:, 1]) - np.outer(d2[:, 0], d1[:, 0])
    Bm = 0.5 * (Bm + Bm.T)
    # orthonormalise w.r.t. G0 + G3 for conditioning
    S = G[0] + G[1] + G[2] + G[3]
    ev, U = np.linalg.eigh(S)
    T = U / np.sqrt(ev)
    return [T.T @ g @ T for g in G], T.T @ Bm @ T


def main():
    mom = json.load(open(sys.argv[1]))
    N = int(sys.argv[2]); nsig = int(sys.argv[3]) if len(sys.argv) > 3 else 60
    smax = float(sys.argv[4]) if len(sys.argv) > 4 else 60.0
    kap = mom['kappa']; cen = mom['central_about_mean']
    m = [cen[2 * j] * kap ** (2 * j) / cen[0] for j in range(4)]
    G, Bm = basis_mats(N)
    a = cp.Variable(4); phi = cp.Variable()
    cons = [a[1] == 1, a[3] >= 0]
    # s = g^2 in (0, pi^2]; the LMI is divided by g^-6: a0 s^3 G0 + a1 s^2 G1 + s (a2 G2 + phi B) + a3 G3 >= 0;
    # s -> 0 (tiny gaps) is included down to 1e-12 and at s = 0 exactly.
    sg = np.concatenate([[0.0], np.pi ** 2 * np.geomspace(10.0 ** -smax, 1, nsig)])
    for s in sg:
        M = a[0] * s ** 3 * G[0] + a[1] * s ** 2 * G[1] + s * (a[2] * G[2] + phi * Bm) + a[3] * G[3]
        cons.append(0.5 * (M + M.T) >> 0)
    obj = a[0] * m[0] + a[1] * m[1] + a[2] * m[2] + a[3] * m[3]
    pr = cp.Problem(cp.Minimize(obj), cons + [cp.abs(a[2]) <= 1e3, phi <= 1e3, phi >= -1e3])
    pr.solve(solver=cp.CLARABEL)
    print('kappa=%.4f m=(%.6f, %.6f, %.6f, %.6f)  Wirtinger value m2-m0 = %.6f' % (kap, m[0], m[1], m[2], m[3], m[1] - m[0]))
    print('N=%d nsig=%d smax=%g status=%s  J=3 min value = %.6g' % (N, nsig, smax, pr.status, pr.value))
    print('  a =', np.round(a.value, 8), 'phi =', phi.value)


if __name__ == '__main__':
    main()
