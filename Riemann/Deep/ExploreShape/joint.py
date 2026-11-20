"""explore53: joint optimisation of the J-th order local criterion and BM's amplifier (midpoint mollifier, r, theta, deg d).
Central forms C_k = <(S - v0)^k> kappa^k (k <= 2J) in the amplifier basis of GapBar's optmid; twist v = v0 + dv.
Loop: data m_2j(y) -> SDP (sdpJ relaxation, N-poly basis) -> form a, phi -> amplifier y = top eigvec of (C0, R),
R = sum_{j>=1} a_j C_2j(v) (h_J = -a0 lambda_max(C0, R) > 1 means Lambda > kappa) -> repeat.
usage: joint.py r theta d kappa J [iters] [N]"""
import sys, os, json
from math import comb, factorial
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import moments as Mo
import gcore as G
import cvxpy as cp
import sdpJ


def central_forms(r, th, d, kap, K, v0, D=40):
    Nb = d + 2 * r + 6
    H = [h.real for h in Mo.arrays_K(kap * np.pi, r, K, D)]
    N0 = Nb + D + 2
    B = G.Basis(th, N0, nq=N0 + 2 * r * r + 12 + K)
    q0 = r * r - 1
    Wc = [(B.V * (B.w * B.X ** q0 * (B.X - v0) ** p)[:, None]).T @ B.V for p in range(K + 1)]
    Ipow = [np.eye(N0)]
    for k in range(D):
        Ipow.append(B.I @ Ipow[-1])
    Ip = [P[:, :Nb] for P in Ipow]
    s0 = 2 * r - 1
    fac = [float(factorial(k)) for k in range(D + s0 + 2)]
    def Gm(p, h):
        out = np.zeros((Nb, Nb)); L = [Wc[p] @ Ip[k] for k in range(D)]
        for i in range(D):
            for j in range(D):
                if h[i, j] != 0:
                    out += fac[i + s0] * fac[j + s0] * h[i, j] * (Ip[i].T @ L[j])
        return 0.5 * (out + out.T)
    Ir = np.linalg.matrix_power(B.I, 2 * r)[:Nb, :d + 1]
    Q, _ = np.linalg.qr(Ir)
    C = []
    for k in range(K + 1):
        A = sum(comb(k, i) * Gm(i, H[k - i]) for i in range(k + 1))
        C.append((Q.T @ A @ Q) * kap ** k)
    return C


def shift(C, dv, kap):
    """forms about v = v0 + dv: (S-v)^k = sum_i binom(k,i)(S-v0)^i (-dv)^{k-i}; scaled units: dv*kap"""
    t = -dv * kap
    return [sum(comb(k, i) * C[i] * t ** (k - i) for i in range(k + 1)) for k in range(len(C))]


def lmax(A, Bm):
    L = np.linalg.cholesky(Bm); Li = np.linalg.inv(L)
    ev, U = np.linalg.eigh(Li @ A @ Li.T)
    return ev[-1], Li.T @ U[:, -1]


def sdp(m, J, N, nsig=30):
    G_, Bd = sdpJ.basis(N, J)
    a = cp.Variable(J + 1); ph = {k: cp.Variable() for k in Bd}
    cons = [a[1] == 1, a[J] >= 0] + [cp.abs(v) <= 1e3 for v in ph.values()] + [cp.abs(a) <= 1e3]
    for s in np.concatenate([[0.0], np.pi ** 2 * np.geomspace(1e-8, 1, nsig)]):
        M = sum(a[j] * s ** (J - j) * G_[j] for j in range(J + 1))
        for (p_, q), Bm in Bd.items():
            M = M + ph[(p_, q)] * s ** ((2 * J - 1 - p_ - q) // 2) * Bm
        cons.append(0.5 * (M + M.T) >> 0)
    pr = cp.Problem(cp.Minimize(sum(a[j] * m[j] for j in range(J + 1))), cons)
    pr.solve(solver=cp.CLARABEL)
    return pr.value, np.array(a.value), {k: float(v.value) for k, v in ph.items()}


def run(r, th, d, kap, J, iters=4, N=14, v0=1.3617, verbose=True):
    C0 = central_forms(r, th, d, kap, 2 * J, v0)
    from scipy.optimize import minimize_scalar
    # Wirtinger start
    fw = lambda dv: -lmax(*[shift(C0, dv, kap)[i] for i in (0, 2)])[0]
    rw = minimize_scalar(fw, bounds=(-0.01, 0.01), method='bounded', options={'xatol': 1e-9})
    dv = rw.x; h1 = -rw.fun
    _, y = lmax(*[shift(C0, dv, kap)[i] for i in (0, 2)])
    hist = [('W', h1, dv)]
    a = None; ph = None; hJ = None
    for it in range(iters):
        C = shift(C0, dv, kap)
        n0 = y @ C[0] @ y
        m = [y @ C[2 * j] @ y / n0 for j in range(J + 1)]
        val, a, ph = sdp(m, J, N)
        def fa(dvv):
            Cs = shift(C0, dvv, kap)
            R = sum(a[j] * Cs[2 * j] for j in range(1, J + 1))
            return -(-a[0]) * lmax(Cs[0], R)[0]
        rr = minimize_scalar(fa, bounds=(dv - 0.003, dv + 0.003), method='bounded', options={'xatol': 1e-9})
        dv = rr.x; hJ = -rr.fun
        Cs = shift(C0, dv, kap)
        _, y = lmax(Cs[0], sum(a[j] * Cs[2 * j] for j in range(1, J + 1)))
        hist.append((it, val, hJ, dv))
        if verbose:
            print('  it %d: SDP value %.6g -> h_J = %.8f (v = %.6f)' % (it, val, hJ, v0 + dv), flush=True)
    return h1, hJ, a, ph, v0 + dv, hist


if __name__ == '__main__':
    r = int(sys.argv[1]); th = float(sys.argv[2]); d = int(sys.argv[3]); kap = float(sys.argv[4]); J = int(sys.argv[5])
    iters = int(sys.argv[6]) if len(sys.argv) > 6 else 4; N = int(sys.argv[7]) if len(sys.argv) > 7 else 14
    h1, hJ, a, ph, v, hist = run(r, th, d, kap, J, iters, N)
    print('r=%d d=%d kappa=%.4f J=%d: h1 (Wirtinger) = %.8f ; h_J = %.8f  a = %s phi = %s' % (r, d, kap, J, h1, hJ, ' '.join('%.6g' % x for x in a), ph))
    json.dump({'r': r, 'theta': th, 'd': d, 'kappa': kap, 'J': J, 'h1': h1, 'hJ': hJ, 'a': list(a), 'phi': {str(k): v for k, v in ph.items()},
               'v': v, 'N': N}, open(os.path.join(HERE, 'logs', 'joint_r%d_d%d_k%s_J%d.json' % (r, d, sys.argv[4], J)), 'w'), indent=1)
