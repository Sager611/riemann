"""float: best COMPLEX amplifier P of degree d for BM's criterion (general r, theta), Hermitian pencil
(orientation: Q = sum_ij Gam_ij conj(I^i phi)^T W_p (I^j phi), Gam Hermitian; P -> conj P gives the other orientation).
usage: copt.py r theta d kappa [scan]  -> logs/copt_r<r>_d<d>_k<kappa>.json ; 'scan' = find kappa_max(d) instead"""
import sys, json, time, resource
import os
HERE = os.path.dirname(os.path.abspath(__file__))
import numpy as np
from scipy.optimize import brentq
sys.path.insert(0, HERE)
import gcore as G
from ctest import mats, best, top

def Bas(th, r, N0):
    B = G.Basis(th, N0, nq=N0 + 2 * r * r + 12)
    return B

def reduced(kap, th, r, d, D=40):
    Nb = d + 2 * r + 6
    M0, MS, MSS, _ = mats(kap * np.pi, th, r, Nb, D, True)
    B = G.Basis(th, Nb + D + 2, nq=Nb + D + 2 + 2 * r * r + 12)
    Ir = np.linalg.matrix_power(B.I, 2 * r)[:Nb, :d + 1]
    Q, R = np.linalg.qr(Ir)
    return [Q.T @ A @ Q for A in (M0, MS, MSS)], Q, R, B

def main():
    t0 = time.time()
    r = int(sys.argv[1]); th = float(sys.argv[2]); d = int(sys.argv[3]); kap = float(sys.argv[4])
    if len(sys.argv) > 5 and sys.argv[5] == 'scan':
        f = lambda k: best(*reduced(k, th, r, d)[0])[0] / k ** 2 - 1
        k = brentq(f, 3.0, 3.42, xtol=1e-6)
        print('r=%d theta=%.4f deg=%d: complex kappa_max = %.6f [%.0fs]' % (r, th, d, k, time.time() - t0), flush=True)
        return
    (m0, ms, mss), Q, R, B = reduced(kap, th, r, d)
    lam, v = best(m0, ms, mss)
    M1 = v * v * m0 - 2 * v * ms + mss
    L = np.linalg.cholesky(M1); Li = np.linalg.inv(L)
    ev, U = np.linalg.eigh(Li @ m0 @ Li.conj().T)
    y = Li.conj().T @ U[:, -1]
    a = np.linalg.solve(R, y)
    An = a * B.norm[:d + 1]
    An = An / An[0]
    print('r=%d theta=%.4f deg=%d kappa=%.5f: float complex h1 = %.8f at v = %.6f [%.0fs]' % (r, th, d, kap, lam / kap ** 2, v, time.time() - t0))
    json.dump({'r': r, 'theta': sys.argv[2], 'd': d, 'kappa': sys.argv[4], 'v': '%.6f' % v, 'h1_float': lam / kap ** 2,
               'Are': ['%.15e' % x.real for x in An], 'Aim': ['%.15e' % x.imag for x in An]},
              open(HERE + '/logs/copt_r%d_d%d_k%s.json' % (r, d, sys.argv[4]), 'w'), indent=1)

if __name__ == '__main__':
    main()
