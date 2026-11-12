"""Mechanism: real amplifier times a phase, f(U) = g(U) exp(i mu U) (= mollifier shifted to M(s - i mu/L) up to orientation).
Hermitian matrices from ctest.mats (orientation L); real g of degree d; scan mu at fixed kappa. usage: mutest.py r kappa d mus"""
import sys, time
import os
HERE = os.path.dirname(os.path.abspath(__file__))
import numpy as np
from numpy.polynomial.legendre import leggauss, legvander
sys.path.insert(0, HERE)
import gcore as G
from ctest import mats, best
t0 = time.time()
r = int(sys.argv[1]); kap = float(sys.argv[2]); d = int(sys.argv[3]); mus = [float(x) * kap * np.pi for x in sys.argv[4].split(',')]
cplx = not (len(sys.argv) > 5 and sys.argv[5] == 'reK')
th = 0.25; K = d + 24; N = K + 2 * r + 2; D = 36
M0, MS, MSS, _ = mats(kap * np.pi, th, r, N, D, cplx)
B = G.Basis(th, N + D + 2, nq=N + D + 2 + 2 * r * r + 12)
Ir = np.linalg.matrix_power(B.I, 2 * r)[:N, :K]          # orthonormal f-coeffs (K) -> phi coeffs (N)
x, w = leggauss(K + 40); U = th * (x + 1) / 2; wq = w * th / 2
Vq = legvander(x, K - 1) * np.sqrt((2 * np.arange(K) + 1) / th)   # orthonormal Legendre on [0, th] at nodes
Vg = legvander(x, d) * np.sqrt((2 * np.arange(d + 1) + 1) / th)
for mu in mus:
    E = (Vq * wq[:, None]).T @ (Vg * np.exp(1j * mu * U)[:, None])   # K x (d+1): coefficients of exp(i mu U) q_n
    C = Ir @ E
    Qc, Rc = np.linalg.qr(np.vstack([C.real, C.imag]))
    Ri = np.linalg.inv(Rc)
    m = [Ri.T @ np.real(C.conj().T @ A @ C) @ Ri for A in (M0, MS, MSS)]
    m = [0.5 * (a + a.T) for a in m]
    lam, v = G.best_v(*m, lo=0.6, hi=2.2)
    print(('K   ' if cplx else 'ReK ') + 'r=%d kappa=%.4f deg g=%d mu=%+.4f (mu/c=%+.4f): h1 = %.7f (v=%.4f) [%.0fs]' % (r, kap, d, mu, mu / (kap * np.pi), lam / kap ** 2, v, time.time() - t0), flush=True)
