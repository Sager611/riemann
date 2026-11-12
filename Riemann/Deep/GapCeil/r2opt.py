"""float: best amplifier P of degree d for BM's criterion with general r at theta (phi-Galerkin on I^{2r} P_d).
Prints the unnormalised shifted-Legendre coefficients A_n of f(U) = P(1-U/theta) = sum A_n P_n(2U/theta-1) and v.
usage: r2opt.py r theta d kappa  -> logs/r2opt_r<r>_d<d>_k<kappa>.json"""
import sys, json
import numpy as np
import gcore as G
import gcore_r as GR
r = int(sys.argv[1]); th = float(sys.argv[2]); d = int(sys.argv[3]); kap = float(sys.argv[4])
Nb = d + 2 * r + 6
M0, MS, MSS, B = GR.mats_r(kap * np.pi, th, r, Nb, 40)
Ir = np.linalg.matrix_power(B.I, 2 * r)[:Nb, :d + 1]
Q, R = np.linalg.qr(Ir)
m0, ms, mss = (Q.T @ A @ Q for A in (M0, MS, MSS))
lam, v = G.best_v(m0, ms, mss, lo=0.8, hi=2.0)
ev, y = G.top(m0, ms, mss, v)
a = np.linalg.solve(R, y)                      # orthonormal coefficients of f
An = a * B.norm[:d + 1]                        # unnormalised P_n(2U/th-1) coefficients
An = An / An[0]
print('r=%d theta=%.4f deg=%d kappa=%.5f: float h1 = %.8f at v = %.6f' % (r, th, d, kap, lam / kap ** 2, v))
json.dump({'r': r, 'theta': sys.argv[2], 'd': d, 'kappa': sys.argv[4], 'v': '%.6f' % v, 'h1_float': lam / kap ** 2,
           'A': ['%.15e' % x for x in An]}, open('logs/r2opt_r%d_d%d_k%s.json' % (r, d, sys.argv[4]), 'w'), indent=1)
