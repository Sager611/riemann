import sys, time, resource
import numpy as np
from numpy.polynomial import legendre as Lg
import gcore as G

def phi_of_P(Pmono, th, B):
    """orthonormal-Legendre coefficients (length N0) of phi = I^2 f, f(U) = P(1-U/th)."""
    X, w = B.X, B.w
    f = np.polyval(Pmono[::-1], 1 - X / th)
    a = B.V.T @ (w * f)          # coefficients of f
    return B.I @ (B.I @ a)

t0 = time.time()
th = 0.25
for kap, v, P, ref in [(3.18, 1.26, [1, -5.8, 6.4], 1.0021778)]:
    c = kap * np.pi
    N = 30; D = 36
    M0, MS, MSS, B = G.mats(c, th, N, D)
    # full-size phi (N0) -> need matrices on N0; recompute with N = N0 small poly: phi degree <= 4, so first N fine
    ph = phi_of_P(np.array(P, float), th, B)[:N]
    Q0 = ph @ M0 @ ph
    Q1 = ph @ (v * v * M0 - 2 * v * MS + MSS) @ ph
    print('kap=%.4f v=%.3f P=%s: h1 = %.9f (BM/cert ref %.7f)  Q0=%.6e Q1=%.6e' % (kap, v, P, Q0 / (kap ** 2 * Q1), ref, Q0, Q1))
    ev, vb = G.best_v(M0, MS, MSS)
    print('   pencil N=%d: sup_v lambda_max = %.6f at v=%.5f ; lambda/kappa^2 = %.6f' % (N, ev, vb, ev / kap ** 2))
print('time %.1fs peak %.0f MB' % (time.time() - t0, resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2 ** 20))
