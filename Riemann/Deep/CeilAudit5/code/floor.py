#!/usr/bin/env python3
# ceilaudit5 floor.py (N1): least eigenvalue of one off-line zero pair's rank-2 term 2Re[G(d+ig) conj G(-d+ig)], G(z) = int g e^{zx},
# on L^2[-a/2, a/2] (discretized), vs sinh(d a)/d; and the same on the cone g >= 0 (two bumps at +-a/2).
import math
import numpy as np
gam = 14.134725; out = []
for d in (0.05, 0.1, 0.2):
    for a in (5.0, 10.0, 20.0, 40.0, 80.0):
        N = 4000; x = (np.arange(N) + 0.5)*a/N - a/2; w = a/N
        u = np.exp((d + 1j*gam)*x); v = np.exp((-d + 1j*gam)*x)
        # q(g) = 2 Re[(u.g w)(conj(v.g w))] = g^* H g with H = w^2 (conj(u) v^T + conj(v) u^T)... use the 2x2 reduction on span{conj u, conj v}
        A = np.array([[np.vdot(u, u), np.vdot(u, v)], [np.vdot(v, u), np.vdot(v, v)]])*w   # Gram of u, v
        M = np.array([[0, 1], [1, 0]], dtype=complex)  # form in coordinates c: q = c^* (Gram M Gram) c over g = sum c_k conj(basis)
        # least eigenvalue of q relative to ||g||^2 on span: generalized eigenproblem Gram M Gram c = lam Gram c  ->  M Gram c = lam c
        lam = np.linalg.eigvals(M @ A).real.min()
        # cone g >= 0: two narrow bumps at +-a/2 (width a/50), phase cos(gam a) chosen worst by shifting a slightly
        out.append((d, a, lam, -(math.sinh(d*a)/d)))
        print('delta %.2f a %5.1f: least eigenvalue %.4e   -(sinh(delta a)/delta) %.4e   ratio %.4f' % (d, a, lam, -(math.sinh(d*a)/d), lam/(-(math.sinh(d*a)/d))))
