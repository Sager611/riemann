# B1 float test: PR 137/141 anchor induction with a Gauss-selection cost 4 - rho per unit slot width (instead of 4).
# Bookkeeping (High/Algebra.lean, ChainData via B2Below/chain.py): r = HR, s = H S_R + u/k, S_r = s + u, payment (K-F+(5-rho)u)+(H+u) = F,
# H = 1 - (6-rho)u. Then z_r = Hz + u(1/k + 1) - (1-H)/2 = Hz + u(1/k - 2 + rho/2); thr (gate 1, G(k) = H F0) does not involve u.
# Continuum: r z' = z - T_rho, T_rho = (1/thr - 2 + rho/2)/(6 - rho), thr = (5r/4 - 5/6)/(r - 1/2 - z); start z = 0 where T_rho = 0.
# mode 'udiv' reproduces B2Below's perturbation (keeps 1/k - 2, divides by 6 - rho). Count window + return: B2Below's ode.py pipeline.
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'B2Below'))
import numpy as np
from system import Broot, Rd_native
from infimum import golden
def make(rho, mode):
    def thr(r, z): return (1.25 * r - 5 / 6) / (r - 0.5 - z)
    def T(r, z):
        if mode == 'udiv': return (1 / thr(r, z) - 2) / (6 - rho)
        return (1 / thr(r, z) - 2 + rho / 2) / (6 - rho)
    # start: T(r0, 0) = 0  <=>  thr(r0, 0) = 1/(2 - c), c = rho/2 (correct) or 0 (udiv)
    c = rho / 2 if mode != 'udiv' else 0.0
    k0 = 1 / (2 - c); r0 = (5 / 6 - 0.5 * k0) / (1.25 - k0)
    f = lambda r, z: (z - T(r, z)) / r
    n = 120000; r1 = 0.725; h = (r1 - r0) / n; rs = [r0]; zs = [0.0]; r, z = r0, 0.0
    for _ in range(n):
        k1 = f(r, z); k2 = f(r + h / 2, z + h * k1 / 2); k3 = f(r + h / 2, z + h * k2 / 2); k4 = f(r + h, z + h * k3)
        z += h * (k1 + 2 * k2 + 2 * k3 + k4) / 6; r += h; rs.append(r); zs.append(z)
    rs = np.array(rs[::-1]); zs = np.array(zs[::-1])
    return (lambda w: np.interp(w, rs, zs)), r0
def contact(G, hi=0.7775):
    def Rc(d, x):
        V = 5 / 6 - d; F = lambda c: (1 - d * (c + 2 * x * ((1 - c) / 2 + G(c)))) - (1 - d + V * (1.5 * c - 1))
        a, b = 0.735, hi
        for _ in range(100):
            m = (a + b) / 2
            if (F(m) > 0) == (F(a) > 0): a = m
            else: b = m
        c = (a + b) / 2
        return min(1 - d + V * (1.5 * c - 1), Rd_native(d, x)), c
    def Bb(b):
        neg = lambda d: -Broot(b, d, 0.5, Rc(d, 0.5)[0])
        d = golden(neg, 0.38, 0.47, 1e-10)
        return -neg(d), d
    bo = golden(lambda b: Bb(b)[0], 0.09, 0.13, 1e-11)
    B, d = Bb(bo); R, c = Rc(d, 0.5)
    return B, bo, d, c, R
if __name__ == '__main__':
    TH = 0.8733458561067185768
    def rmax(p): return 6 * (p - 6) * (2 * p - 3) / (p * p - 193 * p + 198)
    print("rho_max(p): p=2 %.5f p=2.5 %.5f p=3 %.5f p=4.5 %.5f p=6*0.7498 %.5f" % (rmax(2), rmax(2.5), rmax(3), rmax(4.5), rmax(6 * 0.7498)))
    for mode in ('correct', 'udiv'):
        for rho in (0.0, 0.05, rmax(4.5), rmax(2), rmax(3), 1.0):
            G, r0 = make(rho, mode)
            B, b, d, c, R = contact(G)
            print("%-7s rho=%.5f start r0=%.6f z(0.7523)=%.4e  B*=%.13f  B*-theta_inf=%+.3e  b=%.6f delta=%.6f cross=%.6f R=%.8f" % (
                mode, rho, r0, float(G(0.7523)), B, B - TH, b, d, c, R)); sys.stdout.flush()
