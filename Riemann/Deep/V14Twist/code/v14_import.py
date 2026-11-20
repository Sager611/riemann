# V14Twist: import of Liu's V14 unit twist into PR 137/141's marked-anchor chain, continuum limit (copied pipeline: CeilAudit2/b1_twist.py, B2Below/system.py, infimum.py).
# PR 137 Poisson step (Select/Scale.lean selection_exponents; Assembly.lean RestoreChildAt): payment K'-N'-V+Hc' = N'+V (K' cancels), child excess Hc'x.
# Twist a part tau of the selected width u (V14 l.9096-9103): not extended (K' = 2F-1+4(u-tau)), removed from the column (N' = N-u), Cauchy over the free
# twisted tuple loses Z^tau. Hence Hc' = 1-6u+4tau, excess e = u-tau (mode 'free'); V14's own count R_loc+2s gives e = u-2tau (mode 'full').
# Stage (High/Algebra.lean with the excess u -> e, selected width u in S_r = s+u): s = H S_R + e/k, S_r = s + u, z_r = S_r - (1-HR)/2.
#   free:  z_r = Hz + e(1/k-2),          1-H = 6e+2tau  ->  T = (1/thr-2)/(6+2 sig),          sig = tau/e
#   full:  z_r = Hz + e(1/k-2) - 2tau,   1-H = 6e+8tau  ->  T = (1/thr-2-2 sig)/(6+8 sig)
#   fixed (hypothetical, no count loss; not available: PR 137's selected tuple is free): e = u, T = (1/thr-2+2rp)/(6-4rp), rp = tau/u
#   ca2 (CeilAudit2 model: extension cost 4-rho, payment unchanged): T = (1/thr-2+rho/2)/(6-rho)
# gate 1 (G k = H F0) does not involve the excess, so thr = (5r/4-5/6)/(r-1/2-z) in every mode.
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np
from system import Broot, Rd_native
from infimum import golden
def make(mode, a):
    def thr(r, z): return (1.25 * r - 5 / 6) / (r - 0.5 - z)
    if mode == 'ca2':   T = lambda r, z: (1 / thr(r, z) - 2 + a / 2) / (6 - a); c = a / 2
    elif mode == 'free':  T = lambda r, z: (1 / thr(r, z) - 2) / (6 + 2 * a); c = 0.0
    elif mode == 'full':  T = lambda r, z: (1 / thr(r, z) - 2 - 2 * a) / (6 + 8 * a); c = -2 * a
    elif mode == 'fixed': T = lambda r, z: (1 / thr(r, z) - 2 + 2 * a) / (6 - 4 * a); c = 2 * a
    k0 = 1 / (2 - c); r0 = (5 / 6 - 0.5 * k0) / (1.25 - k0)   # T(r0, 0) = 0
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
if __name__ == '__main__' and len(sys.argv) == 1:
    TH = 0.8733458561067185768
    rho = 24 / 289; rp = rho / 3; sig = rp / (1 - rp)
    runs = [('ca2', 0.0), ('ca2', rho), ('free', sig), ('free', 3 * sig), ('free', 0.5), ('full', sig), ('fixed', rp), ('fixed', rho)]
    print("V14 rho_max(9/2) = 24/289; twisted fraction rp = rho/3 = %.6f; sig = tau/e = %.6f" % (rp, sig))
    for mode, a in runs:
        G, r0 = make(mode, a)
        B, b, d, c, R = contact(G)
        print("%-5s a=%.6f r0=%.6f z(0.7523)=%.4e  theta=%.13f  theta-theta_inf=%+.3e  b=%.6f delta=%.6f cross=%.6f" % (
            mode, a, r0, float(G(0.7523)), B, B - TH, b, d, c)); sys.stdout.flush()

def make_gamma(gam):
    # hypothetical twisted-family estimate saving Z^(gamma tau) on the Cauchy count; best rho' in {0, 1}: T = max(T_0, gamma/(2k)); start pinned at 7/9
    def thr(r, z): return (1.25 * r - 5 / 6) / (r - 0.5 - z)
    T = lambda r, z: max((1 / thr(r, z) - 2) / 6, gam / (2 * thr(r, z)))
    f = lambda r, z: (z - T(r, z)) / r
    r0 = 7 / 9; n = 120000; r1 = 0.725; h = (r1 - r0) / n; rs = [r0]; zs = [0.0]; r, z = r0, 0.0
    for _ in range(n):
        k1 = f(r, z); k2 = f(r + h / 2, z + h * k1 / 2); k3 = f(r + h / 2, z + h * k2 / 2); k4 = f(r + h, z + h * k3)
        z += h * (k1 + 2 * k2 + 2 * k3 + k4) / 6; r += h; rs.append(r); zs.append(z)
    rs = np.array(rs[::-1]); zs = np.array(zs[::-1])
    return (lambda w: np.interp(w, rs, zs)), r0
if __name__ == '__main__' and len(sys.argv) > 1 and sys.argv[1] == 'gamma':
    TH = 0.8733458561067185768
    for gam in (0.03, 0.05, 0.1, 0.2):
        G, r0 = make_gamma(gam); B, b, d, c, R = contact(G)
        print("gamma=%.3f z(0.7523)=%.4e theta=%.13f theta-theta_inf=%+.3e cross=%.6f" % (gam, float(G(0.7523)), B, B - TH, c)); sys.stdout.flush()
