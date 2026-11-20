"""Wirtinger (J = 1) optimum of Bui-Hall's criterion over V = v+1 and the amplifier P (degree d), fixed eta (or midpoint)."""
import sys, numpy as np, scipy.linalg as sl, zk
from scipy.optimize import minimize_scalar
k = int(sys.argv[1]); lam = float(sys.argv[2]); eta = sys.argv[3]; d = int(sys.argv[4])
nodes = tuple(int(x) for x in sys.argv[5].split(',')) if len(sys.argv) > 5 else (8, 8, 6, 8, 6, 12, 12)
eta = lam / 4 if eta == 'mid' else float(eta)
v0 = 1.0
M, herm, imag = zk.forms(k, lam, eta, v0, d, 2, nodes=nodes)
def h(V, ret=False):
    t = V - v0
    C2 = M[2] - 2 * t * M[1] + t * t * M[0]
    ev, U = sl.eigh(M[0], lam**2 * C2)
    return (ev[-1], U[:, -1]) if ret else -ev[-1]
Vs = np.linspace(-0.5, 2.5, 61); hs = [-h(V) for V in Vs]
i = int(np.argmax(hs)); r = minimize_scalar(h, bounds=(Vs[max(i-1,0)], Vs[min(i+1,60)]), method='bounded', options={'xatol': 1e-8})
hv, b = h(r.x, True); b = b / b[0]
print('k=%d lam=%.4f eta=%.4f d=%d: max h = %.6f at V=%.5f (v=%.5f), P coeffs (w^i) %s; at BH v=0.2: h=%.6f' % (k, lam, eta, d, hv, r.x, r.x - 1, np.round(b, 4).tolist(), -h(1.2)))
