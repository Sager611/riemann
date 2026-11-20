"""Bui-Hall Table 2 check: max over linear P and V of h_k at lam = C(k), eta = 0.5 (printed) and eta = lam/4 (midpoint)."""
import sys, json, numpy as np, scipy.linalg as sl, zk
from scipy.optimize import minimize_scalar
C = {1: 1.9, 2: 1.606, 3: 1.451, 4: 1.365, 5: 1.306, 6: 1.265, 7: 1.231, 8: 1.205, 9: 1.184, 10: 1.167}
ks = [int(x) for x in sys.argv[1].split(',')]; d = int(sys.argv[2]) if len(sys.argv) > 2 else 1
def best(k, lam, eta, d, nodes):
    M, herm, imag = zk.forms(k, lam, eta, 1.0, d, 2, nodes=nodes)
    def h(V):
        t = V - 1.0; C2 = M[2] - 2 * t * M[1] + t * t * M[0]
        return -sl.eigh(M[0], lam**2 * C2, eigvals_only=True)[-1]
    Vs = np.linspace(0.0, 2.0, 41); hs = [-h(V) for V in Vs]; i = int(np.argmax(hs))
    r = minimize_scalar(h, bounds=(Vs[max(i-1,0)], Vs[min(i+1,40)]), method='bounded', options={'xatol': 1e-8})
    return -r.fun, r.x
for k in ks:
    nodes = (8, 8, 6, 8, 6, 12 + k, 12 + k)
    out = [k, C[k]]
    for eta in (0.5, C[k] / 4):
        hv, V = best(k, C[k], eta, d, nodes); out += [eta, hv, V - 1]
    print('k=%d C(k)=%.3f | eta=%.4f: max h=%.6f (v=%.4f) | eta=%.4f: max h=%.6f (v=%.4f)' % tuple(out), flush=True)
