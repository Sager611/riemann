"""Float values of sup_v lam(kappa)/kappa^2 (r = 1, L2 closure) from the Arb-exact matrices of gcert.Ops at large kappa
(the FFT pipeline of gcore loses accuracy there). Not a certificate. usage: evid_arb.py kappa [N]"""
import sys, time, math, resource
import numpy as np
from scipy.linalg import eigh
from scipy.optimize import minimize_scalar
from flint import arb
import gcert as GC
kap = GC.dec(sys.argv[1]); N = int(sys.argv[2]) if len(sys.argv) > 2 else 16
t0 = time.time()
gam = int(float(2 * GC.A_(kap) * arb.pi().mid() + 3)) + 1
D = max(50, int(2.5 * gam / 4 + 34))
GC.set_prec(128 + int(D * math.log2(max(gam / 4, 2))) + 32)
ops = GC.Ops(kap, N, D, gam, Nw=N)
F = {hp: np.array([[float(x.mid()) for x in row[:N]] for row in ops.A(*hp).tolist()[:N]]) for hp in [('a', 0), ('a', 1), ('a', 2), ('b', 0), ('b', 1), ('c', 0)]}
M0 = F[('a', 0)]; MS = F[('a', 1)] + F[('b', 0)]; MSS = F[('a', 2)] + 2 * F[('b', 1)] + F[('c', 0)]
top = lambda v: eigh(M0, v * v * M0 - 2 * v * MS + MSS, eigvals_only=True)[-1]
rr = minimize_scalar(lambda v: -top(v), bounds=(0.6, 2.0), method='bounded', options={'xatol': 1e-9})
k2 = float(kap) ** 2
print('r=1 kappa=%s N=%d D=%d: sup_v lam = %.6f at v=%.5f ; lam/k^2 = %.7f [%.0fs, peak %.0f MB]' % (sys.argv[1], N, D, -rr.fun, rr.x, -rr.fun / k2, time.time() - t0, resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2 ** 20))
