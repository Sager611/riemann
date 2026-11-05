# test_identity.py: lattice identity at P = 2L, 3L, 5L vs Galerkin (galerkin_fast, independent code) for cosine-basis f.
import sys, math, numpy as np
sys.path.insert(0, '/Users/adriansager/playground/riemann/Riemann/Deep/ChekAudit/evidence')
from weil_lat import *
from flint import ctx
import galerkin_fast as gf
c = int(sys.argv[1]) if len(sys.argv) > 1 else 5
L = math.log(c); at = atoms(c)
ctx.prec = 200
N = 6
M = gf.build(c, N)
Mf = np.array([[float(M[i, j].mid()) for j in range(N + 1)] for i in range(N + 1)])
x = np.array([0.0, 1.0, 0.3, -0.2, 0.1, 0.05, 0.02])
Qg = x @ Mf @ x
om = [2 * math.pi * j / L for j in range(N + 1)]
Nj = [math.sqrt(1 / L)] + [math.sqrt(2 / L)] * N
def F(r):
    r = np.asarray(r, dtype=float); tot = np.zeros_like(r)
    for j in range(N + 1):
        a = r - om[j]; b = r + om[j]
        sa = np.where(np.abs(a) < 1e-12, L / 2, np.sin(a * L / 2) / np.where(np.abs(a) < 1e-12, 1, a))
        sb = np.where(np.abs(b) < 1e-12, L / 2, np.sin(b * L / 2) / np.where(np.abs(b) < 1e-12, 1, b))
        tot += x[j] * Nj[j] * (sa + sb)
    return tot
print('c', c, 'Galerkin Q(f) =', Qg, ' ||x||^2 =', x @ x)
for fac in [2.0, 3.0, 5.0]:
    P = fac * L
    for K in [20000, 200000, 2000000]:
        k = np.arange(-K, K + 1)
        r = 2 * math.pi * k / P
        val = np.sum(m_sym(r, L, at) * F(r) ** 2) / P
        nrm = np.sum(F(r) ** 2) / P
        print(f'  P={fac}L K={K}: lattice Q = {val:.10f}  Parseval ||f||^2 = {nrm:.8f}')
