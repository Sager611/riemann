"""Arb check that the Jacobi basis of gcert_r.OpsR is orthonormal in L2(X^b dX) on [0,1/4]: QT W_b QT^T contains the identity."""
import sys
from flint import arb
import gcert as GC
import gcert_r as GR
GC.set_prec(289)
kap = GC.dec('3.31'); gam = int(float(2 * GC.A_(kap) * arb.pi().mid() + 3)) + 1
ops = GR.OpsR(kap, 72, 50, gam, 2)
n = ops.Nw
E = (ops.QT * GC.sub(ops.Wq[0], n, n) * ops.QT.transpose()).tolist()
ok = all(E[i][j].contains(1 if i == j else 0) for i in range(n) for j in range(n))
mr = max(float(E[i][j].rad()) for i in range(n) for j in range(n))
print('Nw=%d: QT W_b QT^T contains I: %s (max radius %.2e)' % (n, ok, mr))
