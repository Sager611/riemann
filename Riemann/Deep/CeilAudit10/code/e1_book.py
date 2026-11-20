# e1_book.py — ceilaudit10: exponent bookkeeping for e1rate's barrier. LB = max{0, c, (c+2k)/3, h(theta;k)} exceeds Wang's c(theta) iff k > c(theta)
# (identity h(theta; c) = c, h increasing in k). Tabulate e1rate's one-piece kappa(nu) (kappa_opt.opt, deg 4) and compare nu*(theta) = kappa^{-1}(c(theta))
# with the class limits: E1 (3t-1)/4; separated DI (4/7)(2t-1); separated nu1=1: 2t-1; phase kept (Kuznetsov with Bessel weight + Lemma I type I) min(4t/7, (2/3)(2t-1)); SL cap 1/2.
import sys, numpy as np
sys.path.insert(0, 'Riemann/Deep/E1Rate')
import kappa_opt as K
nus = [0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.5714285714, 0.65, 0.75, 0.875, 1.0]
ks = [0.052919, 0.101893, 0.142568, 0.180739, 0.218596, 0.256394, 0.293593, 0.329492, 0.363501, 0.408059, 0.451090, 0.496839, 0.541167, 0.573609]  # from e1_book.log run 1 (K.opt deg 4)
if len(sys.argv) > 1 and sys.argv[1] == 'recompute':
    ks = []
    for nu in nus:
        k, p = K.opt(nu, deg=4); ks.append(k); print('kappa(%.4f) = %.6f  R = %.4f' % (nu, k, p[0]), flush=True)
ks = np.array(ks); nus = np.array(nus)
assert np.all(np.diff(ks) > 0)
kinv = lambda k: np.interp(k, ks, nus, left=np.nan, right=np.inf)
c = K.cth
for th in (5/7, 9/11):
    print('c(%.6f) = %.6f' % (th, c(th)))
# check identity h(th, c) = c
print('identity check max|h(th,c(th)) - c(th)| over th in [0.56,0.99]:', max(abs(K.h(t, c(t)) - c(t)) for t in np.linspace(0.56, 0.99, 200)))
# theta where c = kappa(1/2), kappa(4/7), kappa(1)
from scipy.optimize import brentq
for nm, kv in (('kappa(1/2)', ks[8]), ('kappa(4/7)', ks[9]), ('kappa(1)', ks[-1])):
    print('c(theta) = %s = %.6f at theta = %.6f' % (nm, kv, brentq(lambda t: c(t) - kv, 0.52, 0.999)))
print('theta   c(theta)  nu*=kinv(c)  E1      sepDI   sep(nu1=1)  phase-kept  cap  -> any class with limit > max(nu*, E1)?')
for th in np.arange(0.55, 1.0001, 0.025):
    cc = c(th); ns = kinv(cc) if cc > 0 else 0.0
    E1 = (3*th-1)/4; sepDI = 4/7*(2*th-1); sep1 = 2*th-1; pk = min(4*th/7, 2/3*(2*th-1)); cap = 0.5
    best = max(min(sepDI, cap), min(pk, cap))
    print('%.3f  %8.4f  %9.4f  %7.4f %7.4f %8.4f %9.4f %6.2f   %s' % (th, cc, ns, E1, sepDI, sep1, pk, cap,
          'LB gain' if (best > max(ns, E1)) else ('nu gain only' if best > E1 else 'none')), flush=True)

# part 2: what the two phase-kept escape routes of e1rate (DERIVATION §5 (a), (b)) would give for LB, bookkeeping only (inputs not available)
print('--- LB = max{0,c,(c+2k)/3,h} with nu = E1, (a) max(E1, min((2/3)(2t-1), 4t/7, 1/2)), (b) type II only, phase kept: max(E1, min(4t/7, 1/2))')
kof = lambda nu: np.interp(nu, nus, ks)
for th in (0.51, 0.55, 0.575, 0.6, 0.625, 0.65, 0.66, 0.675, 0.7, 0.714285714, 0.75):
    E1 = (3*th-1)/4; na = max(E1, min(2/3*(2*th-1), 4*th/7, 0.5)); nb = max(E1, min(4*th/7, 0.5))
    print('theta=%.4f  c=%.4f | E1 nu=%.4f LB=%.6f | (a) nu=%.4f LB=%.6f | (b) nu=%.4f LB=%.6f' % (th, c(th), E1, K.lb(th, kof(E1)), na, K.lb(th, kof(na)), nb, K.lb(th, kof(nb))), flush=True)
th3 = brentq(lambda t: 4*t/7 - kinv(c(t)), 0.61, 0.72); print('route (b) helps LB exactly for theta < %.5f (4t/7 = nu*(t))' % th3)
