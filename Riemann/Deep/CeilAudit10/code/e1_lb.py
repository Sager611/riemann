# e1_lb.py — ceilaudit10: (1) check nu_DI(delta) = (4/7)[(3t-1)/2 - (1-t) delta] at theta = 0.75, delta = 0.2 with e1rate's exponents.py;
# (2) LB = max{0,c,(c+2k)/3,h} at nu = E1 vs nu = (2/7)(3t-1) (full phase keeping), kappa from kappa_opt.opt (deg 6).
import sys
sys.dont_write_bytecode = True
sys.path.insert(0, 'Riemann/Deep/E1Rate')
import exponents as X, kappa_opt as K
from scipy.optimize import brentq
th, d = 0.75, 0.2
def eII_d(th_, pv, a, b, n):
    etr, esep, x0 = X.bounds(th_, pv, a+b, n); return etr + d*max(0.0, x0) + X.gDI(a, b, n, pv)
X.eII = eII_d
print('check theta=0.75 delta=0.2: nu_K=%.5f predicted %.5f (E1 %.5f)' % (X.nuK(th), 4/7*((3*th-1)/2-(1-th)*d), (3*th-1)/4), flush=True)
kap = {}
def k(nu):
    if nu not in kap: kap[nu] = K.opt(nu, deg=6)[0]
    return kap[nu]
for t in (0.51, 0.55, 0.58, 0.6, 0.61, 0.62, 0.63):
    e1 = (3*t-1)/4; na = 2*(3*t-1)/7
    print('theta=%.3f c=%.5f | E1 nu=%.4f kappa=%.6f LB=%.6f | kept nu=%.4f kappa=%.6f LB=%.6f' % (t, K.cth(t), e1, k(e1), K.lb(t, k(e1)), na, k(na), K.lb(t, k(na))), flush=True)
f = lambda t: k(round(2*(3*t-1)/7, 6)) - K.cth(t)
print('kappa((2/7)(3t-1)) = c(t) at theta = %.4f' % brentq(f, 0.6, 0.64, xtol=1e-4), flush=True)
