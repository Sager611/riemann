# e1phase lb_dls.py — effect on the zero-proportion lower bound LB = max{0,c,(c+2k)/3,h} (SL/E1 detector, kappa_opt.opt deg 6, read-only import)
# of nu_DLS = max{E1, min(2t-1, t/2)} (double large sieve, conditional on the spacing lemma), vs E1.
import sys
sys.dont_write_bytecode = True
sys.path.insert(0, 'Riemann/Deep/E1Rate')
import kappa_opt as K
from scipy.optimize import brentq
kap = {}
def k(nu):
    nu = round(nu, 6)
    if nu not in kap: kap[nu] = K.opt(nu, deg=6)[0]
    return kap[nu]
for t in (0.6, 0.602, 0.604, 0.606, 0.608, 0.61, 0.615):
    e1 = (3*t-1)/4; nd = max(e1, min(2*t-1, t/2))
    print('theta=%.3f c=%.5f | E1 nu=%.4f kappa=%.6f LB=%.6f | DLS nu=%.4f kappa=%.6f LB=%.6f' % (t, K.cth(t), e1, k(e1), K.lb(t, k(e1)), nd, k(nd), K.lb(t, k(nd))), flush=True)
f = lambda t: k(2*t-1) - K.cth(t)
print('kappa(2t-1) = c(t) at theta = %.4f' % brentq(f, 0.6005, 0.62, xtol=1e-4), flush=True)
