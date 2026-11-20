import sys, numpy as np
sys.path.insert(0, 'Riemann/Deep/E1Rate')
import kappa_opt as K
for th in (0.51, 0.55, 0.6, 0.65, 0.66):
    E1 = (3*th-1)/4; nb = min(4*th/7, 0.5)
    k0, _ = K.opt(E1, deg=4); kb, _ = K.opt(nb, deg=4)
    print('theta=%.3f c=%.6f | E1 nu=%.4f kappa=%.6f LB=%.6f | (b) nu=%.4f kappa=%.6f LB=%.6f' % (th, K.cth(th), E1, k0, K.lb(th, k0), nb, kb, K.lb(th, kb)), flush=True)
