import sys, numpy as np
sys.path.insert(0, 'Riemann/Deep/E1Rate')
import kappa_opt as K
for xi, cap in ((25/256, 0.07005633), (1/8, 0.08971026)):
    for deg in (6, 10):
        k, p = K.opt(xi, deg=deg, R0=None if deg == 6 else 9.87 if xi < 0.1 else None)
        print('xi=%.8f deg=%d kappa=%.10f R=%.5f (certified cap %.8f; below: %s)' % (xi, deg, k, p[0], cap, k <= cap), flush=True)
