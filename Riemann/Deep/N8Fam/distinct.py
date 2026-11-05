# Distinct zeros: 2610.08965 Thm main (first window, c = 67/20, (x1,x2,b) = (0.054350, 0.069433, 6*0.006219)),
# dilated to bandwidth lam:  N_d/N >= (3 - C(f_lam) + x2 - x1 - lam*b)/(2 - 2 x1 + x2)   [Thm general; mu,nu,r_h,r_k
# do not involve b or the scale, so they stay valid]. Compare with Johns Thm 1.4(a): (1 + p(beta;F1))/2.
import numpy as np, sys
sys.path.insert(0, '.')
from cfam import C_model, cos_window
from pbeta import pbeta
A1 = np.array([193305, -23105, 29166, 13952, -21009, 41052]) * 1e-6
x1, x2, al = 0.054350, 0.069433, 0.006219
b = 6 * al
f1 = cos_window(A1)
C1 = C_model(f1, 1.0, 'zeta')
print('first window C(f) = %.11f (paper 1.33270495828); zeta value %.7f (paper 1645064/1960733 = %.7f)'
      % (C1, (3 - 266541 / 200000 + x2 - x1 - b) / (2 - 2 * x1 + x2), 1645064 / 1960733))
johns = {0: 0.966141, 0.25: 0.956287, 0.5: 0.947492, 1: 0.932836, 2: 0.912177, 3: 0.898606, 5: 0.882074, 10: 0.863742}
print('kappa beta    Johns distinct (1+p)/2   paper first window + cert (best lam<=beta)   gain')
for k in [0, 0.5, 1, 2, 3, 5, 10, 20, 100, 1000]:
    beta = (2 + k) / (1 + k)
    p = pbeta(beta)[0]
    lams = np.linspace(max(0.9, beta - 0.4), beta, 41)
    vals = [(3 - C_model(f1, l, 'gue') + x2 - x1 - l * b) / (2 - 2 * x1 + x2) for l in lams]
    i = int(np.argmax(vals))
    print('%-5g %.4f  %.6f %s           %.6f (lam %.4f)                   %+.6f'
          % (k, beta, (1 + p) / 2, '(J %.6f)' % johns[k] if k in johns else '            ', vals[i], lams[i], vals[i] - (1 + p) / 2), flush=True)
