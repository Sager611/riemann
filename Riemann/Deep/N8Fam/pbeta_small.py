# bandwidth lam <= 1 (single prime modulus family of Hua-Yang at T = q^theta: lam = 1/(1+theta); fixed chi: lam = 1)
# plain optimum p(lam; F=alpha) vs dilated paper window with scaled certificate (x, lam*b)
import numpy as np, sys
sys.path.insert(0, '.')
from pbeta import pbeta
from cfam import C_model, cos_window, A2, X_CERT, B_CERT
f2 = cos_window(A2)
print('theta  lam      p(lam) plain-opt   paper-window plain   dilated cert   gain')
for th in [0, 0.01, 0.05, 0.1, 0.25, 0.5]:
    lam = 1 / (1 + th)
    p = pbeta(lam)[0]
    plain = 2 - C_model(f2, lam, 'gue')
    cert = (2 - C_model(f2, lam, 'gue') - lam * B_CERT) / (1 - X_CERT)
    print('%-5g %.5f  %.6f          %.6f            %.6f      %+.6f' % (th, lam, p, plain, cert, cert - p))
