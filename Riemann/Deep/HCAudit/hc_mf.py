# hcaudit: M_f = sup_u u^2 |f^(u)| for the cert768 f (EPS = 1/1000), the constant of FULLPROOF §2(2) and §6(a).
# Usage: python3 hc_mf.py cert768.json 16604 9.90343   (N0 and L from hc_zeros768.log)
import sys, json, numpy as np
from fractions import Fraction as Fr
d = json.load(open(sys.argv[1])); D = d['D']; s = float(Fr(d['s'])); N0 = int(sys.argv[2]); L = float(sys.argv[3])
B = [Fr(b) for b in d['B']]; B[0] += Fr(1, 1000); Bf = np.array([float(b) for b in B]); k = np.sqrt(2*np.pi)
def fh(u):
    v = k*s*np.asarray(u, float); hs = [np.ones_like(v), 2*v]
    for m in range(1, 2*D): hs.append(2*v*hs[-1] - 2*m*hs[-2])
    return s*sum(Bf[n]*(-1)**n*hs[2*n] for n in range(D+1))*np.exp(-v*v/2)
u = np.linspace(0, 12, 120001); y = u*u*np.abs(fh(u)); i = y.argmax()
print('M_f = sup u^2|f^(u)| ~ %.5f at u = %.4f; 4 pi^2 M_f N0/L^2 = %.1f' % (y[i], u[i], 4*np.pi**2*y[i]*N0/L**2))
