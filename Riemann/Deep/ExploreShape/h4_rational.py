"""explore53: h4 of the Arb-validated rational J = 4 form (logs/cert4_k3.34.json) on the degree-64 data at kappa = 3.34
(float mean values; phi terms telescope and do not enter).  usage: h4_rational.py"""
import sys, os, json
from fractions import Fraction as Fr
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import joint as Jt
from scipy.optimize import minimize_scalar
C = json.load(open(os.path.join(HERE, 'logs', 'cert4_k3.34.json')))
a = [float(Fr(x)) for x in C['a']]
kap = 3.34; v0 = 1.3617
C0 = Jt.central_forms(2, 0.2499, 64, kap, 8, v0)
def h(dv):
    Cs = Jt.shift(C0, dv, kap)
    R = sum(a[j] * Cs[2 * j] for j in range(1, 5))
    return -a[0] * Jt.lmax(Cs[0], R)[0]
rr = minimize_scalar(lambda d: -h(d), bounds=(-0.005, 0.005), method='bounded', options={'xatol': 1e-10})
fw = lambda dv: Jt.lmax(*[Jt.shift(C0[:3], dv, kap)[i] for i in (0, 2)])[0]
rw = minimize_scalar(lambda d: -fw(d), bounds=(-0.01, 0.01), method='bounded', options={'xatol': 1e-10})
print('kappa = 3.34, degree 64, r = 2, theta = 0.2499: Wirtinger h1 = %.8f ; rational Arb-valid J=4 form: h4 = %.8f at v = %.6f' % (-rw.fun, -rr.fun, v0 + rr.x))
