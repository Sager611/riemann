# dh_alllaw.py (explore35): all-ordinate law for DH: S_all = S_on + sum_off 2 x^{i gamma} vs its Landau value; residual = -2 E_off (Lemma 1).
import json, math, sys, numpy as np
sys.path.insert(0, '.')
from dh_ft import coeffs
on = np.array(json.load(open(sys.argv[1]))['online']); off = np.array(json.load(open(sys.argv[2]))['offline']); T = 2500.0
bet, gam = off[:, 0], off[:, 1]; a, c = coeffs(40)
print('| x | Landau P | S_all - P (all ordinates, unweighted) | -2 E_off (Lemma 1) | (S_all - P)/P |')
for x in (2, 3, 4, 6, 7, 9, 11, 13, 14, 19, 29, 31, 39):
    L = math.log(x); P = -(T/(2*math.pi))*c[x]/math.sqrt(x)
    S_all = np.exp(1j*on*L).sum().real + (2*np.exp(1j*gam*L)).sum().real
    E2 = -(2*np.exp(1j*gam*L)*(np.cosh((bet - 0.5)*L) - 1)).sum().real
    print('| %d | %+.1f | %+.1f | %+.1f | %+.3f |' % (x, P, S_all - P, E2, (S_all - P)/P))
