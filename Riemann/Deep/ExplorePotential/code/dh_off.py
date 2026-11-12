# dh_off.py: add the known off-line DH zeros (ExploreHelio dh_zeros.json, T <= 420) to the on-line sums:
# x^{-1/2} sum_rho x^rho = S_on + sum_{off right} x^{i gamma} 2 cosh(delta log x); compare with Landau -(T/2pi) c(x)/sqrt(x).
import json, math, sys, numpy as np
sys.path.insert(0, '.')
from dh_ft import coeffs
d = json.load(open(sys.argv[1])); off = json.load(open(sys.argv[2]))['offline']; T = d['TMAX']
on = np.array(d['online']); a, c = coeffs(40)
print('T', T, 'on-line', len(on), 'off-line right members', len(off), [(round(b, 3), round(g, 2)) for b, g in off])
print('| x | c(x) | Landau P | Re S_on - P | Re S_all - P | off share Re(S_all - S_on)/P | 2 E_off (Re) |')
res_on = []; res_all = []
for x in range(2, 41):
    L = math.log(x)
    S_on = np.exp(1j*on*L).sum()
    S_off = sum(np.exp(1j*g*L)*2*math.cosh((b - 0.5)*L) for b, g in off)
    E2 = sum(np.exp(1j*g*L)*2*(math.cosh((b - 0.5)*L) - 1) for b, g in off)
    P = -(T/(2*math.pi))*c[x]/math.sqrt(x)
    S_all = S_on + S_off
    res_on.append(abs(S_on.real - P)); res_all.append(abs(S_all.real - P))
    sh = (S_off.real/P) if abs(P) > 1e-9 else float('nan')
    print('| %d | %+.4f | %+.1f | %+.1f | %+.1f | %s | %+.2f |' % (x, c[x], P, S_on.real - P, S_all.real - P, ('%+.2f' % sh) if abs(P) > 1e-9 else '-', E2.real))
print('mean |residual| on-line only %.2f, with off-line zeros %.2f' % (np.mean(res_on), np.mean(res_all)))
