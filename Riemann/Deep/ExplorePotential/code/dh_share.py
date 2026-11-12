# dh_share.py (explore35): on-line share of DH's Landau sum rule vs T (cumulative), with the located off-line zeros.
import json, math, sys, numpy as np
sys.path.insert(0, '.')
from dh_ft import coeffs
on = np.array(json.load(open(sys.argv[1]))['online']); off = np.array(json.load(open(sys.argv[2]))['offline'])
bet, gam = off[:, 0], off[:, 1]
a, c = coeffs(40)
xs = [2, 3, 4, 6, 7, 9, 11, 13, 14]
print('| T | on-line zeros | off-line zeros (2 x right members) | ' + ' | '.join('share x=%d' % x for x in xs) + ' | max sum-rule residual |')
for T in (420, 1000, 1500, 2000, 2500):
    o = on[on <= T]; m = gam <= T; shares = []; res = []
    for x in xs:
        L = math.log(x); P = -(T/(2*math.pi))*c[x]/math.sqrt(x)
        S_on = np.exp(1j*o*L).sum().real; S_off = (np.exp(1j*gam[m]*L)*2*np.cosh((bet[m] - 0.5)*L)).sum().real
        shares.append(S_on/P); res.append(abs(S_on + S_off - P))
    print('| %d | %d | %d | ' % (T, len(o), 2*m.sum()) + ' | '.join('%.2f' % s for s in shares) + ' | %.1f |' % max(res))
# zeta comparison line: share = 1 by RH verification; residual from zeta_ft logs
