# extra.py (explore23, T2): (1) direct glitch peaks far above the 'Nyquist' prime sqrt(T/2pi); (2) local second-difference statistics
import json, numpy as np
from math import isqrt
PI = np.pi
g = np.load('zeta_zeros.npy')
def win(gs, T1, T2):
    x = (gs - T1)/(T2 - T1)
    return np.where((x > 0) & (x < 1), np.sin(PI*x)**2, 0.0)
out = {}
for (T1, T2) in [(2000.0, 20000.0), (2000.0, 4000.0)]:
    w = win(g, T1, T2); Nw = w.sum(); Lb = (w*np.log(g/(2*PI))).sum()/Nw
    rows = []
    for p in [17, 19, 101, 499, 997, 4999, 10007, 49999, 100003]:
        s = (w*np.exp(1j*np.log(p)*g)).sum()/Nw
        rows.append([p, round(float(np.sqrt(T1/(2*PI))), 1), round(float(s.real), 6), round(float(s.imag), 6), round(float(-np.log(p)/np.sqrt(p)/Lb), 6)])
    out['beyond_nyquist_%d_%d' % (T1, T2)] = rows
# local statistics: normalized gaps and second differences; moved pair = double ordinate (gap 0) regardless of delta
dgap = np.diff(g)
loc = 2*PI/np.log(0.5*(g[1:] + g[:-1])/(2*PI))
ng = dgap/loc
order = np.argsort(ng)[:6]
out['smallest_normalized_gaps_[n, g_n, g_n+1, gap/mean]'] = [[int(i + 1), round(float(g[i]), 6), round(float(g[i + 1]), 6), round(float(ng[i]), 5)] for i in order]
d2 = g[2:] - 2*g[1:-1] + g[:-2]
nd2 = d2/loc[1:]
out['second_diff_normalized_[min,max,rms]'] = [round(float(nd2.min()), 4), round(float(nd2.max()), 4), round(float(np.sqrt(np.mean(nd2**2))), 4)]
json.dump(out, open('extra.json', 'w'), indent=1)
print(json.dumps(out))
