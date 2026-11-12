# Theorem 4: (smoothed mean of R - 7/4) * sum(phi)/log N for fixed bumps as N grows (should stay bounded).
import sys, numpy as np
from math import log
R = np.load(sys.argv[1])
def bump(x, a, b):
    y = np.zeros_like(x); s = (x > a) & (x < b); u = (2*x[s] - a - b)/(b - a); y[s] = np.exp(-1/(1 - u*u)); return y
for (a, b) in ((0.5, 1), (0.9, 1), (0.98, 1), (0.01, 0.02)):
    out = []
    for N in (2000, 5000, 10**4, 2*10**4, 5*10**4, 10**5):
        w = bump(np.arange(1, 100001)/N, a, b)
        if w.sum() < 1: continue
        e = (w*(R[1:] - 1.75)).sum()/log(N); out.append('%d:%+.2f' % (N, e))
    print('(%g,%g) [sum phi (R - 7/4)]/log N: ' % (a, b) + ' '.join(out))
