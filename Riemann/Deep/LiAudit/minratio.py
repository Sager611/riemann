# min/max over all N of the Theorem 1 and Theorem 2 ratios (R from the FFT table).
import sys, numpy as np
from math import pi, sqrt
R = np.load(sys.argv[1]); c1 = sqrt(3/5)*(4/5)**2/(2*pi*sqrt(2))
cs = np.cumsum(R**2); N = np.arange(100, 100001)
q = cs[N]/(c1*N**1.5*np.sqrt(np.log(N)))
print('Thm1 ratio over 100<=N<=1e5: min %.4f at N=%d, max %.4f at N=%d' % (q.min(), N[q.argmin()], q.max(), N[q.argmax()]))
for a in (1000, 3000, 10000):
    s = q[N >= a]; print('  over %d<=N<=1e5: min %.4f at N=%d' % (a, s.min(), N[N >= a][s.argmin()]))
M = np.arange(100, 50001); w = (cs[2*M]-cs[M])/(M**1.5*np.sqrt(np.log(M)))
print('Thm2 window ratio /(N^1.5 sqrt log N) over 100<=N<=5e4: min %.4f (N=%d), max %.4f (N=%d)' % (w.min(), M[w.argmin()], w.max(), M[w.argmax()]))
