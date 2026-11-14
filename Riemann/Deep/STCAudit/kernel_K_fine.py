# Fine scans of (K) in the active region: b in [0.45, 0.5] x t in [-3, 3] (501 x 60001) and the edges b = 1/2, b = 0
# for t in [-60, 60] (1.2e6 points), all six kernels; reuses Hfun of kernel_K_test.py.
import json, math
import numpy as np
src = open('Riemann/Deep/STCAudit/kernel_K_test.py').read()
exec(src.split('res = {}')[0])   # defines K, Hfun
for name, k in K.items():
    m1 = np.inf
    for b0 in np.linspace(0.45, 0.5, 501):
        t = np.linspace(-3, 3, 30001); t = t[t != 0]
        h = Hfun(k, np.full_like(t, b0), t, scaled=True); i = h.argmin()
        if h[i] < m1: m1, a1 = h[i], (b0, t[i])
    T = np.linspace(-60, 60, 1200001); T = T[T != 0]
    m5, t5, m0 = np.inf, 0, np.inf
    for t in np.array_split(T, 24):
        e5 = Hfun(k, np.full_like(t, 0.5), t, scaled=True); e0 = Hfun(k, np.zeros_like(t), t, scaled=True)
        if e5.min() < m5: m5, t5 = e5.min(), t[e5.argmin()]
        m0 = min(m0, e0.min())
    print('%-5s fine active-region min %.4e at b=%.4f t=%.5f | edge b=1/2 min %.4e at t=%.4f | edge b=0 min %.4e' % (name, m1, a1[0], a1[1], m5, t5, m0))
