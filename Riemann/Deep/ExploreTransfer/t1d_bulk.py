# t1d_bulk.py -- explore3: eigenvalue cloud of L_b at sigma = 1/4, 1/2, 3/4 for t = 10, 20, 25:
# how many eigenvalues have |lambda| > 0.5, > 1, and the median modulus of the top-n cloud (semiclassical bulk).
import sys, time
import numpy as np
sys.path.insert(0, '.')
from t1_flint import build_M
LAM = np.pi**2 / (6*np.log(2))
if __name__ == '__main__':
    t0 = time.time()
    for t in (5.0, 10.0, 15.0, 20.0, 25.0):
        K = 100 if t <= 20 else 150
        line = []
        for s in (0.25, 0.5, 0.75, 1.0):
            ev = np.linalg.eigvals(build_M(complex(s, t), K, prec=64 + 3*K))
            m = np.sort(abs(ev))[::-1]
            n05 = int((m > 0.5).sum()); n1 = int((m > 1).sum()); n02 = int((m > 0.2).sum())
            line.append('s=%.2f: #>0.2=%d #>0.5=%d #>1=%d top5=%s' % (s, n02, n05, n1, ','.join('%.2f' % x for x in m[:5])))
        print('t=%.0f | ' % t + ' | '.join(line), flush=True)
    print('exp(Lambda*(1/2-s)) for s=1/4,1/2,3/4,1: ' + ', '.join('%.3f' % np.exp(LAM*(0.5 - s)) for s in (0.25, 0.5, 0.75, 1.0)))
    print('time %.1fs' % (time.time() - t0))
