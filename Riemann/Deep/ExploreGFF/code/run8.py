# run8.py (explore17): T1 margins for Davenport-Heilbronn (odd gamma factor, conductor 5, coefficients b_n of -f'/f; no pole)
import os, sys, math, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gw import *
kap, a, b = dh_coeffs(200)
print('kappa %.6f; b_2..b_9:' % kap, ' '.join('%.3f' % b[n] for n in range(2, 10)))
print('  c     M1_DH(LK+V, primes by norm)   M3-type (w_n half-overlap)   exact lam_min(Q_DH) (float64, N=40)')
for c in (1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 7.0, 10.0):
    out = []
    for sec in ('even', 'odd'):
        B = Basis(c, 40, sec); n = len(B.idx)
        LK = B.jump(kchi_odd, 600); V = Vdh(B.L)*np.eye(n)
        PR = B.primes({m: b[m]/math.sqrt(m) for m in range(2, 200)})
        pn = abs(np.linalg.eigvalsh(PR)).max() if np.any(PR) else 0.0
        a0 = np.linalg.eigvalsh(LK + V)[0]
        sw = sum(abs(b[m])/math.sqrt(m) for m in range(2, int(math.ceil(c))) if m < c)
        out.append((a0 - pn, a0 - sw, np.linalg.eigvalsh(LK + V + PR)[0]))
    m = [min(out[0][k], out[1][k]) for k in range(3)]
    print('  %-5g %+.4f                        %+.4f                       %+.3e' % (c, *m), flush=True)
