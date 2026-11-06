import sys, math, numpy as np
sys.path.insert(0, '/Users/adriansager/playground/riemann/Riemann/Deep/ExploreGFF/code')
from gw import *
kap, a, b = dh_coeffs(400)
W = {m: b[m]/math.sqrt(m) for m in range(2, 400)}
for c in [float(x) for x in sys.argv[1:]]:
    out = []
    for sec in ('even', 'odd'):
        for N in (40, 80):
            B = Basis(c, N, sec); n = len(B.idx)
            Q = B.jump(kchi_odd, 1200) + Vdh(B.L)*np.eye(n) + B.primes(W)
            out.append('%s N=%d %.3e' % (sec, N, np.linalg.eigvalsh(Q)[0]))
    print('c=%g  ' % c + '  '.join(out), flush=True)
