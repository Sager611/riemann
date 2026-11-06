import os, sys, time, math, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gw import *
t0 = time.time()
for c in (2.0, 3.0):
    for sec in ('even', 'odd'):
        for N, nq in ((30, 400), (30, 800), (50, 800)):
            P = weil_parts(c, N, sec, nq)
            Q = P['LK'] + P['V'] + P['POLE'] + P['PR']
            ev = np.linalg.eigvalsh(Q)
            print(c, sec, N, nq, 'lam0..2', ['%.6e' % e for e in ev[:3]], flush=True)
print('time', round(time.time() - t0, 1))
