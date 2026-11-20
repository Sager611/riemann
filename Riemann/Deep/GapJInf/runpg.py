"""usage: runpg.py forms.npz amp.json Jmin Jmax [K] [s1 s2 ...]  (s = dilation of nu: kappa_eff = kappa*s)"""
import sys, os, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from ampmom import load
import pgcut as PG
import json

npz, amj, J1, J2 = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
KM = int(sys.argv[5]) if len(sys.argv) > 5 else 3000
ss = [float(x) for x in sys.argv[6:]] or [1.0]
sv, am, kap, v = load(npz, amj)
m0 = sv[0::2]
for s in ss:
    m = m0 * s ** (2 * np.arange(len(m0)))
    for J in range(J1, J2 + 1):
        t0 = time.time()
        P = PG.PerGap(m, J, KM=KM)
        val, xt, it, hist = P.run()
        a = xt[:J + 1]; ph = xt[J + 1:]
        print('kappa_eff=%.4f (s=%.4f) J=%d: h_J = %.7f  [%d rounds, %d cuts, last worst %.1e, %.1f s] a/a1 = %s phi = %s' % (
            kap * s, s, J, val, it, len(P.cuts), hist[-1][2], time.time() - t0, ' '.join('%.5g' % x for x in a / a[1]), ' '.join('%.4g' % x for x in ph / a[1])), flush=True)
        fn = os.path.join(HERE, 'logs', 'form_k%.4f_J%d.json' % (kap * s, J))
        json.dump({'kappa_eff': kap * s, 'J': J, 'h': val, 'a': list(a / a[1]), 'phi': {'%d,%d' % (2 * mm + 1, 2 * l): float(p / a[1]) for (mm, l), p in zip(P.tkeys, ph)},
                   'rounds': it, 'last_worst': hist[-1][2], 'moments': list(m[:J + 1])}, open(fn, 'w'), indent=1)
