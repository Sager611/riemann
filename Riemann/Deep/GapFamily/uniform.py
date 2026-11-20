"""gapfamily: order-J criterion for the uniform spectral measure (GL(2) mixed second moments, 1410.7765):
m_{2j} = Lam^{2j}/(2j+1).  Bisection on Lam of the SDP-with-cuts value (< 0: criterion J proves Lambda_f > Lam, float)."""
import sys, json, os
import numpy as np
from sdpcut import sdp_cut
HERE = os.path.dirname(os.path.abspath(__file__))
Js = [int(x) for x in sys.argv[1].split(',')]; lo0 = float(sys.argv[2]); hi0 = float(sys.argv[3]); nb = int(sys.argv[4])
res = {}
for J in Js:
    lo, hi = lo0, hi0; best = None
    for it in range(nb):
        lam = (lo + hi) / 2
        m = [lam ** (2 * j) / (2 * j + 1) for j in range(J + 1)]
        val, a, ph, mnv, rds = sdp_cut(m, J, rounds=12)
        ok = val < -1e-7 and mnv > -1e-8
        print('J=%d Lam=%.6f SDP %.4g rounds %d minvalid %.2g -> %s' % (J, lam, val, rds, mnv, 'OK' if ok else 'no'), flush=True)
        if ok:
            lo = lam; best = {'lam': lam, 'val': val, 'a': [float(x) for x in a], 'phi': {str(k): v for k, v in ph.items()}, 'minvalid': mnv}
        else:
            hi = lam
    res[J] = {'lo': lo, 'hi': hi, 'best': best}
    print('J=%d: crossing in [%.6f, %.6f]' % (J, lo, hi), flush=True)
json.dump(res, open(os.path.join(HERE, 'logs', 'uniform_%s.json' % sys.argv[1].replace(',', '')), 'w'), indent=1)
