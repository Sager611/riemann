"""Converged-quadrature scan: for kappa grid, max over v and P (deg <= Dmax-1) of h1; usage bmscan.py theta n Dmax k1,k2,..."""
import sys, time, resource, json
sys.path.insert(0, '.')
import numpy as np
import bmfloat as bf
th = float(sys.argv[1]); n = int(sys.argv[2]); Dmax = int(sys.argv[3]); ks = [float(s) for s in sys.argv[4].split(',')]
r = int(sys.argv[5]) if len(sys.argv) > 5 else 1
t0 = time.time(); out = {}
def best_v(M, k):
    vs = np.linspace(1.0, 1.6, 61)
    hs = [bf.h1max(M, k, v)[0] for v in vs]; i = int(np.argmax(hs)); v, h, dv = vs[i], hs[i], 0.005
    while dv > 1e-7:
        moved = False
        for vv in (v - dv, v + dv):
            hh = bf.h1max(M, k, vv)[0]
            if hh > h: h, v, moved = hh, vv, True
        if not moved: dv /= 2
    return h, v, bf.h1max(M, k, v)[1]
for k in ks:
    M, im = bf.mats(th, r, k, Dmax, n)
    row = {}
    for D in range(1, Dmax + 1):
        h, v, P = best_v([A[:D, :D] for A in M], k)
        row[D - 1] = (h, v, list(P))
    out[k] = row
    print('kappa=%.4f: ' % k + ' '.join('d%d:%.6f' % (d, row[d][0]) for d in row) + '  (%.0fs)' % (time.time() - t0), flush=True)
json.dump({str(k): v for k, v in out.items()}, open('bmscan_th%s_n%d_D%d_r%d.json' % (sys.argv[1], n, Dmax, r), 'w'))
print('peak %.0f MB' % (resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2 ** 20))
