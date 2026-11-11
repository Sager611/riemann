"""Float search: max kappa with max_{P,v} h1 > 1 in Bui-Milinovich's argument. usage: bmopt.py r theta n D1,D2,..."""
import sys, time, resource
sys.path.insert(0, '.')
import numpy as np
import bmfloat as bf
r = int(sys.argv[1]); th = float(sys.argv[2]); n = int(sys.argv[3]); Ds = [int(s) for s in sys.argv[4].split(',')]
t0 = time.time()
Dmax = max(Ds)
cache = {}
def M_of(k):
    if k not in cache:
        cache[k] = bf.mats(th, r, k, Dmax, n)[0]
    return cache[k]
def best(k, D):
    M = [A[:D, :D] for A in M_of(k)]
    vs = np.linspace(0.6, 2.4, 91)
    hs = [bf.h1max(M, k, v)[0] for v in vs]
    i = int(np.argmax(hs)); v = vs[i]; h = hs[i]; dv = 0.02
    while dv > 1e-6:
        moved = False
        for vv in (v - dv, v + dv):
            hh = bf.h1max(M, k, vv)[0]
            if hh > h: h, v, moved = hh, vv, True
        if not moved: dv /= 2
    return h, v, bf.h1max(M, k, v)[1]
for D in Ds:
    lo, hi = 2.9, 3.8
    h_lo = best(lo, D)[0]
    assert h_lo > 1, h_lo
    for _ in range(22):
        mid = 0.5 * (lo + hi)
        if best(mid, D)[0] > 1: lo = mid
        else: hi = mid
    h, v, P = best(lo, D)
    print('r=%d theta=%.4f n=%d deg P=%d: kappa_max %.5f (v=%.4f, h1=%.6f) Pleg=%s (%.0fs, peak %.0f MB)' % (
        r, th, n, D - 1, lo, v, h, np.array2string(P, precision=4, max_line_width=250), time.time() - t0,
        resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2 ** 20), flush=True)
    cache.clear()
