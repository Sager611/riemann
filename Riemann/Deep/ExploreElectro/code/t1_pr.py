# t1_pr.py (explore22): positive-real (passivity) test of the 'xi network' Z(w) = F'/F(1/2 + w), Re w = eps > 0,
# i.e. Foster/Lagarias: Z is PR on Re w > 0 iff F has no zeros right of the line (RH for xi).
# Scans G(eps, t) = Re F'/F(1/2 + eps + i t) for t in [t0, t1]; reports min, argmin and the negative intervals.
# Usage: t1_pr.py <xi|mv|dh> t0 t1 dt eps1,eps2,...
import sys, json, time
sys.path.insert(0, '.')
import mpmath as mp
import ecommon as E
mp.mp.dps = 18
kind = sys.argv[1]; t0 = float(sys.argv[2]); t1 = float(sys.argv[3]); dt = float(sys.argv[4])
epss = [float(e) for e in sys.argv[5].split(',')]
L = E.LOGD[kind]
for eps in epss:
    ts0 = time.time()
    n = int(round((t1 - t0) / dt))
    vals = []
    for k in range(n + 1):
        t = t0 + k * dt
        vals.append((t, float(mp.re(L(mp.mpc(0.5 + eps, t))))))
    gmin = min(vals, key=lambda x: x[1])
    neg = []; cur = None
    for t, g in vals:
        if g < 0:
            if cur is None: cur = [t, t, g]
            else: cur[1] = t; cur[2] = min(cur[2], g)
        elif cur is not None:
            neg.append([round(cur[0], 3), round(cur[1], 3), round(cur[2], 4)]); cur = None
    if cur is not None: neg.append([round(cur[0], 3), round(cur[1], 3), round(cur[2], 4)])
    print(json.dumps({'kind': kind, 'eps': eps, 't0': t0, 't1': t1, 'dt': dt, 'min_G': round(gmin[1], 6), 'argmin_t': round(gmin[0], 3),
                      'negative_intervals[t_a,t_b,minG]': neg, 'secs': round(time.time() - ts0, 1)}), flush=True)
