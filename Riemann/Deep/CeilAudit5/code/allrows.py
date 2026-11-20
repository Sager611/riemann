#!/usr/bin/env python3
# ceilaudit5 allrows.py: per row, LDP-predicted near-band credit of the finite head vs Raymond's credits; then the float chain endpoint.
# credit_i(logN) = torus max of the prime-power head polynomial - Cramer level at logN (near band [X, 1.1X]: logN ~ 27.4-30);
# composites kept at triangle; Raymond's credit = model triangle(t0) - first_hi_i.  Far band x > 1.1X: main slope A log(1.1)/4 must cover it.
import sys, math, json, os
import numpy as np
sys.path.insert(0, 'Riemann/Deep/DBNHigh'); sys.path.insert(0, '.scratch/ceilaudit5')
import model as M, torus as T
arc = sys.argv[1]; rows = M.load_rows(arc)
G = 2048; th = np.arange(G)*2*math.pi/G
def ldp_gap(c, n, logNs):
    pp = {}; comp = 0.0
    for j, m in enumerate(n):
        f = T.FAC[m]
        if len(f) == 1:
            (p, k), = f.items(); pp.setdefault(p, {})[k] = c[j]
        else: comp += abs(c[j])
    tabs = [sum(a*np.cos(k*th) for k, a in ks.items()) for ks in pp.values()]
    m0s = np.array([float(t.max()) for t in tabs]); F = np.array(tabs) - m0s[:, None]
    mx = float(m0s.sum())
    def Lam(lam):
        e = np.exp(lam*F); E = e.mean(axis=1); return float(lam*mx + np.log(E).sum()), float(mx + ((F*e).mean(axis=1)/E).sum())
    out = {}
    for logN in logNs:
        lo, hi = 0.0, 5000.0
        for _ in range(60):
            lam = 0.5*(lo + hi); s, ds = Lam(lam); I = lam*ds - s
            if I < logN: lo = lam
            else: hi = lam
        out[logN] = mx - Lam(lam)[1]
    return mx, comp, out
res = []
for i, row in enumerate(rows):
    n, c = T.coeffs(row, row['t0'], row['t1'])
    tri = float(np.abs(c).sum())
    mx, comp, gaps = ldp_gap(c, n, (27.4, 30.0))
    ray = tri - row['A']
    slope_far = row['main']*math.log(1.1)/4
    res.append(dict(row=i, t0=row['t0'], tri=tri, A=row['A'], ray_credit=ray, torus_max_pp=mx, comp_tri=comp, gap274=gaps[27.4], gap30=gaps[30.0], far_slope=slope_far))
    print('%2d t0 %.5f tri %.5f A %.5f ray %.5f | pp torus max %.5f comp %.5f | LDP gap 27.4 %.5f 30 %.5f | far slope %.3f' % (i, row['t0'], tri, row['A'], ray, mx, comp, gaps[27.4], gaps[30.0], slope_far), flush=True)
json.dump(res, open('.scratch/ceilaudit5/allrows.json', 'w'), indent=0)
# chain endpoint
sys.argv = [sys.argv[0], arc, '/dev/null']
src = open('Riemann/Deep/DBNHigh/chain.py').read(); src = src[:src.index('t0 = time.time()')]
exec(src)
E0 = run()[0]; print('base endpoint', repr(E0), flush=True)
for key in ('gap274', 'gap30'):
    dS = [max(0.0, r[key] - max(0.0, r['ray_credit'])) for r in res]
    E = run(dS=dS)[0]; print(key, 'new credit beyond Raymond: mean dS %.4f; endpoint %r gain %.3e' % (sum(dS)/67, E, E0 - E), flush=True)
    dS2 = [0.5*x for x in dS]; E2 = run(dS=dS2)[0]; print(key, 'half of it: endpoint %r gain %.3e' % (E2, E0 - E2), flush=True)
