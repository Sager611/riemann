#!/usr/bin/env python3
# ceilaudit5 selfcons.py: self-consistent near-band credit per row with Raymond's actual native share beta (ALL_67_ORIGINAL_AFFINE_NATIVE.json):
# requirement for all x >= X: Q(theta(x)) <= Qmax - cr + s*log(x/X), s = (1-beta)*M/4 (M = node mass v(0)); with sup_{[X,x]} Q ~ Qmax - gap(logN(x)),
# logN(x) = log(kappa*(x - X)) (Cramer level of the prime-power polynomial), the best credit is cr* = min_u [gap(log(kappa X (e^u - 1))) + s u].
import sys, math, json
import numpy as np
sys.path.insert(0, 'Riemann/Deep/DBNHigh'); sys.path.insert(0, '.scratch/ceilaudit5')
import model as M, torus as T
arc = sys.argv[1]; rows = M.load_rows(arc); B = json.load(open('.scratch/ceilaudit5/beta.json'))
X = 6000000185827.0
G = 2048; th = np.arange(G)*2*math.pi/G
lams = np.exp(np.linspace(math.log(0.5), math.log(8000.0), 260))
def rate_curve(c, n):
    pp = {}
    for j, m in enumerate(n):
        f = T.FAC[m]
        if len(f) == 1:
            (p, k), = f.items(); pp.setdefault(p, {})[k] = c[j]
    tabs = np.array([sum(a*np.cos(k*th) for k, a in ks.items()) for ks in pp.values()])
    m0 = tabs.max(axis=1); Fm = tabs - m0[:, None]; mx = float(m0.sum())
    I = []; gap = []
    for lam in lams:
        e = np.exp(lam*Fm); E = e.mean(axis=1); dl = ((Fm*e).mean(axis=1)/E).sum()   # Lambda'(lam) - mx
        Lm = lam*mx + np.log(E).sum(); y = mx + dl
        I.append(lam*y - Lm); gap.append(-dl)
    return mx, np.array(I), np.array(gap)
us = np.exp(np.linspace(math.log(1e-4), math.log(5.0), 400))
res = []
for i, row in enumerate(rows):
    n, c = T.coeffs(row, row['t0'], row['t1'])
    tri = float(np.abs(c).sum()); ray = tri - row['A']
    mx, I, gap = rate_curve(c, n)
    beta = B[i]['beta']; mass = B[i]['mass']; s = (1 - beta)*mass/4
    out = {'row': i, 'beta': beta, 'mass': mass, 's': s, 'ray': ray, 'mx': mx}
    for kap in (1.32, 7.0):
        logN = np.log(kap*X*(np.expm1(us)))
        g = np.interp(logN, I, gap, left=gap[0], right=gap[-1])   # I increasing in lam, gap decreasing
        tot = g + s*us; k = int(np.argmin(tot))
        out['cr_%g' % kap] = float(tot[k]); out['u_%g' % kap] = float(us[k]); out['logN_%g' % kap] = float(logN[k])
    res.append(out)
    if i % 6 == 0 or i == 66:
        print('%2d beta %.4f M %.4f s %.4f ray %.5f | cr*(k=1.32) %.5f at x/X = %.3f (logN %.1f) | cr*(k=7) %.5f at x/X = %.3f' % (i, beta, mass, s, ray, out['cr_1.32'], math.exp(out['u_1.32']), out['logN_1.32'], out['cr_7'], math.exp(out['u_7'])), flush=True)
json.dump(res, open('.scratch/ceilaudit5/selfcons.json', 'w'), indent=0)
sys.argv = [sys.argv[0], arc, '/dev/null']
src = open('Riemann/Deep/DBNHigh/chain.py').read(); src = src[:src.index('t0 = time.time()')]
exec(src)
E0 = run()[0]; print('base endpoint', repr(E0), flush=True)
for kap in ('1.32', '7'):
    dS = [max(0.0, r['cr_' + kap] - max(0.0, r['ray'])) for r in res]
    E = run(dS=dS)[0]; print('kappa', kap, 'mean new credit %.4f (min %.4f row %d); endpoint %r gain %.3e' % (sum(dS)/67, min(dS), dS.index(min(dS)), E, E0 - E), flush=True)
    E2 = run(dS=[0.5*x for x in dS])[0]; print('kappa', kap, 'half credit: endpoint %r gain %.3e' % (E2, E0 - E2), flush=True)
