#!/usr/bin/env python3
# ceilaudit5 ldp.py: how close can x in a window [X, X+W] bring row-0's head to its torus max?  (float model)
# phases theta_p = (x/2) log p, p <= 4096, ~ independent uniform over a long window (Kronecker-Weyl);
# sup over the window ~ level y* with P(Q >= y*) = 1/N_eff (Cramer large deviations of a sum of independent f_p(theta_p)),
# N_eff = W * (max frequency log(4096)/2)/pi.  Also: head mass by prime.
import sys, math, json
import numpy as np
sys.path.insert(0, 'Riemann/Deep/DBNHigh'); sys.path.insert(0, '.scratch/ceilaudit5')
import model as M, torus as T
rows = M.load_rows(sys.argv[1]); out = {}
for i in [int(x) for x in sys.argv[2].split(',')]:
    row = rows[i]; t = row['t0']
    n, c = T.coeffs(row, t, row['t1'])
    pp = {}
    for j, m in enumerate(n):
        f = T.FAC[m]
        if len(f) == 1:
            (p, k), = f.items(); pp.setdefault(p, {})[k] = c[j]
    th = np.linspace(0, 2*math.pi, 4097)[:-1]
    fs = {p: sum(a*np.cos(k*th) for k, a in ks.items()) for p, ks in pp.items()}
    mx = sum(float(f.max()) for f in fs.values())
    primes = sorted(fs)
    cum = np.cumsum([float(fs[p].max()) for p in primes])
    marks = {P: float(cum[[k for k, p in enumerate(primes) if p <= P][-1]]/mx) for P in (2, 3, 5, 7, 13, 31, 101, 1009)}
    def Lam(lam):
        s = 0.0; ds = 0.0
        for f in fs.values():
            m0 = float(f.max()); e = np.exp(lam*(f - m0)); E = float(e.mean())
            s += lam*m0 + math.log(E); ds += float((f*e).mean())/E
        return s, ds
    res = {}
    for logN in (20.0, 25.0, 27.4, 30.0, 35.0):
        # find y with I(y) = sup_lam (lam y - Lam(lam)) = logN ; I(y) = lam y - Lam(lam) at lam with Lam'(lam) = y
        lo, hi = 0.0, 2000.0
        for _ in range(80):
            lam = 0.5*(lo + hi); s, ds = Lam(lam); I = lam*ds - s
            if I < logN: lo = lam
            else: hi = lam
        s, ds = Lam(lam); res[logN] = (ds, mx - ds, lam)
    W = 6.6e12 - 6e12; Neff = W*(math.log(4096)/2)/math.pi
    print('row %d: torus max (prime powers) %.6f; mass fraction of primes <= P: %s' % (i, mx, {k: round(v, 4) for k, v in marks.items()}))
    print('   window [6e12, 6.6e12]: log N_eff = %.1f' % math.log(Neff))
    for logN, (y, gap, lam) in res.items():
        print('   log N = %.1f: sup level y* = %.5f, torus max - y* = %.5f (%.2f%% of max), lambda %.1f' % (logN, y, gap, 100*gap/mx, lam))
    out[i] = dict(mx=mx, marks=marks, res={str(k): v for k, v in res.items()})
json.dump(out, open(sys.argv[3], 'w'), indent=1)
