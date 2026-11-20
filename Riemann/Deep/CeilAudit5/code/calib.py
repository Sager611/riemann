#!/usr/bin/env python3
# ceilaudit5 calib.py: actual running sup of row i's head polynomial (primes <= PMAX, prime powers <= 4096) over x in [X0, X0 + w]
# vs the Cramer/LDP level of the same truncated polynomial; base phases (X0/2) log p mod 2pi in mpmath.
import sys, math, json, time
import numpy as np, mpmath as mp
sys.path.insert(0, 'Riemann/Deep/DBNHigh'); sys.path.insert(0, '.scratch/ceilaudit5')
import model as M, torus as T
mp.mp.dps = 40
row_i = int(sys.argv[2]); PMAX = int(sys.argv[3]); X0s = [mp.mpf(s) for s in sys.argv[4].split(',')]; WMAX = float(sys.argv[5]); STEP = 0.03
rows = M.load_rows(sys.argv[1]); row = rows[row_i]
n, c = T.coeffs(row, row['t0'], row['t1'])
pp = {}
for j, m in enumerate(n):
    f = T.FAC[m]
    if len(f) == 1:
        (p, k), = f.items()
        if p <= PMAX: pp.setdefault(p, {})[k] = c[j]
G = 8192
th = np.arange(G)*2*math.pi/G
tab = {p: sum(a*np.cos(k*th) for k, a in ks.items()) for p, ks in pp.items()}
mx = sum(float(t.max()) for t in tab.values())
def Lam(lam):
    s = 0.0; ds = 0.0
    for f in tab.values():
        m0 = float(f.max()); e = np.exp(lam*(f - m0)); E = float(e.mean()); s += lam*m0 + math.log(E); ds += float((f*e).mean())/E
    return s, ds
def level(logN):
    lo, hi = 0.0, 3000.0
    for _ in range(70):
        lam = 0.5*(lo + hi); s, ds = Lam(lam); I = lam*ds - s
        if I < logN: lo = lam
        else: hi = lam
    return Lam(lam)[1]
out = dict(row=row_i, PMAX=PMAX, torus_max=mx, runs=[])
print('row %d, primes <= %d (%d primes): torus max %.6f' % (row_i, PMAX, len(tab), mx), flush=True)
for X0 in X0s:
    ph0 = {p: float(mp.fmod(X0*mp.log(p)/2, 2*mp.pi)) for p in tab}
    lp = {p: math.log(p)/2 for p in tab}
    best = -1e9; marks = [10**k for k in range(1, 9) if 10**k <= WMAX] + ([WMAX] if WMAX not in [10**k for k in range(1, 9)] else [])
    res = []; y0 = 0.0; CH = 1_000_000; t0 = time.time(); mi = 0
    while y0 < WMAX:
        y = y0 + STEP*np.arange(CH)
        Q = np.zeros(CH)
        for p, tb in tab.items():
            idx = (np.mod(ph0[p] + lp[p]*y, 2*math.pi)*(G/(2*math.pi))).astype(np.int64) % G
            Q += tb[idx]
        # running max checkpoints
        while mi < len(marks) and marks[mi] <= y[-1]:
            k = int((marks[mi] - y0)/STEP) + 1
            b = max(best, float(Q[:k].max())) if k > 0 else best
            res.append((marks[mi], b)); mi += 1
        best = max(best, float(Q.max())); y0 += CH*STEP
    for w, b in res:
        out['runs'].append(dict(X0=str(X0), w=w, sup=b, gap=mx - b))
        print('  X0 %s w %.0e: sup %.5f gap %.5f (%.2f%%)' % (mp.nstr(X0, 15), w, b, mx - b, 100*(mx - b)/mx), flush=True)
    print('  (%.0f s)' % (time.time() - t0), flush=True)
out['ldp'] = {}
for logN in (5, 8, 10, 12, 14, 16, 18, 20, 23, 27.4):
    y = level(logN); out['ldp'][logN] = (y, mx - y)
    print('  LDP logN %.1f: level %.5f gap %.5f (%.2f%%)' % (logN, y, mx - y, 100*(mx - y)/mx), flush=True)
json.dump(out, open(sys.argv[6], 'w'), indent=1)
