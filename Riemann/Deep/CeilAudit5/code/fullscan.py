#!/usr/bin/env python3
# ceilaudit5 fullscan.py: actual running sup of row i's FULL head polynomial (all 4095 labels n <= 4096, composites included)
# Q(x) = sum_n c_n cos((x/2) log n) over x in [X, X + W], grid step h, via type-1 NUFFT (finufft) in chunks of J samples;
# compared with the Cramer levels of the prime-power part (independent phases).
import sys, math, json, time
import numpy as np, mpmath as mp, finufft
sys.path.insert(0, 'Riemann/Deep/DBNHigh'); sys.path.insert(0, '.scratch/ceilaudit5')
import model as M, torus as T
mp.mp.dps = 40
row_i = int(sys.argv[2]); W = float(sys.argv[3]); h = float(sys.argv[4]); J = int(sys.argv[5])
rows = M.load_rows(sys.argv[1]); row = rows[row_i]
n, c = T.coeffs(row, row['t0'], row['t1'])
X = mp.mpf(6000000185827)
om = np.log(n.astype(float))/2
ph0 = np.array([float(mp.fmod(X*mp.log(int(k))/2, 2*mp.pi)) for k in n])
xt = h*om
tri = float(np.abs(c).sum())
best = -1e9; bestx = None; marks = [10**k for k in range(1, 12) if 10**k <= W]; mi = 0; res = []
chunk = J*h; m = 0; t0 = time.time()
while m*chunk < W:
    x0off = m*chunk
    phi = ph0 + x0off*om
    cp = c*np.exp(1j*(phi + (J//2)*xt))
    f = finufft.nufft1d1(xt, cp.astype(np.complex128), J, eps=1e-11, isign=1, nthreads=1)
    Q = f.real
    if m == 0:
        for j in (0, 12345, J - 1):
            direct = float(np.sum(c*np.cos(ph0 + j*h*om)))
            print('  check j=%d nufft %.10f direct %.10f' % (j, Q[j], direct), flush=True)
    while mi < len(marks) and marks[mi] <= x0off + chunk:
        k = int((marks[mi] - x0off)/h) + 1
        b = max(best, float(Q[:k].max()))
        res.append((marks[mi], b)); mi += 1
    jm = int(np.argmax(Q))
    if Q[jm] > best: best = float(Q[jm]); bestx = x0off + jm*h
    m += 1
print('row %d: triangle %.6f; scanned W = %.3g at step %g in %.0f s; overall sup %.6f at X + %.1f' % (row_i, tri, W, h, time.time() - t0, best, bestx), flush=True)
# Cramer levels of the prime-power part
G = 2048; th = np.arange(G)*2*math.pi/G
pp = {}
for j, mm in enumerate(n):
    fct = T.FAC[mm]
    if len(fct) == 1:
        (p, k), = fct.items(); pp.setdefault(p, {})[k] = c[j]
tabs = np.array([sum(a*np.cos(k*th) for k, a in ks.items()) for ks in pp.values()])
m0 = tabs.max(axis=1); Fm = tabs - m0[:, None]; mx = float(m0.sum())
lams = np.exp(np.linspace(math.log(0.5), math.log(8000.0), 260)); I = []; gp = []
for lam in lams:
    e = np.exp(lam*Fm); E = e.mean(axis=1); dl = ((Fm*e).mean(axis=1)/E).sum(); Lm = lam*mx + np.log(E).sum(); y = mx + dl
    I.append(lam*y - Lm); gp.append(-dl)
I = np.array(I); gp = np.array(gp)
out = dict(row=row_i, W=W, h=h, tri=tri, pp_torus_max=mx, sup=best, at=bestx, runs=[])
for w, b in res:
    preds = {kap: float(np.interp(math.log(kap*w), I, gp)) for kap in (1.32, 7.0)}
    out['runs'].append(dict(w=w, sup=b, gap=mx - b, ldp132=preds[1.32], ldp7=preds[7.0]))
    print('  w %.0e: sup %.5f  gap to pp torus max %.5f (%.2f%%) | LDP gap at kappa 1.32: %.5f, kappa 7: %.5f' % (w, b, mx - b, 100*(mx - b)/mx, preds[1.32], preds[7.0]), flush=True)
json.dump(out, open(sys.argv[6], 'w'), indent=1)
