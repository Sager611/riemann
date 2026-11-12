# t10_robust.py (explore39): convergence of the flip ranges in n and of Barta-LP in the grid size.
import sys, numpy as np, math, json, time
sys.path.insert(0, '.scratch/explore39'); from xform import Form, barta_lp
t0 = time.time()
def galflip(c, sec, n):
    F = Form(c); Qm, p = F.Q(n, sec); Fm, _ = F.flip(n, sec, m2=600); e, V = np.linalg.eigh(Qm - Fm); return float(e[0]), F, V[:, 0]
for sec, lo, hi in (('even', 1.90, 1.95), ('odd', 2.04, 2.08)):
    for n in [14, 20, 28]:
        a, b = lo, hi
        for _ in range(16):
            m = (a + b)/2
            if galflip(m, sec, n)[0] > 0: a = m
            else: b = m
        print(json.dumps({'sec': sec, 'n': n, 'galerkin_crossing': round((a + b)/2, 5), 't': round(time.time() - t0, 1)}), flush=True)
for c, sec in ((1.9, 'even'), (2.05, 'odd')):
    for npts in [400, 800]:
        g, F, v = galflip(c, sec, 20)
        b, _, _ = barta_lp(F, 20, sec, flip=True, npts=npts, v0=v)
        print(json.dumps({'c': c, 'sec': sec, 'npts': npts, 'galerkin': g, 'barta_lp': b, 't': round(time.time() - t0, 1)}), flush=True)
