import json, math, glob, re
import numpy as np
# moved-pair thresholds
rows = []
for f in sorted(glob.glob('logs/mp_k*_d*.log')):
    if 'N16' in f: continue
    L = [json.loads(l) for l in open(f) if l.strip().startswith('{')]
    head = L[0]; th = [x for x in L if 'threshold_c' in x]
    if not th: continue
    lo, hi = th[0]['threshold_c']
    rows.append((head['k'], head['m'], head['delta'], (lo+hi)/2, (hi-lo)/2))
rows.sort()
print('k m delta c* +- | m/(2pi) m/(4pi)')
for r in rows: print(r[0], round(r[1], 3), r[2], round(r[3], 3), round(r[4], 3), '|', round(r[1]/(2*math.pi), 3), round(r[1]/(4*math.pi), 3))
for d in (0.3, 0.03):
    pts = [(r[1], r[3]) for r in rows if r[2] == d]
    m = np.array([p[0] for p in pts]); c = np.array([p[1] for p in pts])
    A = np.vstack([m, np.ones_like(m)]).T
    (a, b), res, *_ = np.linalg.lstsq(A, c, rcond=None)
    print('delta', d, 'fit c* = %.5f m + %.3f  (1/slope = %.3f = %.4f pi), max resid %.3f' % (a, b, 1/a, 1/a/math.pi, np.max(np.abs(A@np.array([a, b]) - c))))
# DH prediction from the delta = 0.3 fit, slope scaled by conductor 5
pts = [(r[1], r[3]) for r in rows if r[2] == 0.3]
m = np.array([p[0] for p in pts]); c = np.array([p[1] for p in pts])
a, b = np.polyfit(m, c, 1)
print('DH pred (slope x5): %.2f ; (all x5): %.2f ; actual (30.74404, 30.74497] (Elston), ours (30, 31)' % (5*a*85.699348 + b, 5*(a*85.699348 + b)))
print('DH 2nd pair pred (114.163, delta .15, slope x5): %.2f ; ours (42, 44)' % (5*a*114.163343 + b))
# zeta ground state law from WeilSimple U(c)
U = {2:1.335e-3,3:5.967e-8,4:9.637e-13,5:1.222e-17,6:9.575e-23,7:1.006e-27,8:5.211e-33,9:3.229e-38,10:1.944e-43,11:1.371e-48,12:6.271e-54,13:3.721e-59,14:1.642e-64,15:8.446e-70,16:3.743e-75,17:2.095e-80,18:6.775e-86,19:4.223e-91,20:1.657e-96}
cs = np.array(sorted(U)); y = np.array([math.log(U[k]) + 4*math.pi*k for k in cs])
for lo in (2, 5):
    sel = cs >= lo
    p, q = np.polyfit(np.log(cs[sel]), y[sel], 1)
    print('ln U + 4 pi c = %.3f ln c + %.3f (c >= %d), max resid %.3f' % (p, q, lo, np.max(np.abs(p*np.log(cs[sel]) + q - y[sel]))))
# free slope fit
sel = cs >= 5
X = np.vstack([cs[sel], np.log(cs[sel]), np.ones(sel.sum())]).T
co, *_ = np.linalg.lstsq(X, np.log([U[k] for k in cs[sel]]), rcond=None)
print('free fit ln U = %.4f c + %.3f ln c + %.3f   (4 pi = %.4f)' % (co[0], co[1], co[2], 4*math.pi))
# DH ground state slope
dh = {16:4.4152e-14, 24:2.80248e-22, 28:7.80193e-27, 30:9.56938e-29}
ks = sorted(dh)
for i in range(len(ks)-1):
    print('DH slope %d->%d: %.3f (4pi/5 = %.3f)' % (ks[i], ks[i+1], (math.log(dh[ks[i+1]]) - math.log(dh[ks[i]]))/(ks[i+1]-ks[i]), 4*math.pi/5))
