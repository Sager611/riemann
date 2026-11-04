# check_signs_grid.py (was t5): dense-grid sign check of H1, H2 built from the re-solved (near-exact) coefficients, nodes <= 240.
import sys, os, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from fastnp import *
import os as _os
_ROOT = _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..', '..'))
SRC = _os.environ.get('CIRC2D_SRC', _os.path.join(_ROOT, '.scratch', 'circ2d', 'src'))
SCR = _os.environ.get('CIRC2D_SCRATCH', _os.path.join(_ROOT, '.scratch', 'circ2d'))
ex = json.load(open(SCR + '/exact_coeffs.json'))
ns = ex['nodes']
which = sys.argv[1] if len(sys.argv) > 1 else 'exact'
L = {1: {0: FIXED[1][0]}, 2: {0: FIXED[2][0]}}
Cc = {1: FIXED[1][1], 2: FIXED[2][1]}
if which == 'exact':
    for k, n in enumerate(ns):
        L[1][n] = (ex['c1'][k], ex['d1'][k]); L[2][n] = (ex['c2'][k], ex['d2'][k])
else:  # the paper's finite table (zero beyond 100)
    tab = read_table(SRC + '/build/data/coefficients.tex')
    for n, (c1, d1, c2, d2) in tab.items():
        L[1][n] = (c1, d1); L[2][n] = (c2, d2)
Q = Quad(gl=48, sub=4)
dens = {i: Q.density(L[i]) for i in (1, 2)}
N0 = [0] + ns
N0a = np.array(N0, float)
sig = {1: -1.0, 2: 1.0}
SMAX = float(sys.argv[2]) if len(sys.argv) > 2 else 240.0
t0 = time.time()
# grid: uniform step 1e-3 plus geometric refinements near every node
grid = np.arange(0.0, SMAX + 1e-9, 1e-3)
near = []
for n in N0:
    for u in np.geomspace(2e-3, 0.05, 40):
        near += [n + u, n - u]
grid = np.unique(np.concatenate([grid, np.array([g for g in near if 0 <= g <= SMAX])]))
# nearest center and distance
idxn = np.searchsorted(N0a, grid)
lo = N0a[np.clip(idxn - 1, 0, len(N0a) - 1)]; hi = N0a[np.clip(idxn, 0, len(N0a) - 1)]
m = np.where(np.abs(grid - lo) <= np.abs(hi - grid), lo, hi)
m = np.where((grid == 41.5), 43.0, m)
u = grid - m
keep = np.abs(u) >= 2e-3
grid, m, u = grid[keep], m[keep], u[keep]
res = {}
for i in (1, 2):
    p = Pval(grid) * R_val(L[i], Cc[i], grid)
    K = Q.K(dens[3 - i], Cc[3 - i], grid, 0)
    H = p + K
    nu = np.where(m == 0, 0, np.where((m == 1) & (i == 1), 1, 2)).astype(float)
    dom = grid >= (1.0 if i == 1 else 0.0)
    q = sig[i] * H / np.abs(u) ** nu * np.where((nu == 1) & (u < 0), -1, 1)
    q = q[dom]; mm = m[dom]; gg = grid[dom]; HH = H[dom]
    res[i] = (q, mm, gg, HH)
    bad = sig[i] * HH < 0
    print('i=%d: grid pts %d, sign violations (sigma_i H_i<0): %d' % (i, len(gg), bad.sum()), flush=True)
    if bad.any():
        print('   first violations s=', gg[bad][:10], 'H=', HH[bad][:10])
print('grid time', round(time.time() - t0, 1), flush=True)
# minimum quotient per center group
groups = [('m=0', lambda x: x == 0, 0.735), ('m=1', lambda x: x == 1, 0.168), ('m=3,4', lambda x: (x == 3) | (x == 4), 0.017),
          ('4<m<=40', lambda x: (x > 4) & (x <= 40), 0.006), ('43<=m<=100', lambda x: (x >= 43) & (x <= 100), 0.00448),
          ('100<m', lambda x: x > 100, 0.00448)]
for i in (1, 2):
    q, mm, gg, HH = res[i]
    for name, f, claim in groups:
        sel = f(mm) & (gg <= SMAX - 6)
        if sel.any():
            k = np.argmin(q[sel])
            print('i=%d %-11s min quotient %.6g at s=%.4f   (paper lower bound %.5g)' % (i, name, q[sel][k], gg[sel][k], claim))
# node curvatures sigma_i H_i''(m)/2 (limit of the quotient) and H_1'(1)
print('node limits sigma_i*H_i^(nu)(m)/nu!:')
lim = {}
for i in (1, 2):
    Ms = np.array(N0, float)
    d0 = Q.p(dens[i], Cc[i], Ms, 0) + Q.K(dens[3 - i], Cc[3 - i], Ms, 0)
    d1 = Q.p(dens[i], Cc[i], Ms, 1) + Q.K(dens[3 - i], Cc[3 - i], Ms, 1)
    d2 = Q.p(dens[i], Cc[i], Ms, 2) + Q.K(dens[3 - i], Cc[3 - i], Ms, 2)
    lim[i] = (d0, d1, d2)
    vals = sig[i] * d2 / 2
    sel = Ms >= (3 if i == 1 else 1)
    print(' i=%d max|H(n)| over n>=1: %.3e  max|H\'(n)| n>=1 (excl. i=1,n=1): %.3e  min sigma*H\'\'/2 over nodes: %.6g at n=%d' % (
        i, np.abs(d0[1:]).max(), np.abs(np.where((Ms == 1) & (i == 1), 0, d1))[1:].max(), vals[sel].min(), int(Ms[sel][np.argmin(vals[sel])])))
    if i == 1:
        Q1 = node_QD(1)[0]
        print('   H1(1)=%.3e  H1\'(1)=%.12f  -Q1=%.12f   -H1\'(1) (quotient at 1+) ' % (d0[1], d1[1], -Q1))
    print('   H(0)=%.12f' % d0[0])

