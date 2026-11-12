# Part A: pair-kernel 'graph' of zeros at normalized scale s: what a weak-regularity / dense-model step sees.
# For a simple point set X (density 1) and a dense model g (cell counts at resolution r=s), the pair kernel
# W_s = 1[|x-y|<=s] differs from the model's by (i) the atomic (diagonal) term: exactly N pairs (x,x) counted by g*g,
# relative size 1/(2s); (ii) cell-count fluctuations: sqrt(Sigma2(s))/s relative (Sigma2 = number variance).
# Also: the sparse cut deviation of the s-graph from its dense model, over unions of cells (exact max over S=T sets
# of cells via the sign trick is NP-hard; we report the l1 row-sum lower bound and the cell-diagonal deficit).
import json, sys, numpy as np
rng = np.random.default_rng(1)
import os
R = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'ExploreHelio', 'results') + '/'
def unfold_zeta(t):
    return t/(2*np.pi)*np.log(t/(2*np.pi*np.e)) + 7/8
def unfold_dh(t):  # DH: L(s,chi) mod 5 functional eq.: N(T) ~ T/2pi log(5T/2pi e)
    return t/(2*np.pi)*np.log(5*t/(2*np.pi*np.e))
z = np.load(R+'zeta_zeros.npy'); z = z[z > 100]
dh = np.array(json.load(open(R+'dh_zeros.json'))['online'])
sets = {}
sets['zeta (100<t<2e4)'] = unfold_zeta(z)
sets['DH on-line (5<t<420)'] = unfold_dh(dh)
# GUE bulk via Dumitriu-Edelman tridiagonal beta=2, unfolded with the semicircle
from scipy.linalg import eigh_tridiagonal
n = 20000
d = rng.normal(0, np.sqrt(2), n)
e = np.sqrt(rng.chisquare(2*np.arange(n-1, 0, -1)))
ev = eigh_tridiagonal(d/np.sqrt(2), e/np.sqrt(2), eigvals_only=True)
x = ev/(np.max(np.abs(ev))*1.0000001)  # semicircle on [-1,1] (radius = max |eigenvalue|)
F = lambda x: n*(0.5 + (x*np.sqrt(1-x*x) + np.arcsin(x))/np.pi)
u = F(x); u = u[(u > 0.25*n) & (u < 0.75*n)]
sets['GUE bulk (n=2e4)'] = u
sets['Poisson'] = np.sort(rng.uniform(0, 10000, 10000))
sets['lattice'] = np.arange(10000) + 0.5
print('name, N, mean gap')
for k, v in sets.items():
    print('  %-22s N=%6d mean gap %.4f' % (k, len(v), np.mean(np.diff(v))))
S = [0.5, 1, 2, 4, 8, 16]
def numvar(v, s, m=20000):
    v = np.sort(v); v = (v - v[0])/np.mean(np.diff(v))  # renormalize to mean gap 1
    L = v[-1]
    if L < 3*s: return np.nan
    x0 = rng.uniform(0, L - s, m)
    c = np.searchsorted(v, x0 + s) - np.searchsorted(v, x0)
    return np.var(c)
print('\nscale s | atomic 1/(2s) | fluctuation sqrt(Sigma2(s))/s for each set')
hdr = '%6s %8s ' % ('s', 'atomic') + ' '.join('%14s' % k[:14] for k in sets)
print(hdr)
for s in S:
    row = '%6.1f %8.3f ' % (s, 1/(2*s))
    for k, v in sets.items():
        nv = numvar(v, s)
        row += '%14.3f' % (np.sqrt(nv)/s)
    print(row)
print('\nSigma2(s) itself (GUE theory (1/pi^2)(log 2pi s + 1.5772 - pi^2/8) for s>=1):')
for s in S:
    th = (np.log(2*np.pi*s) + 0.5772156649 + 1 - np.pi**2/8)/np.pi**2
    print('  s=%5.1f  theory %.3f  ' % (s, th) + '  '.join('%s %.3f' % (k[:6], numvar(v, s)) for k, v in sets.items()))
