# Theorem 4 test (smoothed mean of R = 7/4) and the Merlini-Sala-Sala conjecture |lambda_tiny(n)| <= a log n, a ~ 2
# (arXiv 1903.06233, p. 2-3); lambda_tiny = lambda_n - lambda_bar_n, lambda_bar = Li coefficients of pi^{-s/2}Gamma(1+s/2).
import math, numpy as np
def rd(f):
    d = {}
    for l in open(f):
        if l.startswith('#'): continue
        a, b = l.split()[:2]; d[int(a)] = float(b)
    return d
lam = rd('logs/lam_r999.txt'); bar = rd('logs/trend_r999.txt')
NM = max(lam); n = np.arange(1, NM + 1, dtype=float)
L = np.array([lam[k] for k in range(1, NM + 1)]); B = np.array([bar[k] for k in range(1, NM + 1)])
g = 0.57721566490153286; kap = 1 + math.log(2 * math.pi) - g
t = n / 2 * (np.log(n) - kap); R = L - t; tiny = L - B
for k in (10, 100, 1000, 100000):
    print(f'lambda_bar_n - t(n) - 3/4 at n={k}: {B[k-1]-t[k-1]-0.75:+.3e}   (-1/(24n) = {-1/(24*k):+.3e})')
def bump(u): return np.where((u > 0.5) & (u < 1), np.exp(-1 / ((u - 0.5) * (1 - u) * 16 + 1e-300)), 0.0)
print('\nsmoothed mean of R(n), weight phi(n/N), phi = bump on (1/2,1):  [Theorem 4: 7/4 + O(log N/N)]')
for N in (1000, 3000, 10000, 30000, 100000):
    w = bump(n / N); m = np.sum(w * R) / np.sum(w); mt = np.sum(w * tiny) / np.sum(w)
    print(f'  N={N:6d}: mean R = {m:.6f}  (7/4 = 1.75)   mean lambda_tiny = {mt:.6f}  (pred 1)')
print('sharp means: mean R n<=1e5 =', round(R.mean(), 4), ' mean lambda_tiny =', round(tiny.mean(), 4))
viol = np.nonzero(np.abs(tiny) > 2 * np.log(n))[0]
print(f'\n|lambda_tiny(n)| > 2 log n: first n = {int(viol[0])+1} (lambda_tiny = {tiny[viol[0]]:.6f}, 2 log n = {2*math.log(viol[0]+1):.6f}); count n<=1e5: {len(viol)}')
k = np.arange(30, NM + 1)
rat = np.abs(tiny[k - 1]) / np.log(k)
i = int(np.argmax(rat)) + 30
print(f'max_{{30<=n<=1e5}} |lambda_tiny(n)|/log n = {rat.max():.4f} at n = {i} (lambda_tiny = {tiny[i-1]:.6f})')
for a in (2, 3, 4):
    v = np.nonzero(np.abs(tiny[29:]) > a * np.log(k))[0]
    print(f'  a={a}: first violation n = {int(v[0])+30 if len(v) else None}, count {len(v)}')
