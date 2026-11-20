# ceilaudit7 D: is the hairpin weight X^{2d} of Theorem U' Step 6e' attained? Davenport-Heilbronn F = sum a_n n^{-s} (a periodic mod 5),
# truncated inverse ("mollifier") M_X(s) = sum_{n<X} b_n (1 - log n/log X)^3 n^{-s}, b = Dirichlet inverse of a. On the hairpin of an off-line zero
# rho = 1/2 + d + i g: s+ = 1/2 + y + i g, s- = 1/2 - y + i g (0 < y < d); Lemma Z: |M(s+)M(s-)| <~ X^{2d}. Fit the slope in log X.
import numpy as np
k = (np.sqrt(10-2*np.sqrt(5))-2)/(np.sqrt(5)-1)
Xmax = 200000
n = np.arange(Xmax+1); a = np.zeros(Xmax+1); r = n % 5
a[r == 1] = 1; a[r == 2] = k; a[r == 3] = -k; a[r == 4] = -1
b = np.zeros(Xmax+1); b[1] = 1.0
for m in range(1, Xmax//2+1):
    if b[m] != 0.0:
        d = np.arange(2, Xmax//m+1); b[d*m] -= a[d]*b[m]
print(f'kappa={k:.6f}; max|b_n| n<=1e3: {np.abs(b[1:1001]).max():.3g}, n<=Xmax: {np.abs(b[1:]).max():.3g}')
def M(s, X):
    nn = np.arange(1, int(X)); w = (1-np.log(nn)/np.log(X))**3
    return np.sum(b[1:int(X)]*w*np.exp(-s*np.log(nn)))
zs = [(0.8085171824566374, 85.6993484853776), (0.6508300806097371, 114.16334273075698), (0.574356050450806, 166.47930591316816), (0.5, 94.5)]
Xs = [1e3, 3e3, 1e4, 3e4, 1e5, 2e5]
for (beta, g) in zs:
    dd = beta-0.5
    for frac in (0.5, 0.9):
        y = frac*dd if dd > 0 else 0.1
        v = [abs(M(0.5+y+1j*g, X)*M(0.5-y+1j*g, X)) for X in Xs]
        sl = np.polyfit(np.log(Xs[2:]), np.log(v[2:]), 1)[0]
        print(f'zero beta={beta:.4f} g={g:.3f} d={dd:.4f} y={y:.4f}: |M+M-| = {[f"{x:.3g}" for x in v]}; slope (X>=1e4) = {sl:.3f}; 2d = {2*dd:.3f}', flush=True)
