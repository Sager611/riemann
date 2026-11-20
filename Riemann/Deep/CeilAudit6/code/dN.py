# ceilaudit6: d_N^2 for every N <= NMAX from a C Gram file via one in-place Cholesky. usage: dN.py mode N gramfile out.npy [q S...]
import sys, math, resource, numpy as np
from scipy.linalg import cholesky, solve_triangular
from scipy.special import digamma
mode, N, gf, out = sys.argv[1], int(sys.argv[2]), sys.argv[3], sys.argv[4]
G = np.fromfile(gf, dtype=np.float64).reshape(N, N)
ks = np.arange(1, N + 1)
if mode == 'zeta':
    EG = 0.57721566490153286061; b = (np.log(ks) + 1 - EG) / ks
else:
    q = int(sys.argv[5]); S = [float(x) for x in sys.argv[6:6 + q]]
    beta = sum(S[r] * (digamma((r + 1) / q) - digamma(r / q)) / q for r in range(1, q)); b = beta / ks
print('sym err', float(np.max(np.abs(G[:50, :50] - G[:50, :50].T))), 'G11', G[0, 0], flush=True)
A = G.T  # F-contiguous view, symmetric
Lc = cholesky(A, lower=True, overwrite_a=True, check_finite=False)
y = solve_triangular(Lc, b, lower=True, check_finite=False)
d2 = 1.0 - np.cumsum(y * y)
np.save(out, d2)
for n in [10, 40, 100, 400, 600, 1000, 1280, 2000, 2560, 3000, 4000]:
    if n <= N: print(f'N={n:5d} d2={d2[n-1]:.12e} d2L={d2[n-1]*math.log(n):.8f}')
print('maxrss MB', resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6)
