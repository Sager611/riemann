# ceilaudit6: precision check of d_N^2 from a C Gram file: (a) relative entrywise symmetric noise eps on G -> change of d2;
# (b) norm of the optimal coefficient vector x = G_N^{-1} b at chosen N. usage: dN_check.py mode N gramfile d2ref.npy eps [q S...]
import sys, math, resource, numpy as np
from scipy.linalg import cholesky, solve_triangular
from scipy.special import digamma
mode, N, gf, ref, eps = sys.argv[1], int(sys.argv[2]), sys.argv[3], sys.argv[4], float(sys.argv[5])
ks = np.arange(1, N + 1)
if mode == 'zeta':
    EG = 0.57721566490153286061; b = (np.log(ks) + 1 - EG) / ks
else:
    q = int(sys.argv[6]); S = [float(x) for x in sys.argv[7:7 + q]]
    beta = sum(S[r] * (digamma((r + 1) / q) - digamma(r / q)) / q for r in range(1, q)); b = beta / ks
d2ref = np.load(ref)[:N]
rng = np.random.default_rng(1)
G = np.fromfile(gf, dtype=np.float64).reshape(N, N)
# noise: G_ij *= (1 + eps*u_ij), u symmetric uniform[-1,1]; applied row-block-wise to save memory
for r0 in range(0, N, 256):
    r1 = min(N, r0 + 256)
    U = rng.uniform(-1, 1, size=(r1 - r0, N))
    G[r0:r1] *= (1 + eps * U)
G = 0.5 * (G + G.T) if N <= 1500 else G  # symmetrise (small N only; large N uses lower triangle only below)
A = G.T
Lc = cholesky(A, lower=True, overwrite_a=True, check_finite=False)
y = solve_triangular(Lc, b, lower=True, check_finite=False)
d2p = 1.0 - np.cumsum(y * y)
piv = np.abs(np.diag(Lc))
print(f'eps={eps:g}: max|d2p-d2| over N<=1000: {np.max(np.abs(d2p[:1000]-d2ref[:1000])):.3e}; over all: {np.max(np.abs(d2p-d2ref)):.3e} at N={int(np.argmax(np.abs(d2p-d2ref)))+1}')
for n in [1000, 1280, 2000, 2560, 3000, 4000]:
    if n <= N:
        x = solve_triangular(Lc[:n, :n], y[:n], lower=True, trans='T', check_finite=False)
        print(f'N={n}: d2ref={d2ref[n-1]:.10e} d2pert={d2p[n-1]:.10e} diff={d2p[n-1]-d2ref[n-1]:+.2e} |x|_2={np.linalg.norm(x):.3e} |x|_inf={np.max(np.abs(x)):.3e} minpiv^2={np.min(piv[:n])**2:.3e}')
print('maxrss MB', resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6)
