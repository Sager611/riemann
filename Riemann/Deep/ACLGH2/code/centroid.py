# centroid.py: Mellin-frequency location of the large-angle modes of the Laurent model: top kE eigenpairs at each end of H (float32),
# for each mode with kappa > kmin the |tau|-centroid and the |tau|-quartiles of |M e_n|^2 (FFT of the input vector).
# Usage: centroid.py primes tmax depth kE kmin
import os
for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'): os.environ[k] = '1'
import sys, json, math, numpy as np, scipy.linalg as sla
sys.argv = [sys.argv[0], sys.argv[1], sys.argv[2], sys.argv[3]] + sys.argv[4:]
primes = [] if sys.argv[1] == '-' else [int(q) for q in sys.argv[1].split(',')]
tmax = float(sys.argv[2]); depth = float(sys.argv[3]); kE = int(sys.argv[4]); kmin = float(sys.argv[5])
src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'laurent_count.py')).read()
src = src.split("H = np.empty((N, N), dtype=np.float32, order='F')")[0]          # reuse the symbol / coefficient code
sys.argv = [sys.argv[0], sys.argv[1], sys.argv[2], sys.argv[3]]
exec(src)
def build():
    H = np.empty((N, N), dtype=np.float32, order='F'); ar = np.arange(N)
    for i0 in range(0, N, 128):
        i1 = min(N, i0 + 128); H[i0:i1] = c[1 + np.arange(i0, i1)[:, None] + ar[None, :]]
    return H
H = build(); mu1, v1 = sla.eigh(H, subset_by_index=[N - kE, N - 1], driver='evr', overwrite_a=True, check_finite=False); del H
H = build(); mu2, v2 = sla.eigh(H, subset_by_index=[0, kE - 1], driver='evr', overwrite_a=True, check_finite=False); del H
mu = np.concatenate([mu1, mu2]); V = np.concatenate([v1, v2], axis=1).astype(np.float64)
o = np.argsort(-np.abs(mu)); mu = mu[o]; V = V[:, o]
nf = 1 << int(math.ceil(math.log2(2*N)))
F = np.abs(np.fft.rfft(V, n=nf, axis=0))**2                     # |M e(i tau)|^2 up to phase, tau_k = 2 pi k/(nf dx)
tau = 2*math.pi*np.arange(F.shape[0])/(nf*dx)
out = []
for j in range(mu.size):
    if abs(mu[j]) < kmin: break
    w = F[:, j]/F[:, j].sum(); cw = np.cumsum(w)
    out.append([round(float(abs(mu[j])), 4), round(float((tau*w).sum()), 1), round(float(tau[np.searchsorted(cw, 0.25)]), 1),
                round(float(tau[np.searchsorted(cw, 0.5)]), 1), round(float(tau[np.searchsorted(cw, 0.75)]), 1)])
print(json.dumps({'primes': primes, 'tmax': tmax, 'tau1': b1*tmax, 'N': N, 'modes_kappa_centroid_q25_q50_q75': out}))
