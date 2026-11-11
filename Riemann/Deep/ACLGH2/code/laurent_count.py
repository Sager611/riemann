# laurent_count.py: counting function N(kappa) = #{n : kappa_n > kappa} of the Laurent (band-limited, unitary) model, from all
# eigenvalues of the symmetric Hankel matrix H (eigenvalues only, float32, in place). Same model as laurent.py.
# Usage: laurent_count.py primes tmax depth [t] [b1]
import os
for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'): os.environ[k] = '1'
import sys, json, math, time, numpy as np
import scipy.linalg as sla
from scipy.special import loggamma
t0 = time.time()
primes = [] if sys.argv[1] == '-' else [int(q) for q in sys.argv[1].split(',')]
tmax = float(sys.argv[2]); depth = float(sys.argv[3])
tsc = float(sys.argv[4]) if len(sys.argv) > 4 else 1.0
b1 = float(sys.argv[5]) if len(sys.argv) > 5 else 0.75
dx = math.pi/tmax; N = int(round(depth/dx))
def psi(tau):
    th = np.imag(loggamma(0.25 + 0.5j*tau)) - 0.5*tau*math.log(math.pi)
    out = -2*th + 2*tau*math.log(tsc)
    for p in primes:
        a = p**-0.5; L = math.log(p)
        out = out + 2*np.arctan2(a*np.sin(tau*L), 1 - a*np.cos(tau*L))
    return out
tau1 = b1*tmax; tau2 = 0.5*(1 + b1)*tmax
p1 = float(psi(np.array([tau1]))[0]); d1 = float((psi(np.array([tau1 + 1e-4])) - psi(np.array([tau1 - 1e-4])))[0]/2e-4)
K = round((p1 + min(d1, -0.5)*(tmax - tau1))/math.pi)
slope = (math.pi*K - p1)/(tmax - tau1)
def step(x):
    x = np.clip(x, 0, 1); a = np.exp(-1/np.maximum(x, 1e-300)); b = np.exp(-1/np.maximum(1 - x, 1e-300))
    return np.where(x <= 0, 0.0, np.where(x >= 1, 1.0, a/(a + b)))
def psit(tau):
    s = np.sign(tau); a = np.abs(tau)
    bl = step((a - tau1)/(tau2 - tau1))
    return s*((1 - bl)*psi(a) + bl*(p1 + slope*(a - tau1)))
Nth = 1 << int(math.ceil(math.log2(4*N + 4096)))
k = np.arange(Nth); th = 2*math.pi*k/Nth; th = np.where(th >= math.pi, th - 2*math.pi, th)
c = (np.fft.fft(np.exp(1j*psit(th/dx)))/Nth).real[:2*N + 2].copy()
H = np.empty((N, N), dtype=np.float32, order='F'); ar = np.arange(N)
for i0 in range(0, N, 128):
    i1 = min(N, i0 + 128)
    H[i0:i1] = c[1 + np.arange(i0, i1)[:, None] + ar[None, :]]
mu = sla.eigvalsh(H, overwrite_a=True, check_finite=False, driver='evr')
del H
S = np.sort(np.abs(mu))[::-1]
KS = [0.5, 0.3, 0.2, 0.15, 0.1, 0.07, 0.05, 0.04, 0.03, 0.02, 0.01]
print(json.dumps({'info': {'model': 'laurent_count', 'primes': primes, 't': tsc, 'tmax': tmax, 'depth': depth, 'N': N, 'sec': round(time.time() - t0, 1)},
                  'N': {str(x): int((S > x).sum()) for x in KS}, 'kappa_top': [float(x) for x in S[:2000]]}))
