# explore31 s3b: head-only refinement theta in [1/5,1] at larger N (natural block from G2880, cross block on the fly).
# usage: uv run --no-project --with numpy --with scipy python3 s3b_big.py N Q   (N <= 2880; peak footprint ~200 MB at N = 2880)
import sys, math, time
import numpy as np
from scipy.linalg import lapack
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from memcap import check, mem
from s3_hybrid import vas
OUT = '/Users/adriansager/playground/riemann/.scratch/explore31/'
__import__('os').makedirs(OUT, exist_ok=True)
EG = 0.57721566490153286061
C = 2 + EG - math.log(4 * math.pi)
N, Q = int(sys.argv[1]), int(sys.argv[2])
K0 = 5
t0 = time.time()
head = [k for k in range(Q, Q * K0 + 1) if k % Q]  # non-natural head points
nh = len(head)
n = N + nh
A = np.empty((n, n))
A[:N, :N] = np.load(OUT + 'G2880.npy', mmap_mode='r')[:N, :N]  # <rho_{Q/(Qi)}, rho_{Q/(Qj)}> = Q G_{Qi,Qj} = G_{ij}
for a, k in enumerate(head):
    for b in range(a, nh):
        A[N + a, N + b] = A[N + b, N + a] = Q * vas(k, head[b])
    for j in range(1, N + 1):
        A[N + a, j - 1] = A[j - 1, N + a] = Q * vas(k, Q * j)
    if a % 16 == 0:
        check(250, f'head {a}'); print(f'row {a}/{nh} t={time.time()-t0:.0f}s', flush=True)
x = np.concatenate((np.arange(1, N + 1) * 1.0, np.array(head) / Q))  # dilation 1/x
t = (np.log(x) + 1 - EG) / x  # <chi, rho_{1/x}>
L = math.log(N)
def dist2_inplace(Mf, tt):  # Mf Fortran-ordered symmetric matrix, factorised in place (no copy)
    c, info = lapack.dpotrf(Mf, lower=1, overwrite_a=1, clean=0)
    assert info == 0, info
    y, info = lapack.dpotrs(c, tt, lower=1)
    return 1 - float(tt @ y)
W = np.asfortranarray(A[:N, :N])  # one copy (natural block)
check(280, 'natural copy')
v0 = dist2_inplace(W, t[:N])
del W
v = dist2_inplace(A.T, t)  # C-order A: A.T is the Fortran view of the same symmetric matrix
check(280, 'full')
print(f'N={N} Q={Q}: natural d2L={v0*L:.6f} C2={L*(v0*L-C):+.5f}; +head[1/5,1]: d2L={v*L:.6f} C2={L*(v*L-C):+.5f}; gain L^2(d2_nat-d2)={L*L*(v0-v):.5f}  t={time.time()-t0:.0f}s footprint/resident {mem()}')
