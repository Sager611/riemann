# explore31 s1: q-refined dilation designs D_{q,N} = {q/k : q <= k <= qN} (all theta in [1/N,1] with numerator q),
# nested along divisibility, D_{1,N} = Baez-Duarte's {1/k}.  d_{q,N}^2 = 1 - q t^T G_sub^{-1} t, G_sub = Vasyunin G[q..qN],
# t_k = (log(k/q) + 1 - gamma)/k.  As q -> inf: delta(log N)^2 = distance with all dilations in [1/N, 1] (continuous relaxation).
# usage: uv run --no-project --with numpy --with scipy python3 s1_qdesign.py [N ...]   (needs .scratch/explore31/G2880.npy)
# memory: G (66 MB) + one preallocated work matrix (66 MB); LAPACK in place (F-order view of the C-order symmetric block).
import sys, math, time
import numpy as np
from scipy.linalg import lapack
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from memcap import check, mem
OUT = '/Users/adriansager/playground/riemann/.scratch/explore31/'
__import__('os').makedirs(OUT, exist_ok=True)
EG = 0.57721566490153286061
C = 2 + EG - math.log(4 * math.pi)
G = np.load(OUT + 'G2880.npy')
MM = G.shape[0]
W = np.empty(MM * MM)
def d2(N, q):
    M = q * (N - 1) + 1
    ks = np.arange(q, q * N + 1)
    A = W[:M * M].reshape(M, M)
    A[...] = G[q - 1:q * N, q - 1:q * N]
    t = (np.log(ks / q) + 1 - EG) / ks
    c, info = lapack.dpotrf(A.T, lower=1, overwrite_a=1, clean=0)
    assert info == 0, info
    y, info = lapack.dpotrs(c, t, lower=1)
    check(250, f'N={N} q={q}')
    return 1 - q * float(t @ y)
QS = [1, 2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128, 192, 256, 288, 384, 512, 576, 768, 960, 1024, 1440]
NS = [int(x) for x in sys.argv[1:]] or [2, 3, 5, 10, 20, 40, 80, 160, 320, 640, 1440, 2880]
t0 = time.time()
rows = []
for N in NS:
    L = math.log(N)
    for q in QS:
        if q * N > MM:
            break
        v = d2(N, q)
        rows.append((N, q, v))
        print(f'N={N:5d} q={q:5d} M={q*(N-1)+1:5d} d2={v:.12f} d2*L={v*L:.8f} L(d2L-C)={L*(v*L-C):+.6f}  t={time.time()-t0:.0f}s', flush=True)
np.save(OUT + 's1_res_' + '_'.join(map(str, NS)) + '.npy', np.array(rows))
print('footprint/resident MB', mem())
