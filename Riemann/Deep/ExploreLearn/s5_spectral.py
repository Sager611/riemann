# explore31 s5: "spectral bias" picture of the Baez-Duarte fit.  Eigen-decomposition of the natural Gram G_M (M = 1440):
# lambda_j (descending), target coefficients theta_j^2 = (u_j . b)^2 / lambda_j (explained variance per eigenmode), sum = 1 - d_M^2.
# Spectral-truncation learning curve E_spec(J) = 1 - sum_{j<=J} theta_j^2 (top-J eigenmodes, the KRR/PCA order) vs d_N^2 (first N dilations).
# usage: uv run --no-project --with numpy --with scipy python3 s5_spectral.py
import sys, math
import numpy as np
from scipy.linalg import eigh
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from memcap import check, mem
OUT = '/Users/adriansager/playground/riemann/.scratch/explore31/'
__import__('os').makedirs(OUT, exist_ok=True)
EG = 0.57721566490153286061
C = 2 + EG - math.log(4 * math.pi)
M = 1440
G = np.array(np.load(OUT + 'G2880.npy', mmap_mode='r')[:M, :M])
ks = np.arange(1, M + 1); b = (np.log(ks) + 1 - EG) / ks
lam, U = eigh(G, check_finite=False)
lam = lam[::-1]; U = U[:, ::-1]
th2 = (U.T @ b) ** 2 / lam
check(250, 'eigh')
print(f'M={M}: lambda_max {lam[0]:.4f}, lambda_min {lam[-1]:.3e}, cond {lam[0]/lam[-1]:.3e}; 1 - sum theta^2 = {1-th2.sum():.8f} (d_M^2)')
for j in [1, 2, 3, 5, 10, 20, 50, 100, 200, 500, 1000, 1440]:
    print(f'  j={j:5d} lambda_j={lam[j-1]:.4e}  j^2 lambda_j={j*j*lam[j-1]:.4f}  theta_j^2={th2[j-1]:.3e}  E_spec(j)={1-th2[:j].sum():.6f}')
cs = 1 - np.cumsum(th2)
print('compare first-N dilations d_N^2 with top-N eigenmodes E_spec(N):')
for N in [2, 5, 10, 20, 50, 100, 200, 500, 1000, 1440]:
    bN = b[:N]; dN = 1 - float(bN @ np.linalg.solve(G[:N, :N], bN))
    print(f'  N={N:5d}  d_N^2={dN:.6f}  E_spec(N)={cs[N-1]:.6f}  ratio {cs[N-1]/dN:.2f}')
# where does the explained variance sit?  mass of theta^2 per eigenvalue decade
edges = 10.0 ** np.arange(1, -12, -1)
print('theta^2 mass by eigenvalue decade:')
for hi, lo in zip(edges[:-1], edges[1:]):
    m = (lam <= hi) & (lam > lo)
    if m.any():
        print(f'  ({lo:.0e},{hi:.0e}]: {m.sum():5d} modes, theta^2 mass {th2[m].sum():.6f}')
print('footprint/resident MB', mem())
