# explore61 s3: sampling defect of the integer dilations at single non-integer dilations theta = a/b, for every N <= NMAX.
#   pi_N(theta)^2 = dist(rho_theta, V_N)^2 (prefix sums of one Cholesky of G_NMAX: Cholesky is nested in N)
#   om(theta)^2   = dist(rho_theta, W1)^2  (N-independent floor, closed form; s1_split.w1_perp_gram)
#   pi_N^2 - om^2 = dist(P_W1 rho_theta, V_N)^2  (in-W1 part; -> 0 under RH since V_inf = W1 <=> RH)
# usage: uv run --no-project --with numpy --with scipy python3 s3_rate.py NMAX
import sys, os, math, time
import numpy as np
from scipy.linalg import cholesky, solve_triangular
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from memcap import check, mem
from s1_split import vas, w1_perp_gram, GF, NMAXF, EG
t0 = time.time()
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 2000
G = np.array(np.memmap(GF, dtype=np.float64, mode='r', shape=(NMAXF, NMAXF))[:NMAX, :NMAX])
Lc = cholesky(G, lower=True, overwrite_a=True, check_finite=False); del G
check(280, 'chol')
TH = [(2, 3), (3, 4), (3, 5), (4, 5), (5, 6), (4, 7), (5, 7), (5, 8), (2, 5), (3, 7), (3, 8), (2, 7), (2, 9), (3, 10)]
ks = np.arange(1, NMAX + 1)
y = solve_triangular(Lc, (np.log(ks) + 1 - EG) / ks, lower=True, check_finite=False)
d2 = 1 - np.cumsum(y * y)
Ns = [10, 40, 160, 640, 1280, 2000]
Ns = [n for n in Ns if n <= NMAX]
out = {}
print('theta  om^2   pi_N^2/om^2 - 1 at N =', Ns, '  [L*(pi^2-om^2)/om^2 at same N]', flush=True)
for a, b in TH:
    Q, k = a, b  # theta = Q/k
    cv = np.array([Q * vas(Q * j, k) for j in range(1, NMAX + 1)])
    z = solve_triangular(Lc, cv, lower=True, check_finite=False)
    nrm = Q * vas(k, k)
    pi2 = nrm - np.cumsum(z * z)
    Om, B, Cp = w1_perp_gram([k], Q, 40000)
    Om2, B2, _ = w1_perp_gram([k], Q, 20000)
    om2 = float((2 * Om - Om2)[0, 0] - (2 * B - B2)[0] ** 2 / Cp)
    ex = pi2 / om2 - 1
    out[f'{a}/{b}'] = dict(pi2=pi2, om2=om2, nrm=nrm)
    print(f'{a}/{b}  {om2:.6f} (|rho|^2 {nrm:.4f})  ' + ' '.join(f'{ex[n-1]:.2e}' for n in Ns) + '   [' + ' '.join(f'{ex[n-1]*math.log(n):.3f}' for n in Ns) + ']', flush=True)
np.savez('/Users/adriansager/playground/riemann/Riemann/Deep/ExploreSparse/logs/s3_rate.npz', d2=d2, **{k.replace('/', '_'): v['pi2'] for k, v in out.items()},
         om2=np.array([v['om2'] for v in out.values()]), th=np.array([a / b for a, b in TH]))
print(f'chi: d_N^2 log N at N =', Ns, ' '.join(f'{d2[n-1]*math.log(n):.4f}' for n in Ns))
print(f'done {time.time()-t0:.0f}s mem {mem()}')
