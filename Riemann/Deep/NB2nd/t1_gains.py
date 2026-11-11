# nb2nd t1: decomposition of the optimal gain S_N = L^2 (I_N - d_N^2) = L^2 r^T G^{-1} r (r = G c_ramp - b, ramp c_k = -mu(k)(1 - log k/L))
# by the support of the correction; theoretical l_n = <zeta'/zeta, zeta n^{-s}> vs L r_n; Sigma_K = l_K^T G_K^{-1} l_K.
# usage: uv run --no-project --with numpy --with mpmath python3 t1_gains.py   (needs .scratch/nb2nd/G640.npz from t0_gram.py)
import math, time, numpy as np, mpmath as mp
t0 = time.time()
O = '/Users/adriansager/playground/riemann/.scratch/nb2nd/'
Z = np.load(O + 'G640.npz'); G = Z['G']; b = Z['b']; mu = Z['mu'].astype(float); NMAX = len(b)
EG = 0.57721566490153286061; C = 2 + EG - math.log(4 * math.pi); J = 3.5626634565363677
mp.mp.dps = 30
a1 = float(mp.zeta(0, derivative=2) / mp.zeta(0) - (mp.zeta(0, derivative=1) / mp.zeta(0)) ** 2)
z = [0.0] + [float(mp.zeta(2 * k + 1)) for k in range(1, 200)]
def ell(n):  # residues at s = 0 and s = -2k of (zeta'/zeta)(s) zeta(1-s) n^(s-1)/(s(1-s))
    s0 = -(math.log(2 * math.pi) * (1 + math.log(n) - EG) + a1) / n
    if n == 1:
        tail = sum((z[k] - 1) / (2 * k * (2 * k + 1)) for k in range(1, 199)) + (1 - math.log(2))
    else:
        x = 1.0 / n; tail = sum(z[k] * x ** (2 * k + 1) / (2 * k * (2 * k + 1)) for k in range(1, 60))
    return s0 - tail
ks = np.arange(1, NMAX + 1)
lth = np.array([ell(int(n)) for n in ks])
print(f'alpha1 = {a1:.15f};  l_1..l_6 (theory) =', np.round(lth[:6], 6).tolist())
# formal check: sum mu(n) log n l_n (partial sums, conditionally -> J)
ps = np.cumsum(mu * np.log(ks) * lth)
print('partial sums sum_{n<=K} mu(n) log n l_n  (formal -> J = %.6f):' % J, [(K, round(ps[K - 1], 4)) for K in (10, 100, 300, 640)])
def solve_q(Gs, v):
    Lc = np.linalg.cholesky(Gs); y = np.linalg.solve(Lc, v); return float(y @ y)
print('Sigma_K = l_K^T G_K^{-1} l_K (theory l, N-free):')
for K in [1, 2, 3, 5, 10, 20, 40, 80, 160, 320, 640]:
    print(f'  K={K:4d}  Sigma_K = {solve_q(G[:K,:K], lth[:K]):.6f}', flush=True)
for N in [40, 80, 160, 320, 640]:
    L = math.log(N); c = -mu[:N] * (1 - np.log(ks[:N]) / L)
    GN = G[:N, :N]; bN = b[:N]
    I = 1 - 2 * c @ bN + c @ GN @ c; r = GN @ c - bN
    d2 = 1 - solve_q(GN, bN); S = L * L * (I - d2)
    print(f'N={N}: L^2 I_N - C L = {L*L*I - C*L:.5f} (J = {J:.5f}); L^2 d^2 - C L = {L*L*d2 - C*L:+.5f}; S_N = L^2(I-d^2) = {S:.5f}; check L^2 r^T G^-1 r = {L*L*solve_q(GN, r):.5f}')
    print('   L r_n vs l_n (n=1..6):', [(round(L * r[n - 1], 4), round(lth[n - 1], 4)) for n in range(1, 7)])
    out = []
    for K in [1, 2, 5, 10, 20, 50, 100, 200]:
        if K <= N: out.append((K, round(L * L * solve_q(GN[:K, :K], r[:K]), 4)))
    print('   small-n gain S(K) (support 1..K):', out)
    out = []
    for kap in [0.1, 0.25, 0.5, 1.0, 1.5, 2.0, 3.0]:
        lo = int(math.ceil(N * math.exp(-kap)))
        idx = np.arange(lo - 1, N); out.append((kap, len(idx), round(L * L * solve_q(GN[np.ix_(idx, idx)], r[idx]), 4)))
    print('   endpoint gain (support [N e^-kap, N]): (kap, size, gain):', out)
    out = []
    for K, kap in [(20, 1.0), (50, 1.0), (50, 2.0), (100, 1.0), (100, 2.0)]:
        lo = int(math.ceil(N * math.exp(-kap)))
        if K < lo:
            idx = np.concatenate([np.arange(K), np.arange(lo - 1, N)])
            out.append((K, kap, round(L * L * solve_q(GN[np.ix_(idx, idx)], r[idx]), 4)))
    print('   small+end gain (K, kap, gain):', out, flush=True)
print(f'done {time.time()-t0:.1f}s')
