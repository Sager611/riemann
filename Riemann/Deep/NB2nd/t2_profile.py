# nb2nd t2: (a) head-problem optimizer beta* = G_K^{-1} l_K vs L (c_opt - c_ramp) at finite N (small n);
# (b) endpoint profile of the optimum: w(y) = L*(P_opt - ramp) on squarefree n near N, y = log(N/n), binned;
# (c) head/tail decoupling at N <= 640: Delta(K, kap) = S(K u tail_kap) - S(K);
# (d) every N <= 640: L^2 d_N^2 - C L (double, nested Cholesky) and L^2 I_N - C L - J.
import math, numpy as np, mpmath as mp
O = '/Users/adriansager/playground/riemann/.scratch/nb2nd/'
Z = np.load(O + 'G640.npz'); G = Z['G']; b = Z['b']; mu = Z['mu'].astype(float); NMAX = len(b)
EG = 0.57721566490153286061; C = 2 + EG - math.log(4 * math.pi); J = 3.5626634565363677
mp.mp.dps = 20
a1 = float(mp.zeta(0, derivative=2) / mp.zeta(0) - (mp.zeta(0, derivative=1) / mp.zeta(0)) ** 2)
z = [0.0] + [float(mp.zeta(2 * k + 1)) for k in range(1, 200)]
def ell(n):
    s0 = -(math.log(2 * math.pi) * (1 + math.log(n) - EG) + a1) / n
    if n == 1: tail = sum((z[k] - 1) / (2 * k * (2 * k + 1)) for k in range(1, 199)) + (1 - math.log(2))
    else:
        x = 1.0 / n; tail = sum(z[k] * x ** (2 * k + 1) / (2 * k * (2 * k + 1)) for k in range(1, 60))
    return s0 - tail
ks = np.arange(1, NMAX + 1); lth = np.array([ell(int(n)) for n in ks])
def q(Gs, v):
    Lc = np.linalg.cholesky(Gs); y = np.linalg.solve(Lc, v); return float(y @ y)
# (a) head optimizer
bstar = np.linalg.solve(G, lth)   # K = 640, N-free: optimal B = sum beta_n n^{-s}, J(B) = J - Sigma
print('(a) head optimizer beta*_n (K=640):', np.round(bstar[:12], 4).tolist())
print('    sum_n beta*_n/n =', round(float(bstar @ (1.0 / ks)), 5), ';  Sigma_640 =', round(float(lth @ bstar), 6))
for N in [160, 320, 640]:
    L = math.log(N); c_r = -mu[:N] * (1 - np.log(ks[:N]) / L); c_o = np.linalg.solve(G[:N, :N], b[:N])
    # A = -sum c_k k^-s; A_opt = V_N + B/L with B = -L (c_o - c_r): beta_n = -L (c_o - c_r)_n
    beta = -L * (c_o - c_r)
    print(f'    N={N}: beta_n (opt, finite N) =', np.round(beta[:12], 4).tolist())
# (b) endpoint profile at N=640: P_opt(n) = -c_o/mu, ramp 1 - log n/L, w = L (P_opt - ramp) for squarefree n
for N in [320, 640]:
    L = math.log(N); c_o = np.linalg.solve(G[:N, :N], b[:N]); n = ks[:N]; sf = mu[:N] != 0
    P = np.where(sf, -c_o / np.where(sf, mu[:N], 1), np.nan); w = L * (P - (1 - np.log(n) / L)); y = np.log(N / n)
    edges = [0, 0.02, 0.05, 0.1, 0.2, 0.35, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0]
    out = []
    for i in range(len(edges) - 1):
        m = sf & (y >= edges[i]) & (y < edges[i + 1])
        if m.any(): out.append((edges[i], round(float(np.nanmean(w[m])), 3), round(float(np.nanstd(w[m])), 3), int(m.sum())))
    nsf = ~sf
    print(f'(b) N={N}: endpoint profile w(y)=L(P_opt-ramp) binned (y_lo, mean, sd, count):', out)
    print(f'    non-squarefree L*c_opt near the end (y<0.5): rms {L*np.sqrt(np.mean(c_o[nsf & (y<0.5)]**2)):.3f}; y>2: rms {L*np.sqrt(np.mean(c_o[nsf & (y>2)]**2)):.3f}')
# (c) decoupling
for N in [320, 640]:
    L = math.log(N); c = -mu[:N] * (1 - np.log(ks[:N]) / L); GN = G[:N, :N]; r = GN @ c - b[:N]
    S = L * L * q(GN, r)
    print(f'(c) N={N}: S_N = {S:.5f}')
    for K in [10, 30, 60, 100]:
        SK = L * L * q(GN[:K, :K], r[:K]); row = []
        for kap in [0.05, 0.1, 0.25, 0.5, 1.0, 1.5]:
            lo = int(math.ceil(N * math.exp(-kap)))
            if lo > K + 1:
                idx = np.concatenate([np.arange(K), np.arange(lo - 1, N)])
                row.append((kap, round(L * L * q(GN[np.ix_(idx, idx)], r[idx]) - SK, 4)))
        print(f'    K={K}: S(K)={SK:.4f}; Delta(K,kap)=S(K u [Ne^-kap,N]) - S(K):', row, f' ; S_N - S(K) = {S - SK:.4f}')
# (d) all N
Lc = np.linalg.cholesky(G); y = np.linalg.solve(Lc, b); d2 = 1 - np.cumsum(y * y)
Ns = np.arange(2, NMAX + 1); LL = np.log(Ns)
V1 = LL * LL * d2[Ns - 1] - C * LL
Iall = np.array([(lambda c: 1 - 2 * c @ b[:N] + c @ G[:N, :N] @ c)(-mu[:N] * (1 - np.log(ks[:N]) / math.log(N))) for N in Ns])
R = LL * LL * Iall - C * LL - J
np.savez(O + 't2_allN.npz', N=Ns, V1=V1, R=R)
for lo, hi in [(50, 100), (100, 200), (200, 400), (400, 640)]:
    m = (Ns >= lo) & (Ns <= hi)
    print(f'(d) N in [{lo},{hi}]: L^2 d^2 - CL mean {V1[m].mean():+.5f} sd {V1[m].std():.5f} | ramp R=L^2 I - CL - J mean {R[m].mean():+.5f} sd {R[m].std():.5f}')
