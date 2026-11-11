# nbrate: second-order analysis (exploratory, double precision except J which is mpmath).
import numpy as np, math, mpmath as mp
O = '/Users/adriansager/playground/riemann/.scratch/nbrate/'
G = np.load(O + 'G_double.npy'); b = np.load(O + 'b_double.npy'); NMAX = len(b)
EG = 0.57721566490153286061; C = 2 + EG - math.log(4 * math.pi)
# (1) J = (zeta'/zeta)'(0) + (1+gamma) log 2pi - sum_{n>=1} (zeta'/zeta)(2n+1)/(2n(2n+1))
mp.mp.dps = 30
z0, z1, z2 = mp.zeta(0), mp.zeta(0, derivative=1), mp.zeta(0, derivative=2)
a1 = z2 / z0 - (z1 / z0) ** 2
lz = lambda s: mp.zeta(s, derivative=1) / mp.zeta(s)
tail = mp.nsum(lambda n: lz(2 * n + 1) / (2 * n * (2 * n + 1)), [1, mp.inf])
J = a1 + (1 + mp.euler) * mp.log(2 * mp.pi) - tail
print('alpha1 = (zeta\'/zeta)\'(0) =', mp.nstr(a1, 20), ' (1+g)log2pi =', mp.nstr((1 + mp.euler) * mp.log(2 * mp.pi), 20), ' sum =', mp.nstr(tail, 20))
print('J =', mp.nstr(J, 20))
# (2) optimal coefficients at N
mu = np.ones(NMAX + 1, dtype=np.int64); mu[0] = 0
for p in range(2, NMAX + 1):
    if all(p % q for q in range(2, int(p ** 0.5) + 1)):
        mu[p::p] *= -1; mu[p * p::p * p] = 0
ks = np.arange(1, NMAX + 1)
for N in [160, 1280]:
    L = math.log(N)
    c = np.linalg.solve(G[:N, :N], b[:N])
    sf = mu[1:N + 1] != 0
    prof = np.where(sf, -c / np.where(sf, mu[1:N + 1], 1), np.nan)
    x = np.log(ks[:N]) / L
    print(f'N={N}: c_1..c_12 =', np.round(c[:12], 4).tolist())
    print(f'  profile P_opt(x) = -c_k/mu(k) at small k (k, x, P, 1-x):', [(int(k), round(x[k-1], 3), round(prof[k-1], 4), round(1 - x[k-1], 4)) for k in [1, 2, 3, 5, 6, 7, 10, 11, 13, 30]])
    nsf = ~sf
    print(f'  max |c_k| over non-squarefree k: {np.abs(c[nsf]).max():.4f}; rms {np.sqrt(np.mean(c[nsf]**2)):.4f}; squarefree rms deviation from 1-x: {np.sqrt(np.nanmean((prof - (1 - x))**2)):.4f}')
    # binned profile
    edges = np.linspace(0, 1, 11)
    bins = []
    for i in range(10):
        m = sf & (x >= edges[i]) & (x < edges[i + 1])
        if m.any(): bins.append((round(edges[i], 1), round(float(np.mean(prof[m])), 4), round(float(np.mean(1 - x[m])), 4), int(m.sum())))
    print('  binned profile (x_lo, mean P_opt, mean 1-x, count):', bins)
# (3) profile families: c_k = -mu(k) P(x_k), P polynomial deg <= d (free), plus free c_n for n <= m
def fam(N, d, m):
    L = math.log(N); x = np.log(ks[:N]) / L
    cols = [(-mu[1:N + 1] * x ** j).astype(float) for j in range(d + 1)]
    for n in range(1, m + 1):
        e = np.zeros(N); e[n - 1] = 1; cols.append(e)
    B = np.array(cols).T
    M = B.T @ G[:N, :N] @ B; be = B.T @ b[:N]
    return 1 - be @ np.linalg.solve(M, be)
print('N, d2_opt, BCF(1-x), fam(d=3,m=0), fam(d=6,m=0), fam(d=6,m=6), fam(d=10,m=12); each as L^2 d2 - C L')
d1 = np.load(O + 'd1_double.npy'); bcf = np.load(O + 'bcf_rig_mid.npy')
for N in [40, 80, 160, 320, 640, 1280]:
    L = math.log(N); f = lambda q: L * L * q - C * L
    vals = [d1[N - 1], bcf[N], fam(N, 3, 0), fam(N, 6, 0), fam(N, 6, 6), fam(N, 10, 12)]
    print(N, ' '.join(f'{v:.6f}({f(v):+.4f})' for v in vals), flush=True)
