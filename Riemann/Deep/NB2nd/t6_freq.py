# nb2nd t6: (a) e-coupling: S_N - Sigma_N (Sigma_N = l^T G_N^{-1} l, N-free head problem with K = N) for N in [100,640];
# (b) frequency content of the optimum's C2(N) = L^2 d^2 - C L: least squares on N in [100,640] with
#     pair frequencies (gamma_b - gamma_a) and single-zero terms N^{-1/2} cos/sin(gamma log N).
import math, numpy as np, mpmath as mp
O = '/Users/adriansager/playground/riemann/.scratch/nb2nd/'
Z = np.load(O + 'G640.npz'); G = Z['G']; b = Z['b']; mu = Z['mu'].astype(float); NMAX = len(b)
EG = 0.57721566490153286061; C = 2 + EG - math.log(4 * math.pi); J = 3.5626634565363677
Sig = np.load(O + 't5_sigma.npy'); A2 = np.load(O + 't2_allN.npz'); Ns = A2['N']; V1 = A2['V1']; R = A2['R']
ks = np.arange(1, NMAX + 1)
from scipy.linalg import solve_triangular
LcF = np.linalg.cholesky(G)  # leading blocks = Cholesky factors of G_N
SN = []
for N in range(100, NMAX + 1, 5):
    L = math.log(N); c = -mu[:N] * (1 - np.log(ks[:N]) / L); r = G[:N, :N] @ c - b[:N]
    y = solve_triangular(LcF[:N, :N], r, lower=True); SN.append((N, L * L * float(y @ y), Sig[N - 1]))
SN = np.array(SN); cpl = SN[:, 1] - SN[:, 2]
for lo, hi in [(100, 250), (250, 450), (450, 640)]:
    m = (SN[:, 0] >= lo) & (SN[:, 0] <= hi)
    print(f'(a) N in [{lo},{hi}]: S_N mean {SN[m,1].mean():.4f}; Sigma_N mean {SN[m,2].mean():.4f}; coupling S_N - Sigma_N mean {cpl[m].mean():+.4f} sd {cpl[m].std():.4f}')
gam = [float(mp.zetazero(k).imag) for k in range(1, 7)]
m = (Ns >= 100); x = np.log(Ns[m].astype(float)); yv = V1[m]; yr = R[m]
pairs = [round(gam[1] - gam[0], 3), round(gam[2] - gam[1], 3), round(gam[2] - gam[0], 3)]  # 6.887, 3.989, 10.876 (log N range [4.6, 6.46] cannot resolve more)
def fit(cols, y):
    X = np.vstack(cols).T; coef, *_ = np.linalg.lstsq(X, y, rcond=None); res = y - X @ coef
    return coef, 1 - res.var() / y.var()
base = [np.ones_like(x), 1 / x]
for name, y in (('optimum C2(N)=L^2d^2-CL', yv), ('ramp R(N)=L^2I-CL-J', yr)):
    c0, r0 = fit(base, y)
    cols_p = base + [f(w * x) for w in pairs for f in (np.cos, np.sin)]; cp, rp = fit(cols_p, y)
    cols_z = base + [np.exp(-x / 2) * f(g * x) for g in gam[:3] for f in (np.cos, np.sin)]; cz, rz = fit(cols_z, y)
    cols_b = cols_p + cols_z[2:]; cb, rb = fit(cols_b, y)
    ampp = [round(float(np.hypot(cp[2 + 2 * i], cp[3 + 2 * i])), 5) for i in range(len(pairs))]
    ampz = [round(float(np.hypot(cz[2 + 2 * i], cz[3 + 2 * i])), 4) for i in range(3)]
    print(f'(b) {name}: sd {y.std():.5f}; const+1/L fit: C2 = {c0[0]:+.4f} + {c0[1]:+.4f}/L (R2 {r0:.2f}); '
          f'+pair freqs R2 {rp:.2f} (amps {list(zip(pairs, ampp))}); +N^-1/2 zero terms R2 {rz:.2f} (amps x N^-1/2: {list(zip([round(g,2) for g in gam[:3]], ampz))}); both R2 {rb:.2f}')
