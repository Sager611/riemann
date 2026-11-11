# nb2nd t8: (a) kappa(N) = S_N - Sigma_N for every N in [100,640] and the end gain beyond the fixed head (K=60) vs L;
# (b) A_th = sum_rho |a_rho|^2/|rho|^2, a_rho = zeta''/(2 zeta')(rho): candidate constant of the slow convergence Sigma_K ~ Sigma_inf - A/log K.
import math, numpy as np, mpmath as mp
from scipy.linalg import solve_triangular
O = '/Users/adriansager/playground/riemann/.scratch/nb2nd/'
Z = np.load(O + 'G640.npz'); G = Z['G']; b = Z['b']; mu = Z['mu'].astype(float); NMAX = len(b)
EG = 0.57721566490153286061; C = 2 + EG - math.log(4 * math.pi); J = 3.5626634565363677
Sig = np.load(O + 't5_sigma.npy'); A2 = np.load(O + 't2_allN.npz'); NsA = A2['N']; V1 = A2['V1']; RA = A2['R']
ks = np.arange(1, NMAX + 1); LcF = np.linalg.cholesky(G)
rows = []
for N in range(100, NMAX + 1):
    L = math.log(N); c = -mu[:N] * (1 - np.log(ks[:N]) / L); r = G[:N, :N] @ c - b[:N]
    y = solve_triangular(LcF[:N, :N], r, lower=True); S = L * L * float(y @ y)
    rows.append((N, L, S - Sig[N - 1], J - Sig[N - 1], RA[N - 2], V1[N - 2]))
R = np.array(rows)
for lo, hi in [(100, 250), (250, 450), (450, 640)]:
    m = (R[:, 0] >= lo) & (R[:, 0] <= hi)
    print(f'(a) N in [{lo},{hi}]: kappa mean {R[m,2].mean():.4f} (sd {R[m,2].std():.4f}, s.e. ~{R[m,2].std()/math.sqrt(m.sum()/20):.4f}); J-Sigma_N {R[m,3].mean():.4f}; R(N) {R[m,4].mean():+.4f}; C2 {R[m,5].mean():+.4f}')
X = np.vstack([np.ones(len(R)), 1 / R[:, 1]]).T; coef, *_ = np.linalg.lstsq(X, R[:, 2], rcond=None)
print(f'(a) kappa(N) ~ k_inf - a/L fit on [100,640]: k_inf = {coef[0]:.4f}, a = {-coef[1]:.4f}')
mp.mp.dps = 15
tot = 0.0; totre = 0.0
for k in range(1, 101):
    rho = mp.zetazero(k); z1 = mp.zeta(rho, derivative=1); z2 = mp.zeta(rho, derivative=2); a = z2 / (2 * z1)
    w = 2 / abs(rho) ** 2; tot += w * abs(a) ** 2; totre += w * float(mp.re(a)) ** 2
    if k in (1, 2, 3, 10, 30, 100):
        print(f'(b) k={k}: gamma={float(rho.imag):.3f} a={complex(a):.4f} Re a + log(gamma/2pi)/2 = {float(mp.re(a)) + math.log(float(rho.imag)/(2*math.pi))/2:+.2e}; partial A_th = {tot:.5f} (Re-part only {totre:.5f})')
