# explore38 s8: DH vs zeta: deficit delta(N) = L^2 (d_N^2 - B(N)) vs L, detrended oscillation correlation.
import numpy as np, math, sys
OUT = '/Users/adriansager/playground/riemann/.scratch/explore38/'
EG = 0.57721566490153286061; C = 2 + EG - math.log(4 * math.pi)
Cf = 0.19527  # explore13 DH constant sum 1/(rho(1-rho))
fn = sys.argv[1] if len(sys.argv) > 1 else 'dh_compare.npy'
R = np.load(OUT + fn); N = R[:, 0]; L = R[:, 1]; d2 = R[:, 2]; Bon = R[:, 3]; B = R[:, 4]
dl = L * L * (d2 - B)
print('DH deficit by window (L^2 (d2 - B_all)):')
for a, b in [(5, 20), (20, 50), (50, 100), (100, 200), (200, 400), (400, 700), (700, 1000)]:
    m = (N >= a) & (N <= b)
    if m.sum() > 2:
        sl = np.polyfit(L[m], dl[m], 1)[0]
        print(f'  N in [{a},{b}] L in [{L[m].min():.2f},{L[m].max():.2f}]: deficit mean {dl[m].mean():.4f}, local slope d(deficit)/dL {sl:.4f}, deficit/L {np.mean(dl[m]/L[m]):.4f}; B/d2 {np.mean(B[m]/d2[m]):.4f}; d2 L {np.mean(d2[m]*L[m]):.4f}; B L {np.mean(B[m]*L[m]):.4f}')
m = N >= 50
p1 = np.polyfit(L[m], dl[m], 1); p2 = np.polyfit(1 / L[m], dl[m], 1)
r1 = dl[m] - np.polyval(p1, L[m]); r2 = dl[m] - np.polyval(p2, 1 / L[m])
print(f'fit deficit = {p1[1]:.4f} + {p1[0]:.4f} L: resid sd {r1.std():.5f};  fit = {p2[1]:.4f} + {p2[0]:.4f}/L: resid sd {r2.std():.5f}')
def detr(y, x):
    X = np.vstack([np.ones_like(x), 1 / x, 1 / x ** 2, x]).T
    co = np.linalg.lstsq(X, y, rcond=None)[0]; return y - X @ co
yb = detr(L[m] ** 2 * B[m] - Cf * L[m], L[m]); yo = detr(L[m] ** 2 * d2[m] - Cf * L[m], L[m])
print(f'DH detrended corr(L^2 B, L^2 d2) on N>=50: {np.corrcoef(yb, yo)[0,1]:.3f}; sd {yb.std():.5f} {yo.std():.5f}')
Z = np.load(OUT + 'compare.npz'); Nz = Z['N']; Lz = Z['L']; dz = Lz * Lz * (Z['D'] - Z['B'])
mz = Nz >= 50
p1z = np.polyfit(Lz[mz], dz[mz], 1)
print(f'zeta: fit deficit = {p1z[1]:.4f} + {p1z[0]:.4f} L on N in [50,1280] (resid sd {np.std(dz[mz]-np.polyval(p1z, Lz[mz])):.5f})')
for a, b in [(5, 20), (20, 50), (50, 100)]:
    mm = (Nz >= a) & (Nz <= b)
    print(f'  zeta N in [{a},{b}]: deficit mean {dz[mm].mean():.4f}')
