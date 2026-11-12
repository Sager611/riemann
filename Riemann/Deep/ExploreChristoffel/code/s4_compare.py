# explore38 s4: compare Burnol/Christoffel lower bound B(N) (s3 scan, K = 400 pairs + diagonal tail) with d_N^2 (s1) for N <= 1280.
import numpy as np, math
OUT = '/Users/adriansager/playground/riemann/.scratch/explore38/'
EG = 0.57721566490153286061; C = 2 + EG - math.log(4 * math.pi)
S = np.load(OUT + 'scan_all400.npy'); d2 = np.load(OUT + 'dN1280.npz')['d2']
N = S[:, 0].astype(int); L = S[:, 1]; B = S[:, 2] + S[:, 4]; BK = S[:, 2]
D = d2[N - 1]
assert np.all(BK <= D), 'bound violated'
print('B_K <= d_N^2 for every N in [2,1280]:', bool(np.all(BK <= D)), ' with tail:', bool(np.all(B <= D)))
r = B / D
c2B = L * L * B - C * L; c2o = L * L * D - C * L; dl = c2o - c2B
for (a, b) in [(2, 10), (10, 100), (100, 250), (250, 450), (450, 640), (640, 1280), (800, 1280)]:
    m = (N >= a) & (N <= b)
    print(f'N in [{a:4d},{b:4d}]: B/d2 mean {r[m].mean():.4f} [{r[m].min():.4f},{r[m].max():.4f}]  L^2B-CL {c2B[m].mean():+.5f} (sd {c2B[m].std():.5f})  L^2d2-CL {c2o[m].mean():+.5f} (sd {c2o[m].std():.5f})  deficit {dl[m].mean():+.5f} (sd {dl[m].std():.5f})')
for n0 in [2, 3, 5, 10, 41, 100, 640, 1280]:
    i = n0 - 2
    print(f'N={n0:5d} d2={D[i]:.6e} B={B[i]:.6e} B/d2={r[i]:.4f} L^2B-CL={c2B[i]:+.5f} L^2d2-CL={c2o[i]:+.5f}')
# oscillation analysis on N in [100,1280]: detrend by a + b/L + c/L^2, then LSQ sinusoids at zero-difference frequencies in L
zc = np.load(OUT + 'zeros_cache.npz'); g = zc['g']
freqs = []
for i in range(6):
    for j in range(i + 1, 7):
        freqs.append((g[j] - g[i], f'g{j+1}-g{i+1}'))
freqs.sort()
m = (N >= 100)
x = L[m]
def fit(y, fl):
    cols = [np.ones_like(x), 1 / x, 1 / x ** 2]
    for f, _ in fl:
        cols += [np.cos(f * x), np.sin(f * x)]
    Mx = np.array(cols).T
    co, res, rk, sv = np.linalg.lstsq(Mx, y, rcond=None)
    yhat = Mx @ co
    R2 = 1 - np.sum((y - yhat) ** 2) / np.sum((y - y.mean()) ** 2)
    amps = [(fl[k][1], fl[k][0], math.hypot(co[3 + 2 * k], co[4 + 2 * k]), math.atan2(-co[4 + 2 * k], co[3 + 2 * k])) for k in range(len(fl))]
    return co, R2, amps, y - Mx[:, :3] @ co[:3]
sel = [f for f in freqs if f[0] < 13]
for name, y in [('L^2B-CL', c2B[m]), ('L^2d2-CL', c2o[m]), ('deficit', dl[m])]:
    co, R2, amps, ydt = fit(y, sel)
    print(f'{name}: trend a + b/L + c/L^2 = {co[0]:+.4f} {co[1]:+.4f}/L {co[2]:+.4f}/L^2; R2 with {len(sel)} freqs = {R2:.3f}')
    print('   ' + '; '.join(f'{nm} {f:.3f}: amp {a:.4f} ph {p:+.2f}' for nm, f, a, p in amps if a > 4e-4))
# correlation of detrended oscillations
_, _, _, yb = fit(c2B[m], [])
_, _, _, yo = fit(c2o[m], [])
print('corr(detrended L^2B, detrended L^2d2) on [100,1280] =', round(float(np.corrcoef(yb, yo)[0, 1]), 3), ' sd:', round(float(yb.std()), 5), round(float(yo.std()), 5))
for a, b in [(100, 400), (400, 1280)]:
    mm = (N[m] >= a) & (N[m] <= b)
    print(f'  corr on [{a},{b}] =', round(float(np.corrcoef(yb[mm], yo[mm])[0, 1]), 3))
np.savez(OUT + 'compare.npz', N=N, L=L, B=B, D=D, c2B=c2B, c2o=c2o)
