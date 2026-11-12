# T2 at scale: finite-section Riesz bounds (normalized sinc Gram on L^2(-pi,pi)) of K consecutive unfolded points;
# Kadec (range of deviations < 1/2 after centering) and Avdonin (block means) statistics. zeta vs CUE vs Poisson vs DH.
import numpy as np, sys
from common import zeta_unfolded, dh_unfolded, cue, riesz, PI
rng = np.random.default_rng(41)
g, lam = zeta_unfolded()
dev = lam - np.arange(1, len(lam) + 1)          # = -S(gamma_n)
print('zeta: N =', len(g), 'dev mean %.4f std %.4f min %.3f max %.3f' % (dev.mean(), dev.std(), dev.min(), dev.max()))
for lo, hi in [(0, 1000), (1000, 5000), (5000, 22491)]:
    d = dev[lo:hi]; print(f'  n in [{lo},{hi}) t in [{g[lo]:.0f},{g[hi-1]:.0f}]: std(S) = {d.std():.4f}')
# reference sets of the same size (density 1)
cu = np.concatenate([cue(1024, rng) + 1024*k for k in range(12)])     # 12 independent CUE(1024), unfolded
po = np.cumsum(rng.exponential(1.0, 12288))
def stats(y, K, maxblocks=400, start=0):
    lm, lM, kad = [], [], []
    nb = min(maxblocks, (len(y) - start)//K)
    idx = np.linspace(start, len(y) - K, nb).astype(int)
    for i in idx:
        b = y[i:i+K]
        a, B = riesz(b); lm.append(a); lM.append(B)
        dv = np.real(b) - np.arange(K)
        kad.append(np.ptp(dv - np.round(np.mean(dv))) < 0.5)
    lm = np.array(lm); lM = np.array(lM)
    return np.median(lm), np.quantile(lm, 0.1), lm.min(), np.median(lM), np.mean(kad), nb
print('\nK | set | median lmin | 10% lmin | min lmin | median lmax | Kadec-applicable fraction | blocks')
for K in [2, 4, 8, 16, 32, 64, 128]:
    for name, y, st in [('zeta n>1000', lam, 1000), ('CUE', cu, 0), ('Poisson', po, 0)]:
        r = stats(y, K, start=st)
        print(f'{K:4d} {name:12s} {r[0]:.4e} {r[1]:.4e} {r[2]:.4e} {r[3]:.4f} {r[4]:.3f} {r[5]}')
    sys.stdout.flush()
# Avdonin: block means of deviations over N consecutive points (after removing the global mean)
print('\nAvdonin block means: N | zeta: max|mean-c|, frac<1/4 | CUE (one CUE(12288)) | Poisson')
def avd(devs, N):
    # devs: list of deviation arrays (each relative to its own lattice); block means over disjoint blocks of N
    ms = []
    for d in devs:
        d = d - d.mean()
        nb = len(d)//N
        if nb == 0: continue
        ms.append(d[:nb*N].reshape(nb, N).mean(1))
    ms = np.concatenate(ms)
    return np.abs(ms).max(), np.mean(np.abs(ms) < 0.25), ms.std(), len(ms)
dz = [dev[1000:]]
dc = [cu[k*1024:(k+1)*1024] - (np.arange(1024) + 1024*k) for k in range(12)]
dp = [po - np.arange(1, len(po) + 1)]
for N in [1, 2, 4, 8, 16, 32, 64, 128, 256, 512]:
    a = avd(dz, N); b = avd(dc, N); c = avd(dp, N)
    print(f'{N:4d} zeta {a[0]:.3f} {a[1]:.3f} sd {a[2]:.4f} | CUE {b[0]:.3f} {b[1]:.3f} sd {b[2]:.4f} | Poisson {c[0]:.2f} {c[1]:.3f} sd {c[2]:.3f}')
