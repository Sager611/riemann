# T3 adversarial: Davenport-Heilbronn (304 on-line zeros + 9 off-line pairs to T = 420) vs zeta in the same height range.
# Riesz bounds of K-blocks with the off-line pairs as complex points x +- ib, removed, or replaced by an on-line pair x +- b.
import numpy as np, json
from common import dh_unfolded, zeta_unfolded, riesz, gram, PI, dvartheta_dh, vartheta_dh
pts, d = dh_unfolded(offline=True)
print('DH points', len(pts), 'off-line pairs (b = (beta-1/2)*L):')
for b, g in d['offline']:
    L = dvartheta_dh([g])[0]/PI; print(f'  rho = {b:.4f} + {g:.3f}i  L = {L:.4f}  b = {(b-0.5)*L:.4f}')
real_rep = np.where(np.abs(pts.imag) > 0, pts.real + np.sign(pts.imag)*np.abs(pts.imag), pts.real).real
real_rep = np.sort(real_rep)
dev = pts.real - np.arange(1, len(pts) + 1)
print('DH dev (real parts, all zeros counted): mean %.4f std %.4f range [%.3f, %.3f]' % (dev.mean(), dev.std(), dev.min(), dev.max()))
g, lam = zeta_unfolded()
zl = lam[g < 420]; zd = zl - np.arange(1, len(zl) + 1)
print('zeta t<420: N = %d, dev std %.4f range [%.3f, %.3f]' % (len(zl), zd.std(), zd.min(), zd.max()))
# blocks containing each off-line pair: centered window of K points
print('\nper off-line pair, centered K-blocks: lmin(complex pair) | lmin(pair -> on-line x+-b) | lmin(pair removed)')
offidx = np.nonzero(np.abs(pts.imag) > 0)[0]
for K in [8, 16, 32]:
    rows = []
    for i in offidx[::2]:
        lo = max(0, min(i - K//2, len(pts) - K)); blk = pts[lo:lo+K]
        a = riesz(blk)[0]
        blk2 = np.sort(np.where(np.abs(blk.imag) > 0, blk.real + np.sign(blk.imag)*np.abs(blk.imag), blk.real).real)
        c = riesz(blk2)[0]
        blk3 = blk[np.abs(blk.imag) == 0].real
        e = riesz(blk3)[0]
        rows.append((a, c, e))
    rows = np.array(rows)
    print(f'K={K:3d} ' + ' '.join(f'({r[0]:.3f},{r[1]:.3f},{r[2]:.3f})' for r in rows))
    print(f'      median ratio complex/on-line {np.median(rows[:,0]/rows[:,1]):.4f}')
print('\nall disjoint-ish blocks: median / 10% / min lmin; Kadec fraction')
def st(y, K):
    idx = np.arange(0, len(y) - K + 1, max(1, K//2)); lm = np.array([riesz(y[i:i+K])[0] for i in idx])
    return np.median(lm), np.quantile(lm, 0.1), lm.min(), len(idx)
for K in [4, 8, 16, 32, 64]:
    a = st(pts, K); b = st(pts[np.abs(pts.imag) == 0].real, K); c = st(zl, K); e = st(lam[1000:1000+330], K)
    print(f'K={K:3d} DH(all) {a[0]:.4f} {a[1]:.4f} {a[2]:.4f} | DH on-line only {b[0]:.4f} {b[1]:.4f} {b[2]:.4f} | zeta t<420 {c[0]:.4f} {c[1]:.4f} {c[2]:.4f} | zeta n=1000.. {e[0]:.4f} {e[1]:.4f} {e[2]:.4f}')
# Avdonin block means
def avd(dv, N):
    dv = dv - dv.mean(); nb = len(dv)//N; m = dv[:nb*N].reshape(nb, N).mean(1); return np.abs(m).max(), m.std(), nb
print('\nAvdonin block means (real parts, all zeros): N | DH max, sd | DH on-line only (count deficit) | zeta t<420')
don = pts[np.abs(pts.imag) == 0].real; don = don - np.arange(1, len(don) + 1)
for N in [2, 4, 8, 16, 32, 64]:
    a = avd(dev, N); b = avd(don, N); c = avd(zd, N)
    print(f'{N:3d} DH {a[0]:.3f} {a[1]:.4f} | DH-on {b[0]:.3f} {b[1]:.4f} | zeta {c[0]:.3f} {c[1]:.4f}')
