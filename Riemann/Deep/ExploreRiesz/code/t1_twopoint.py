# T1 smallest case: Riesz signature of an off-line pair x +- ib vs an on-line pair x +- g/2 (window [-pi, pi], density 1).
import numpy as np
from common import gram, riesz, PI
print('b   lmin(off x+-ib)  lmin(on gap 2b)  ratio   series(2pi b)^2/6')
for b in [0.01, 0.05, 0.1, 0.2, 0.25, 0.3, 0.4, 0.5, 0.75, 1.0, 2.0, 4.0]:
    a = riesz([1j*b, -1j*b])[0]; c = riesz([-b, b])[0]
    print(f'{b:5.2f} {a:.6e} {c:.6e} {a/c:8.5f} {(2*PI*b)**2/6:.6e}')
# pair replacing two lattice points inside a lattice block of K points: x = n0 + 1/2 between n0, n0+1
print('\nblock K=41: lattice with points n0, n0+1 replaced by x+-ib (off) or x+-b (on), x = n0+1/2')
K = 41; n0 = 20
base = np.arange(K, dtype=float)
for b in [0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0, 2.0, 4.0, 8.0]:
    yo = base.astype(complex); yo[n0] = n0 + 0.5 + 1j*b; yo[n0+1] = n0 + 0.5 - 1j*b
    yn = base.copy(); yn[n0] = n0 + 0.5 - b; yn[n0+1] = n0 + 0.5 + b
    lo = riesz(yo); ln_ = riesz(np.sort(yn)) if b <= 0.5 else (np.nan, np.nan)
    print(f'b={b:5.2f} off: lmin={lo[0]:.4e} lmax={lo[1]:.4f} | on(gap 2b): lmin={ln_[0]:.4e} lmax={ln_[1]:.4f} | b*lmin_off={b*lo[0]:.4f}')
