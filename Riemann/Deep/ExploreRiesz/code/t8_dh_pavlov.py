# T8 adversarial: long-window Pavlov constant [w]_A2 = <|F(s+it)|^2><|F(s+it)|^-2> on Re s = s, t in [60, 420],
# zeta (Carlson/RH: zeta(2s)^2/zeta(4s)) vs Davenport-Heilbronn f (dh_copy.py = explore23's dh.py, Euler-Maclaurin Hurwitz).
import numpy as np, mpmath as mp
from dh_copy import hurwitz, f as fdh
ts = np.arange(60.0, 420.0, 0.005)
print(' sigma | zeta: <w> (zeta(2s)) <1/w> (zeta(2s)/zeta(4s)) A2 (pred) | DH: <w> (sum|a_n|^2 n^-2s) <1/w> A2 | DH max 1/w')
k = (np.sqrt(10 - 2*np.sqrt(5)) - 2)/(np.sqrt(5) - 1); A = {1: 1.0, 2: k, 3: -k, 4: -1.0}
for s in [1.0, 0.9, 0.85, 0.8, 0.75, 0.7, 0.65, 0.6]:
    S = s + 1j*ts
    z = hurwitz(S, 1.0, 260); fz = fdh(S, M=260)
    wz = np.abs(z)**2; wd = np.abs(fz)**2
    pz2 = float(mp.zeta(2*s)); pz4 = float(mp.zeta(2*s)/mp.zeta(4*s))
    pd2 = float(sum(A[r]**2*mp.zeta(2*s, mp.mpf(r)/5) for r in A)*mp.power(5, -2*s))
    print(f'{s:5.2f} | {wz.mean():.4f} ({pz2:.4f}) {(1/wz).mean():.4f} ({pz4:.4f}) {wz.mean()*(1/wz).mean():.4f} ({pz2*pz4:.4f}) | '
          f'{wd.mean():.4f} ({pd2:.4f}) {(1/wd).mean():.4f} {wd.mean()*(1/wd).mean():.4f} | {(1/wd).max():.3e} at t={ts[np.argmax(1/wd)]:.2f}', flush=True)
