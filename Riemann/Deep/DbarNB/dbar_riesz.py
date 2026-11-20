# explore59 dbar_riesz.py : Riesz constants of the normalised zero evaluators k_γ/||k_γ|| in V_L (evidence for Lemma R)
import numpy as np, time
from dbar_U import Kblock, zz
t0 = time.time()
K = 800
g = zz[:K]; Zs = np.concatenate([g, -g])
gaps = np.diff(g); print(f'K={K} gamma_K={g[-1]:.1f} min gap={gaps.min():.4f} at gamma={g[np.argmin(gaps)]:.2f}')
for L in (6.0, 7.155, 10.0, 20.0, 40.0, 80.0):
    G = Kblock(Zs, Zs, L)
    d = np.sqrt(np.real(np.diag(G)))
    A = G / d[:, None] / d[None, :]
    ev = np.linalg.eigvalsh((A + A.conj().T) / 2)
    dens = np.log(g[-1] / (2 * np.pi)) / L
    print(f'L={L:6.3f}  zero density/Nyquist at top = {dens:.3f}  Riesz bounds [{ev[0]:.4f}, {ev[-1]:.4f}]')
    del G, A
print('time', round(time.time() - t0, 1), 's')
