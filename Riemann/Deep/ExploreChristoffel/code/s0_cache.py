# explore38 s0: cache zeta zeros (n <= 300 exact from mpmath, rest from ExploreHelio file, |err| < 1e-9 for n > 100) and theta, theta' for n <= 3000.
import numpy as np, mpmath as mp, time
mp.mp.dps = 30
t0 = time.time()
zf = np.load('/Users/adriansager/playground/riemann/Riemann/Deep/ExploreHelio/results/zeta_zeros.npy')
g = zf.copy()
for n in range(1, 301):
    g[n - 1] = float(mp.zetazero(n).imag)
K = 3000
th = np.array([float(mp.siegeltheta(t)) for t in g[:K]])
th1 = np.array([float(mp.siegeltheta(t, derivative=1)) for t in g[:K]])
np.savez('/Users/adriansager/playground/riemann/.scratch/explore38/zeros_cache.npz', g=g, th=th, th1=th1)
print('cached', len(g), 'zeros; theta for', K, f'{time.time()-t0:.1f}s')
