# explore44 c7: rotation invariance. D_{t0}(t) = sum_{n<=N} n^{-1/2 - i(t - t0)}: every t0 gives the SAME law on T^inf (z_p -> p^{i t0} z_p is
# a Haar-preserving rotation), so every T^inf functional inequality treats them alike; on [T,2T] their 4th moments differ once N > T^{1/2}.
import numpy as np, time
t0c = time.time()
T = 2000.0; t = np.arange(T, 2*T, 0.01)
print('T = %g; columns: N, theta, Haar E|D|^4, [T,2T] E|D|^4 for t0 = 0 (arithmetic, ~ zeta partial sum), t0 = 1.5T (rotated)' % T)
for N in [20, 45, 100, 200, 400, 1000]:
    n = np.arange(1, N + 1, dtype=float); ln = np.log(n); w = n**-0.5
    r = {}
    for a_ in range(1, N + 1):
        for b_ in range(1, N + 1):
            r[a_*b_] = r.get(a_*b_, 0) + 1
    h4 = float(sum(v*v/m for m, v in r.items()))
    res = []
    for tz in [0.0, 1.5*T]:
        m4 = 0.0
        for i in range(0, len(t), 2000):
            D = (w[None, :]*np.exp(-1j*np.outer(t[i:i+2000] - tz, ln))).sum(1)
            m4 += (np.abs(D)**4).sum()
        res.append(m4/len(t))
    print('  N=%5d theta=%.2f Haar %9.2f | t0=0: %9.2f (ratio %.2f) | t0=1.5T: %9.2f (ratio %.2f)' % (N, np.log(N)/np.log(T), h4, res[0], res[0]/h4, res[1], res[1]/h4), flush=True)
print('done', round(time.time() - t0c, 1))
