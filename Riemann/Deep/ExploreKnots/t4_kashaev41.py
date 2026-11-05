# t4_kashaev41.py — Kashaev invariant of the figure-eight knot: <4_1>_N = sum_{k<N} |(q;q)_k|^2, q = e^{2 pi i/N} (a sum of squares,
# so positive); growth exp(N V/(2 pi)) with V = vol(4_1) = (3 sqrt3/2) L(2, chi_{-3}) (an L-value at s = 2, inside absolute convergence).
import numpy as np, mpmath as mp
mp.mp.dps = 30
L2 = mp.dirichlet(2, [0, 1, -1])
V = 3*mp.sqrt(3)/2*L2
print('V = (3 sqrt3/2) L(2,chi_-3) =', mp.nstr(V, 16), '; 2 D(e^{i pi/3}) =', mp.nstr(2*mp.clsin(2, mp.pi/3), 16))
for N in [50, 100, 200, 400, 800, 1600, 3200]:
    j = np.arange(1, N)
    lg = np.concatenate([[0.0], np.cumsum(2*np.log(2*np.abs(np.sin(np.pi*j/N))))])  # log |(q;q)_k|^2, k = 0..N-1
    m = lg.max(); logK = m + np.log(np.sum(np.exp(lg - m)))
    pred = float(N*V/(2*mp.pi)) + 1.5*np.log(N) - 0.25*np.log(3)
    print(f'N={N:5d}: log<4_1>_N = {logK:.6f}; minus [N V/2pi + 1.5 log N - log3/4] = {logK - pred:+.6f}; x N = {(logK-pred)*N:+.4f}')
