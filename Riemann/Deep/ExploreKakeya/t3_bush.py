
# t3_bush.py -- Bourgain-Kakeya (fat arithmetic progression) examples are honest Dirichlet polynomials (|b_n| = 1).
# Split [N,2N) into J ~ T^{1/2} intervals I of length L = N T^{-1/2}; D_I(t) = sum_{n in I} n^{it} is ~L on the fat AP
# I* = {t: dist(t, 2 pi n_I Z) <~ T^{1/2}}. Bush: all translations 0, random signs eps_I:
#   D(t) = sum_I eps_I D_I(t), |D|^2 ~ sum_I |D_I|^2 <= J L^2 = N^2 T^{-1/2}  =>  reaches V = N^s only if tau <= 4 - 4s.
# Count large values |W(s)| on 1-separated integer t in [-T/2, T/2] and compare with N^{2-2s} (Montgomery), MVT T N^{1-2s},
# and Guth-Maynard Thm 1.1 (N^{2-2s} + N^{18/5-4s} + T N^{12/5-4s}).
import numpy as np, sys
N = int(sys.argv[1]); tau = float(sys.argv[2])
T = int(round(N**tau)); L = max(1, int(round(N*T**-0.5))); J = N//L
rng = np.random.default_rng(7)
eps = rng.choice([-1.0, 1.0], size=J)
n = np.arange(N, N + J*L, dtype=float); ln = np.log(n)
b = np.repeat(eps, L)
t = np.arange(-T//2, T//2 + 1, dtype=float)
D = np.empty(len(t), dtype=complex)
for a in range(0, len(t), 256):
    D[a:a+256] = np.exp(1j*np.outer(t[a:a+256], ln)) @ b
A = np.abs(D)
kak = N*T**-0.25
print('N=%d tau=%.2f T=%d L=%d J=%d  Kakeya line s*=1-tau/4=%.3f  predicted bush |D|~N T^{-1/4}=%.0f (=N^%.3f)' % (
    N, tau, T, L, J, 1 - tau/4, kak, np.log(kak)/np.log(N)))
print('max|D| = %.0f = N^%.3f at t=%d;  max|D| over |t|>3T^{1/2}: N^%.3f' % (A.max(), np.log(A.max())/np.log(N), t[A.argmax()],
      np.log(A[np.abs(t) > 3*T**0.5].max())/np.log(N)))
print(' s     |W|   N^{2-2s}   MVT T N^{1-2s}   GM Thm1.1   |W| in bush |t|<=2T^{1/2}')
for s in np.arange(0.60, 0.80, 0.02):
    W = A >= N**s
    gm = N**(2-2*s) + N**(18/5-4*s) + T*N**(12/5-4*s)
    print('%.2f  %5d  %8.1f  %10.1f  %10.1f  %5d' % (s, W.sum(), N**(2-2*s), T*N**(1-2*s), gm, (W & (np.abs(t) <= 2*T**0.5)).sum()))
