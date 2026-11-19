# nbaudit b2: Gram entries by DIRECT exact integration, independent of Vasyunin's cotangent formula:
# G_hk = int_0^inf {u/h}{u/k} du/u^2, exact piecewise on [0, J L] (L = lcm(h,k)), periodic tail sum_{j>=J} int_0^L f(v)/(jL+v)^2 dv
# = sum_n (n+1)(-1)^n m_n zeta(n+2, J)/L^(n+2), m_n = int_0^L f v^n dv exact rationals. All pairs h <= k <= 20 (non-coprime too, so homogeneity
# is tested) -> d_N^2 for N <= 20 vs the s2-copy certified intervals; plus random coprime pairs up to 640 vs Vasyunin (mpmath) and s1's G_double.
# usage: uv run --no-project --with mpmath --with numpy python3 b2_direct.py
import sys, math, time, random
import mpmath as mp
import numpy as np
from fractions import Fraction as Fr
t0 = time.time(); mp.mp.dps = 40
def gram_direct(h, k, J=6, NT=62):
    L = h * k // math.gcd(h, k)
    bps = sorted(set(range(0, L + 1, h)) | set(range(0, L + 1, k)))
    pcs = [(bps[i], bps[i + 1], bps[i] // h, bps[i] // k) for i in range(len(bps) - 1)]
    Lh, Lk = L // h, L // k; tot = mp.mpf(0)
    for j in range(J):
        for (a, b, p, q) in pcs:
            A, B, P, Q = a + j * L, b + j * L, p + j * Lh, q + j * Lk
            if A == 0: tot += mp.mpf(B) / (h * k); continue          # P = Q = 0: integrand 1/(hk)
            tot += mp.mpf(B - A) / (h * k) - mp.mpf(Q * k + P * h) / (h * k) * mp.log(mp.mpf(B) / A) + P * Q * (mp.mpf(1) / A - mp.mpf(1) / B)
    tail = mp.mpf(0); hk = h * k
    for n in range(NT):
        S3 = sum(b ** (n + 3) - a ** (n + 3) for (a, b, p, q) in pcs)
        S2 = sum((q * k + p * h) * (b ** (n + 2) - a ** (n + 2)) for (a, b, p, q) in pcs)
        S1 = sum(p * q * (b ** (n + 1) - a ** (n + 1)) for (a, b, p, q) in pcs)
        mn = Fr(S3 * (n + 2) * (n + 1) - S2 * (n + 3) * (n + 1) + S1 * (n + 3) * (n + 2) * hk, (n + 1) * (n + 2) * (n + 3) * hk)
        tail += (n + 1) * (-1) ** n * (mp.mpf(mn.numerator) / mn.denominator) * mp.zeta(n + 2, J) / mp.mpf(L) ** (n + 2)
    return tot + tail
LG = mp.log(2 * mp.pi) - mp.euler
def vas(h, k):   # Vasyunin, coprime h < k
    V = lambda a, c: mp.fsum(mp.mpf((m * a) % c) / c * mp.cot(mp.pi * m / c) for m in range(1, c))
    return LG if h == k else (LG * (h + k) + (k - h) * mp.log(mp.mpf(h) / k) - mp.pi * (V(h, k) + V(k, h))) / (2 * h * k)
NM = 20
G = [[None] * NM for _ in range(NM)]
dmax = mp.mpf(0)
for i in range(1, NM + 1):
    for j in range(i, NM + 1):
        G[i - 1][j - 1] = G[j - 1][i - 1] = gram_direct(i, j)
        d = math.gcd(i, j); dmax = max(dmax, abs(G[i - 1][j - 1] - vas(i // d, j // d) / d))
print(f'[{time.time()-t0:.1f}s] direct Gram N<=20 done; max |direct - Vasyunin/d| over all 210 pairs = {mp.nstr(dmax, 5)}', flush=True)
bv = [(mp.log(k) + 1 - mp.euler) / k for k in range(1, NM + 1)]
L = [[mp.mpf(0)] * NM for _ in range(NM)]; y = []; acc = mp.mpf(1); d1 = []
for i in range(NM):
    for j in range(i + 1):
        s = G[i][j] - mp.fsum(L[i][t] * L[j][t] for t in range(j))
        L[i][j] = mp.sqrt(s) if i == j else s / L[j][j]
    yi = (bv[i] - mp.fsum(L[i][t] * y[t] for t in range(i))) / L[i][i]; y.append(yi); acc -= yi * yi; d1.append(acc)
loh = {}
for line in open('/tmp/nbaudit/loh_V1.txt'):
    N, a, b_, c, d = map(int, line.split())
    if N <= NM: loh[N] = (mp.mpf(a) / b_, mp.mpf(c) / d)
print('direct-Gram d_N^2 inside certified [lo,hi]:', sum(loh[N][0] <= d1[N - 1] <= loh[N][1] for N in range(1, NM + 1)), 'of', NM,
      '; d_20^2 =', mp.nstr(d1[19], 25), '; d_10^2 =', mp.nstr(d1[9], 25))
random.seed(5); Gd = np.load('/tmp/nbaudit/G_double.npy'); worst_v = mp.mpf(0); worst_d = 0.0; npair = 0
while npair < 24:
    h, k = sorted(random.sample(range(2, 641), 2))
    if math.gcd(h, k) != 1 or h * k > 120000: continue
    gd = gram_direct(h, k, J=4, NT=80); gv = vas(h, k)
    worst_v = max(worst_v, abs(gd - gv)); worst_d = max(worst_d, abs(float(gd) - Gd[h - 1, k - 1])); npair += 1
print(f'[{time.time()-t0:.1f}s] 24 random coprime pairs up to 640 (hk <= 1.2e5): max |direct - Vasyunin| = {mp.nstr(worst_v, 5)}; max |direct - G_double(s1)| = {worst_d:.3e}')
