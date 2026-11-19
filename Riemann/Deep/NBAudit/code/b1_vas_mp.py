# nbaudit b1: independent mpmath recomputation of d_N^2 (Baez-Duarte V1, L^2(0,inf)) for N <= NM from Vasyunin's formula,
# written from the formula alone (direct cot sum over m = 1..k-1, no symmetry trick, no fixed point); 60 digits.
# Also V2 (constraint sum c_k/k = 0) for n <= NM. Compares with the s2-copy certified exact intervals /tmp/nbaudit/loh_V1.txt.
# usage: uv run --no-project --with mpmath --with numpy python3 b1_vas_mp.py 100
import sys, math, time
import mpmath as mp
import numpy as np
t0 = time.time()
NM = int(sys.argv[1]) if len(sys.argv) > 1 else 100
mp.mp.dps = 60
LG = mp.log(2 * mp.pi) - mp.euler
cot = {}
def V(h, k):   # sum_{m=1}^{k-1} {m h / k} cot(pi m / k)
    s = mp.mpf(0)
    for m in range(1, k):
        if (k, m) not in cot: cot[(k, m)] = mp.cot(mp.pi * m / k)
        s += mp.mpf((m * h) % k) / k * cot[(k, m)]
    return s
g = {}
for b in range(1, NM + 1):
    for a in range(1, b + 1):
        if math.gcd(a, b) != 1: continue
        g[(a, b)] = LG if a == b else (LG * (a + b) + (b - a) * mp.log(mp.mpf(a) / b) - mp.pi * (V(a, b) + V(b, a))) / (2 * a * b)
G = [[None] * NM for _ in range(NM)]
for i in range(1, NM + 1):
    for j in range(1, NM + 1):
        d = math.gcd(i, j); G[i - 1][j - 1] = g[(min(i, j) // d, max(i, j) // d)] / d
bv = [(mp.log(k) + 1 - mp.euler) / k for k in range(1, NM + 1)]
def chol_d2(A, rhs, n):   # Cholesky A = L L^T, y = L^-1 rhs; returns cumulative 1 - sum y_i^2 (all leading n at once)
    L = [[mp.mpf(0)] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1):
            s = A[i][j] - mp.fsum(L[i][t] * L[j][t] for t in range(j))
            L[i][j] = mp.sqrt(s) if i == j else s / L[j][j]
    y = []; out = []; acc = mp.mpf(1)
    for i in range(n):
        yi = (rhs[i] - mp.fsum(L[i][t] * y[t] for t in range(i))) / L[i][i]; y.append(yi); acc -= yi * yi; out.append(acc)
    return out, min(L[i][i] ** 2 for i in range(n))
d1, piv1 = chol_d2(G, bv, NM)          # d1[N-1] = d_N^2 (V1)
print(f'[{time.time()-t0:.1f}s] V1 done; min pivot {mp.nstr(piv1, 5)}', flush=True)
# V2: functions rho_k - rho_1/k, k = 2..NM
A = [[G[j][k] - G[j][0] / (k + 1) - G[0][k] / (j + 1) + G[0][0] / ((j + 1) * (k + 1)) for k in range(1, NM)] for j in range(1, NM)]
b2 = [bv[k] - bv[0] / (k + 1) for k in range(1, NM)]
d2, _ = chol_d2(A, b2, NM - 1)         # d2[n-2] = V2 at n
print(f'[{time.time()-t0:.1f}s] V2 done', flush=True)
# compare with the s2-copy certified exact intervals
loh = {}
for line in open('/tmp/nbaudit/loh_V1.txt'):
    N, a, b_, c, d = map(int, line.split()); loh[N] = (mp.mpf(a) / b_, mp.mpf(c) / d)
inside = 0; worst = mp.mpf(0)
for N in range(1, NM + 1):
    lo, hi = loh[N]; v = d1[N - 1]
    ok = lo <= v <= hi; inside += ok
    worst = max(worst, abs(v - (lo + hi) / 2) / ((hi - lo) / 2))
    if not ok: print('OUTSIDE', N, mp.nstr(v, 30), mp.nstr(lo, 30), mp.nstr(hi, 30))
print(f'V1 mpmath value inside certified [lo,hi]: {inside} of {NM}; max |v - mid|/halfwidth = {mp.nstr(worst, 4)}')
print('min (hi - v)/(hi - lo):', mp.nstr(min((loh[N][1] - d1[N-1]) / (loh[N][1] - loh[N][0]) for N in range(1, NM + 1)), 5),
      ' min (v - lo)/(hi - lo):', mp.nstr(min((d1[N-1] - loh[N][0]) / (loh[N][1] - loh[N][0]) for N in range(1, NM + 1)), 5),
      ' max (hi - lo):', mp.nstr(max(loh[N][1] - loh[N][0] for N in range(1, NM + 1)), 5))
Gd = np.load('/tmp/nbaudit/G_double.npy')
print('max |G_mp - G_double(s1 copy)| over N <= NM:', float(max(abs(G[i][j] - Gd[i, j]) for i in range(NM) for j in range(NM))))
C = 2 + mp.euler - mp.log(4 * mp.pi)
sg = [0, 0] + [1 if d1[N - 1] * mp.log(N) > C else -1 for N in range(2, NM + 1)]
ch = sum(1 for N in range(3, NM + 1) if sg[N] != sg[N - 1])
print('V1 sign changes of d^2 logN - C on [2,NM]:', ch, ' first below:', next((N for N in range(2, NM + 1) if sg[N] < 0), None),
      ' min |margin|:', mp.nstr(min(abs(d1[N - 1] * mp.log(N) - C) for N in range(2, NM + 1)), 5))
for N in (10, 20, 40, 41, 48, 49, 80, 100):
    print(f'N={N}: V1 {mp.nstr(d1[N-1], 25)}  V2 {mp.nstr(d2[N-2], 20)}  V1 logN - C {mp.nstr(d1[N-1]*mp.log(N) - C, 8)}')
print(f'[{time.time()-t0:.1f}s] done')
