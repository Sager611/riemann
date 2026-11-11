# Tests of the LiRate predictions on lambda_n, n <= NMAX (FFT data, logs/lam_r999.txt).
# usage: python3 analyze.py logs/lam_r999.txt > logs/analyze.txt
import sys, math
import numpy as np
import mpmath as mp

f = sys.argv[1]
lam = {}
for l in open(f):
    if l.startswith('#'): continue
    n, v = l.split()[:2]; lam[int(n)] = float(v)
NMAX = max(lam)
n = np.arange(1, NMAX + 1, dtype=float)
L = np.array([lam[k] for k in range(1, NMAX + 1)])
g = 0.57721566490153286
kap = 1 + math.log(2 * math.pi) - g
R = L - n / 2 * (np.log(n) - kap)
print(f'# data {f}, n <= {NMAX}')
print('R(n) at n=100,500,1000,10^4,5*10^4,10^5:', [round(R[k - 1], 6) for k in (100, 500, 1000, 10000, 50000, 100000) if k <= NMAX])

# --- c0: lambda_n^(0) - t(n) -> c0 = -1/4 (proofs.md, Lemma A); predicted mean R = c0 + 2 = 7/4 -------------
mp.mp.dps = 30
def Ap(t):  # A'(t) = theta_RS'(t)/pi
    return (mp.re(mp.digamma(mp.mpf(1) / 4 + 1j * t / 2)) / 2 - mp.log(mp.pi) / 2) / mp.pi
def a(v):
    t = mp.cot(v / 2) / 2
    return Ap(t) / (4 * mp.sin(v / 2) ** 2)
for nn in (50, 200):
    pts = [mp.pi * k / (4 * nn) for k in range(0, 4 * nn + 1)]
    lam0 = mp.quad(lambda v: 2 * (1 - mp.cos(nn * v)) * a(v), pts)
    print(f'  n={nn}: lambda0 - t(n) = {mp.nstr(lam0 - mp.mpf(nn) / 2 * (mp.log(nn) - kap), 12)}  (-1/4 - 1/(24n) = {-0.25 - 1/(24*nn):.12f})')

# --- mean and mean square ---------------------------------------------------------------------
def T1L1(N):  # T^2 log(T/2pi) = N
    T = math.sqrt(N)
    for _ in range(60): T = math.sqrt(N / math.log(T / (2 * math.pi)))
    return T, math.log(T / (2 * math.pi))
print('\n N      mean(R)   mean(R^2)  M/sqrt(NlogN)  Mc/sqrt(NlogN)  pred(4/3pi)sqrt(N L1)  ratio   ThmI ratio  max|R|/(nlogn)^.25')
cs = np.cumsum(R); cs2 = np.cumsum(R ** 2)
rows = []
for N in (500, 1000, 2000, 5000, 10000, 20000, 50000, 100000):
    if N > NMAX: break
    m1 = cs[N - 1] / N; m2 = cs2[N - 1] / N; mc = m2 - m1 ** 2
    T1, L1 = T1L1(N); pr = 4 / (3 * math.pi) * math.sqrt(N * L1)
    thm = cs2[N - 1] / (0.0558 * N ** 1.5 * math.log(N) ** 0.5)
    mx = np.max(np.abs(R[9:N]) / (n[9:N] * np.log(n[9:N])) ** 0.25)
    rows.append((N, m2, mc))
    print(f'{N:6d} {m1:9.4f} {m2:10.3f} {m2/math.sqrt(N*math.log(N)):10.4f} {mc/math.sqrt(N*math.log(N)):12.4f} {pr:16.3f} {mc/pr:12.3f} {thm:9.2f} {mx:12.3f}')
# dyadic blocks: centred mean square on (N/2, N]
print('\n dyadic block (N/2,N]:  mean   centred mean square  /sqrt(N log N)   exponent fit')
prev = None
for N in (1000, 2000, 4000, 8000, 16000, 32000, 64000, 100000):
    if N > NMAX: break
    seg = R[N // 2:N]; mu = seg.mean(); v = ((seg - mu) ** 2).mean()
    ex = '' if prev is None else f'{math.log(v / prev[1]) / math.log(N / prev[0]):.3f}'
    print(f'{N:6d} {mu:8.3f} {v:12.3f} {v/math.sqrt(N*math.log(N)):10.4f}   {ex}')
    prev = (N, v)

# pointwise bound check  |R| <= 0.170 sqrt(n) log n / sqrt(loglog n)  (asymptotic, RH)
k = np.arange(100, NMAX + 1)
pb = np.max(np.abs(R[k - 1]) / (0.170 * np.sqrt(k) * np.log(k) / np.sqrt(np.log(np.log(k)))))
print(f'\nmax_{{100<=n<=NMAX}} |R(n)| / (0.170 sqrt n log n / sqrt(loglog n)) = {pb:.4f}')
i = int(np.argmax(np.abs(R[99:]) / (k * np.log(k)) ** 0.25)) + 100
print(f'argmax |R(n)|/(n log n)^(1/4) over n>=100: n={i}, R={R[i-1]:.4f}, ratio={abs(R[i-1])/(i*math.log(i))**0.25:.4f}')

# --- spectral lines in u = sqrt n -------------------------------------------------------------
mu = R[999:].mean()
n0 = 100
nn = n[n0 - 1:]; u = np.sqrt(nn); F = (R[n0 - 1:] - mu) / nn ** 0.25
du = 1 / (2 * u); U0, U1 = u[0], u[-1]
w = np.sin(np.pi * (u - U0) / (U1 - U0)) ** 2
W = np.sum(w * du)
def C(om):
    return np.sum(F * w * du * np.exp(-1j * om * u)) / W
def Lam(m):
    for p in range(2, m + 1):
        if m % p == 0:
            q = m
            while q % p == 0: q //= p
            return math.log(p) if q == 1 else 0.0
    return 0.0
print('\n spectral lines of F(u) = (R(u^2) - mean)/sqrt(u), Hann window on u in [%.1f, %.1f]' % (U0, U1))
print(' m   omega_m   |C| data   |C| pred   arg data   arg pred(pi/4=0.785)')
for m in range(2, 32):
    if Lam(m) == 0: continue
    om = 2 * math.sqrt(math.log(m))
    bm = Lam(m) / math.sqrt(m) * math.log(m) ** -0.75
    c = C(om)
    print(f'{m:3d} {om:8.4f} {abs(c):10.4f} {bm/(2*math.sqrt(math.pi)):10.4f} {np.angle(c):9.3f}')
# periodogram peaks on a grid
grid = np.arange(0.3, 4.0, 0.002)
P = np.array([abs(C(o)) for o in grid])
pk = [i for i in range(1, len(grid) - 1) if P[i] > P[i - 1] and P[i] > P[i + 1] and P[i] > 0.03]
print('\n periodogram local maxima with |C| > 0.03 (omega, |C|, nearest omega_m):')
oms = [(2 * math.sqrt(math.log(m)), m) for m in range(2, 200) if Lam(m) > 0]
for i in pk:
    om, m = min(oms, key=lambda x: abs(x[0] - grid[i]))
    print(f'  {grid[i]:.3f} {P[i]:.4f}   m={m} ({om:.3f})')
print(f' median |C| on grid: {np.median(P):.4f}')
