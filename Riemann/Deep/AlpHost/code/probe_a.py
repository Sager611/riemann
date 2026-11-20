# probe (a): quadratic twists xi*chi_d, d prime <= 1e6: max_y |psi(y)|/sqrt(y log y), smoothed truncated S(Z,y), L'/L(1) (Riesz, Z=1e6)
# usage: python3 probe_a.py   (numpy; streams over d; < 50 MB)
import numpy as np, math, sys, time
t0 = time.time()
Z = 10**6
sv = np.ones(Z+1, dtype=bool); sv[:2] = False
for i in range(2, int(Z**0.5)+1):
    if sv[i]: sv[i*i::i] = False
P = np.nonzero(sv)[0].astype(np.int64)            # primes <= Z
# prime powers <= Z: n, log p, base index, exponent
pp_n, pp_lp, pp_b, pp_k = [], [], [], []
for idx, p in enumerate(P):
    k, n = 1, int(p)
    while n <= Z:
        pp_n.append(n); pp_lp.append(math.log(p)); pp_b.append(idx); pp_k.append(k)
        n *= int(p); k += 1
o = np.argsort(pp_n)
N = np.array(pp_n, dtype=np.float64)[o]; LP = np.array(pp_lp)[o]; B = np.array(pp_b)[o]; K = np.array(pp_k)[o]
Ni = N.astype(np.int64)
xi4 = np.where(Ni % 2 == 0, 0.0, np.where(Ni % 4 == 1, 1.0, -1.0))   # chi_{-4}
grid = np.unique(np.round(np.logspace(2, 5, 241))).astype(np.float64)
gs = np.searchsorted(N, grid, side='right')   # count of prime powers <= y
gridS = grid[grid <= 1e4]; gsS = gs[:len(gridS)]
def legendre_primes(d):
    e = (d - 1) // 2
    base = P % d; res = np.ones_like(base); b = base.copy()
    while e:
        if e & 1: res = (res * b) % d
        b = (b * b) % d; e >>= 1
    chi = np.where(res == 1, 1.0, np.where(res == d - 1, -1.0, 0.0))
    return chi
def stats(a):
    cP = np.concatenate(([0.0], np.cumsum(a)))         # psi(y) = cP[count]
    cR = np.concatenate(([0.0], np.cumsum(a / N)))
    psi = cP[gs]
    r1 = np.max(np.abs(psi) / np.sqrt(grid * np.log(grid)))
    RZ, PZ = cR[-1], cP[-1]
    S = cP[gsS] + gridS * (RZ - cR[gsS])             # truncated smoothed sum, Z = 1e6
    r2 = np.max(np.abs(S) / np.sqrt(gridS * np.log(gridS)))
    ld = -(RZ - PZ / Z)                                # Riesz approx of L'/L(1, chi), error ~ Z^-1/2 log d (GRH heuristic)
    # margin test: |psi(y)| <= y^(1/2+eps/2), eps = 1/52
    ok = bool(np.all(np.abs(psi) <= grid ** (0.5 + 1/104)))
    okS = bool(np.all(np.abs(S) <= gridS ** (0.5 + 1/104) / 2))
    return r1, r2, ld, ok, okS, S
mode = sys.argv[1] if len(sys.argv) > 1 else 'sample'
if mode == 'sample':
    D = P[(P > 1000)][::10]
else:   # WinPair-like: d = 7 mod 16, d = -1 mod 105 (z = 7)
    D = P[(P % 16 == 7) & ((P + 1) % 105 == 0)]
print('mode', mode, '#d', len(D), 'd range', int(D[0]), int(D[-1]), 'setup %.1fs' % (time.time()-t0), flush=True)
for M in (1, 4):
    R1, R2, LD, OK, OKS = [], [], [], 0, 0
    msq = np.zeros(len(gridS)); cnt = 0
    for d in D:
        d = int(d)
        chi = legendre_primes(d)
        a = LP * chi[B] ** K
        if M == 4: a = a * xi4
        r1, r2, ld, ok, okS, S = stats(a)
        R1.append(r1); R2.append(r2); LD.append(ld); OK += ok; OKS += okS
        msq += S**2; cnt += 1
    R1, R2, LD = map(np.array, (R1, R2, LD))
    q = lambda v: ' '.join('%.3f' % x for x in np.quantile(v, [0, .1, .5, .9, 1]))
    print('M=%d: max_y|psi(y)|/sqrt(y log y), y in [1e2,1e5], quantiles 0/.1/.5/.9/1: %s' % (M, q(R1)))
    print('M=%d: max_y|S(1e6,y)|/sqrt(y log y), y in [1e2,1e4]: %s' % (M, q(R2)))
    print("M=%d: L'/L(1) (Riesz Z=1e6): %s ; max|.|/log(1e2) = %.3f" % (M, q(LD), np.max(np.abs(LD))/math.log(100)))
    print('M=%d: share with |psi(y)| <= y^(1/2+1/104) on whole grid: %d/%d ; share with |S| <= y^(1/2+1/104)/2 (archive box_of_S test): %d/%d' % (M, OK, cnt, OKS, cnt))
    ms = msq / cnt / (gridS * np.log(gridS))
    print('M=%d: mean over d of |S(1e6,y)|^2/(y log y) at y=%s: %s' % (M, ' '.join('%d' % v for v in gridS[::20]), ' '.join('%.3f' % v for v in ms[::20])))
    print('time %.1fs' % (time.time()-t0), flush=True)
