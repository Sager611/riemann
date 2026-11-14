
# Independent float64 recomputation of C_zeta = sum_p M_p (Prop. 4.5 of Casana-Giner) for the six kernels:
# p <= 1000: phase maximum on a 2^16 grid (no derivative correction, so this is ~ the true value, slightly below the
# certified bound); 1000 < p: Euclidean majorant H_{p,k}; n > N0 = 1e8: the paper's tail formula (4.10) with G = the
# non-increasing hull of g.
import json, math, sys, time
import numpy as np
K = json.load(open('Riemann/Deep/STCAudit/kernels.json'))
N0 = 10**8
t0 = time.time()
sv = np.ones(N0//2, dtype=bool); sv[0] = False   # sv[i] <-> 2i+1
for i in range(1, int(math.isqrt(N0))//2 + 1):
    if sv[i]:
        p = 2*i+1; sv[p*p//2::p] = False
pr = np.concatenate([[2], 2*np.nonzero(sv)[0]+1]).astype(np.int32)
del sv
print('primes', len(pr), 'time %.1f' % (time.time()-t0)); sys.stdout.flush()
lp = np.log(pr.astype(np.float64))
# psi(N0)
psiN0 = 0.0
psiN0 = float(lp.sum())
for k in range(2, 30):
    q = pr[:np.searchsorted(pr, int(N0**(1/k))+2)].astype(np.float64)
    q = q[q**k <= N0]
    if len(q) == 0: break
    psiN0 += float(np.log(q).sum())
print('psi(1e8) = %.6f  (theta-check: psi/N0 = %.8f)' % (psiN0, psiN0/N0))
def polys(k):
    c = dict(((X, j), cc) for (X, j), cc in zip(k['B'], k['c']))
    pa = lambda L: -sum(cc*L**(j-1)/math.factorial(j-1) for (X, j), cc in c.items() if X == 'R')/math.pi
    pb = lambda L: ((np.exp(-k['d']*L) - 2)/L - sum(cc*L**(j-1)/math.factorial(j-1) for (X, j), cc in c.items() if X == 'I'))/math.pi
    return pa, pb
table2 = {'L5': 0.70360002, 'L20': 1.46807534, 'L60': 2.10426949, 'L150': 2.25578346, 'C400': 4.46374163, 'D6000': 4.58951035}
res = {}
for name, k in K.items():
    d = k['d']; s1 = 0.5 + d
    a, b = polys(k)
    h = lambda L: np.sqrt(a(L)**2 + b(L)**2)
    # small primes: phase max
    x = np.linspace(0, 2*np.pi, 2**16, endpoint=False)
    small = 0.0
    for p in pr[pr <= 1000]:
        f = np.zeros_like(x); kk = 1
        while p**kk <= N0:
            L = kk*math.log(p); lam = math.log(p)*p**(-kk*s1)
            f += lam*(a(L)*np.cos(kk*x) + b(L)*np.sin(kk*x)); kk += 1
        small += f.max()
    # medium primes: Euclidean, all k with p^k <= N0
    med = 0.0
    big = np.nonzero(pr > 1000)[0]
    for ch in np.array_split(big, 40):
        L = lp[ch]; med += float(np.sum(L*h(L)*np.exp(-s1*L)))
    sm = pr[(pr > 1000) & (pr <= 10**4)]
    for kk in range(2, 30):
        q = sm[sm.astype(np.float64)**kk <= N0]
        if len(q) == 0: break
        L = kk*np.log(q.astype(np.float64)); med += float(np.sum(np.log(q.astype(np.float64))*h(L)*np.exp(-s1*L)))
    # tail: g(x) = h(log x) x^{-s1}; G = non-increasing hull; Tail = G(N0)(1.03883 N0 - psi(N0)) + 1.03883 int_{N0}^inf G
    L = np.linspace(math.log(N0), 400.0, 400_001)
    lg = np.log(h(L)) - s1*L                     # log g
    hull = np.maximum.accumulate(lg[::-1])[::-1]  # log G (non-increasing hull)
    integrand = np.exp(hull + L)                  # G(e^L) e^L
    I = np.trapezoid(integrand, L)
    tail = math.exp(hull[0])*(1.03883*N0 - psiN0) + 1.03883*I
    mono = bool(np.all(np.diff(lg) <= 1e-15))
    Cz = small + med + tail
    res[name] = dict(small=small, med=med, tail=tail, tail_int=I, g_monotone_from_N0=mono, Cz_est=Cz, Cz_table2=table2[name])
    print('%-5s small(p<=1000, phase max)=%.8f  med=%.8f  tail=%.8f (int %.8f, g nonincr from N0: %s)  Cz_est=%.8f  Table2 Cz<%.8f  diff=%.2e' % (name, small, med, tail, I, mono, Cz, table2[name], table2[name]-Cz))
    sys.stdout.flush()
json.dump(res, open('Riemann/Deep/STCAudit/cz_check.json', 'w'), indent=1)
