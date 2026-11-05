#!/usr/bin/env python3
# sanity.py (not part of the proof; check_all.py is) -- (1) rigorous (arb) check of Nicolas' inequality for every primorial with p_k <= 1e7;
# (2) float sanity check of the explicit-formula bounds (signs) at sample primes;
# (3) float sanity check of the Robin reduction chain at sample colossally abundant numbers.
# Run: uv run --no-project --with python-flint --with numpy --with mpmath python3 sanity.py   (about 1 minute)
import numpy as np, math, time
from flint import arb, ctx
ctx.prec = 128
import certify as C
t0 = time.time()
NMAX = 10**7
sv = np.ones(NMAX + 1, dtype=bool); sv[:2] = False
for i in range(2, int(NMAX**0.5) + 1):
    if sv[i]: sv[i*i::i] = False
primes = np.nonzero(sv)[0].tolist()
print('primes', len(primes), 'largest', primes[-1], f'{time.time()-t0:.1f}s', flush=True)

EG = arb.const_euler()
S = arb(0); TH = arb(0)
worst = None; samples = {}
want = {p for p in primes if p in (2003, 10007, 100003, 1000003, 2000003, 3000017, 5000011, 9999991)}
for k, p in enumerate(primes, 1):
    S += -(arb(-1)/p).log1p()
    TH += arb(p).log()
    if k == 1:
        continue          # N_1 = 2: log log 2 < 0 < log(2/phi(2)); holds trivially
    m = S - EG - TH.log().log()
    if not (m > 0):
        raise SystemExit(f'NICOLAS FAILS OR UNDECIDED at p = {p}: {m}')
    nm = m*math.sqrt(p)*math.log(p)
    if worst is None or nm.lower() < worst[0]:
        worst = (float(nm.lower()), p)
    if p in want:
        samples[p] = (S, TH)
print(f'(1) Nicolas holds (certified) for all primorials with p_k <= {primes[-1]}; min of margin*sqrt(p)log p = {worst}', f'{time.time()-t0:.1f}s', flush=True)

# (2) explicit-formula sanity at sample primes: -I(x) = TrueN - D + S2 must lie in [lower, upper]
print('(2) x, TrueN*n, -I*n, lower*n, upper*n   (n = sqrt x log x)')
b0 = arb(7)/8
for p in sorted(samples):
    if p < 2.1e6: continue
    Sv, THv = samples[p]
    x = arb(p); L = x.log()
    trueN = Sv - EG - THv.log().log()
    D = L.log() - THv.log().log() + (THv - x)/(x*L)
    S2 = arb(1)/(2*x*L)                     # approximation of S_2(x) (sanity only)
    negI = trueN - D + S2
    n = x.sqrt()*L
    lower = C.nMlow(L) - C.nZup(L, b0) + C.LOG2PI*(-L/2).exp()
    upper = C.nMup(L) + C.nZup(L, b0) + C.LOG2PI*(-L/2).exp()
    print(p, float(trueN*n), float(negI*n), float(lower), float(upper), 'OK' if float(lower) <= float(negI*n) <= float(upper) else 'VIOLATION')

# (3) Robin chain at CA numbers N_eps with eps = F1(q) for sample primes q
print('(3) q, x2, true Robin LHS*n, chain(true -I)*n, certified bound*n')
pr = np.array(primes, dtype=np.float64)
for q in (100003, 1000003, 3000017, 9999991):
    eps = math.log1p(1.0/q)/math.log(q)
    # exponents a_p = #{k: F(p,k) >= eps}
    lsig = 0.0; logN = 0.0; x2 = 0.0
    for p in primes:
        if p > q: break
        a = 0; s = 0.0
        while True:
            s += float(p)**(a + 1)
            if math.log1p(1.0/s)/math.log(p) >= eps: a += 1
            else: break
        if a >= 2: x2 = p
        lsig += math.log1p(-float(p)**(-(a + 1))) - math.log1p(-1.0/p)
        logN += a*math.log(p)
    lhs = lsig - float(EG.mid()) - math.log(math.log(logN))
    x = arb(q); L = x.log(); n = float((x.sqrt()*L).mid())
    Sv, THv = samples.get(q, (None, None))
    print(q, x2, lhs*n, 'cert-bound:', end=' ')
    m = C.robin_margin(L, L, b0)
    print(-float(m.mid()), '(certified upper bound for the normalized LHS; must be >= true LHS*n)')

# (4) Level-2 bounds (certify2.py) must also contain the true value of n(x)*(-I(x))
import certify2 as K
P = np.array(primes, dtype=np.float64)
cth = np.cumsum(np.log(P)); clp = np.cumsum(-np.log1p(-1.0/P))
EGf = 0.5772156649015329
b0 = arb(7)/8; e1 = K.eta1_of(b0)
for x in (2100001, 3000017, 5000011, 9999991):
    j = np.searchsorted(P, x, side='right'); th = cth[j-1]; lP = clp[j-1]; L = math.log(x)
    D = math.log(L) - math.log(math.log(th)) + (th - x)/(x*L)
    pg = P[j:]; S2 = float(np.sum(-np.log1p(-1.0/pg) - 1.0/pg)) + 1.0/(2*NMAX*math.log(NMAX))
    negI = (lP - EGf - math.log(math.log(th))) - D + S2
    n = math.sqrt(x)*L
    La = arb(x).log()
    eta = K.best_eta(La, b0, arb(0)); lam = arb(K.lam_opt(La, b0, eta))
    lo = (La/2).exp()*La*K.Mlow3(La + lam, e1) - K.nZphi(La, b0, lam) - K.nWin(La, eta, lam) - (-5*La/2).exp()
    hi = (La/2).exp()*La*K.Mup3(La, e1) + K.nZphi(La, b0, lam) + C.LOG2PI*(-La/2).exp() + K.nWin(La, eta, lam)
    print(x, 'true', round(negI*n, 6), 'L2 bounds', round(float(lo.mid()), 6), round(float(hi.mid()), 6), 'OK' if float(lo.mid()) <= negI*n <= float(hi.mid()) else 'VIOLATION')

print(f'done {time.time()-t0:.1f}s')
