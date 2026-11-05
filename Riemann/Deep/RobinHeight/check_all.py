#!/usr/bin/env python3
# check_all.py -- rechecks every inequality behind REPORT.md (single core, about 7 minutes).
#   uv run --no-project --with python-flint --with mpmath --with numpy python3 check_all.py
# (1) Nicolas' inequality for every primorial with p_k <= 9,999,991, in arb ball arithmetic;
# (2) the analytic certificates (certify.py) on the stated ranges, three half-plane inputs;
# (3) a table of the explicit tail bound and of Lemma 10's eta;
# (4) the Level-2 certificates (certify2.py: ramp-smoothed cutoff + order-k averages) on the stated, larger ranges.
import sys, time, math
import numpy as np
from flint import arb, ctx
import certify as C
t0 = time.time()
ok_all = True

# ---------- (1) direct Nicolas check
ctx.prec = 128
NMAX = 10**7
sv = np.ones(NMAX + 1, dtype=bool); sv[:2] = False
for i in range(2, int(NMAX**0.5) + 1):
    if sv[i]: sv[i*i::i] = False
primes = np.nonzero(sv)[0].tolist()
EG = arb.const_euler(); Ssum = arb(0); TH = arb(0); worst = None
for k, p in enumerate(primes, 1):
    Ssum += -(arb(-1)/p).log1p(); TH += arb(p).log()
    if k == 1: continue                      # N_1 = 2: log log 2 < 0
    m = Ssum - EG - TH.log().log()
    if not (m > 0):
        print('FAIL direct Nicolas at', p); ok_all = False; break
    v = float((m*math.sqrt(p)*math.log(p)).lower())
    if worst is None or v < worst[0]: worst = (v, p)
print(f'(1) PASS: Nicolas holds for all {len(primes)} primorials with p_k <= {primes[-1]} '
      f'(min of margin*sqrt(p)log p = {worst[0]:.4f} at p = {worst[1]}); {time.time()-t0:.0f}s', flush=True)

# ---------- (2) analytic certificates on the stated ranges
ctx.prec = 192
LN10 = arb(10).log()
claims = [  # (label, b0, delta for Lemma 10, X_N, X_R)
    ('7/8',                 arb(7)/8,                       arb(10)**-9,     '1.93e32', '1.56e31'),
    ('0.87338196661320752', arb('0.87338196661320752'),     arb(10)**-9,     '2.67e32', '2.13e31'),
    ('none (b0 = 1)',       arb(1),                         arb('1.25e-8'),  '1.61e24', '2.45e23'),
]
for label, b0, dlt, XN, XR in claims:
    C.setdelta(dlt)
    okN, LN, nbN, mN = C.run('N', b0, 7*LN10, arb(XN).log())
    okR, LR, nbR, mR = C.run('R', b0, 10*LN10, arb(XR).log())
    ok_all &= okN and okR
    eta = C.eta_of(arb(XR).log(), b0)
    log10n = float((arb(XR)*(1 - eta.upper())/LN10).lower())
    print(f'(2) b0 = {label}: Nicolas on [1e7, {XN}]: {"PASS" if okN else "FAIL"} ({nbN} boxes, min normalized margin {mN:.3g}); '
          f'Robin on CA with largest prime in [1e10, {XR}]: {"PASS" if okR else "FAIL"} ({nbR} boxes, min {mR:.3g}); '
          f'=> Robin for 5041 <= n <= 10^(10^{math.log10(log10n):.4f}); {time.time()-t0:.0f}s', flush=True)

# ---------- (3) explicit tail bound table (b0 = 7/8)
C.setdelta(arb(10)**-9)
b0 = arb(7)/8
print('(3) n(x) = sqrt(x) log x;  bounds for n(x) * int_x^inf (t - theta(t)) w(t) dt, and Lemma 10 eta(x):')
for e in (10, 15, 20, 25, 30, 31, 32):
    L = e*LN10
    lo = C.nMlow(L) - C.nZup(L, b0) - (-5*L/2).exp()
    hi = C.nMup(L) + C.nZup(L, b0) + C.LOG2PI*(-L/2).exp()
    tail = C.nZup(L, b0) - (1 + 1/L)*C.BETA0*(1 + (1 + 2*(L+2)/(L*(L+1)))/C.GAM1)
    print(f'    x = 1e{e}: [{float(lo.lower()):.6f}, {float(hi.upper()):.6f}]   zeros above T0 contribute <= {float(tail.upper()):.3e}'
          f'   eta(x) = {float(C.eta_theta(L, b0).upper()):.3e}')
# ---------- (4) Level 2 (smoothed cutoff)
import certify2 as K
claims2 = [  # (label, b0, X_N, X_R)
    ('7/8',                 arb(7)/8,                   '1.04e40', '9.37e38'),
    ('0.87338196661320752', arb('0.87338196661320752'), '1.31e40', '1.19e39'),
    ('none (b0 = 1)',       arb(1),                     '1.42e34', '2.18e33'),
]
for label, b0, XN, XR in claims2:
    e1 = K.eta1_of(b0)
    okN, LN, nbN, mN = K.run2('N', b0, 7*LN10, arb(XN).log(), e1)
    okR, LR, nbR, mR = K.run2('R', b0, 10*LN10, arb(XR).log(), e1)
    ok_all &= okN and okR
    LX = arb(XR).log()
    eta = K.best_eta(LX, b0, e1 if LX > K.LY1 else arb(0))
    log10n = float((arb(XR)*(1 - eta.upper())/LN10).lower())
    print(f'(4) Level 2, b0 = {label}: Nicolas on [1e7, {XN}]: {"PASS" if okN else "FAIL"} ({nbN} boxes, min {mN:.3g}); '
          f'Robin on CA with largest prime in [1e10, {XR}]: {"PASS" if okR else "FAIL"} ({nbR} boxes, min {mR:.3g}); '
          f'=> Robin for 5041 <= n <= 10^(10^{math.log10(log10n):.4f}); eta_1 = {e1.str(4)}; {time.time()-t0:.0f}s', flush=True)
print('ALL PASS' if ok_all else 'SOME CHECK FAILED', f'total {time.time()-t0:.0f}s')
sys.exit(0 if ok_all else 1)
