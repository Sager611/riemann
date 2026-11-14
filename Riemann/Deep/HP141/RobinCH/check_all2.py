#!/usr/bin/env python3
# check_all2.py -- rechecks every inequality behind RobinCH/REPORT.md on one core.
#   uv run --no-project --with python-flint --with mpmath --with numpy python3 check_all2.py
# (0) the Chirre-Helfgott constants (C_T <= 113.67, the pi/(T-1) step) and eta_1;
# (1) Nicolas' inequality for every primorial with p_k <= 9,999,991, in arb ball arithmetic (as RobinHeight/check_all.py);
# (2) Level 2-CH: RobinHeight's Level 2 with eta from Lemma 12 (Chirre-Helfgott), three zero inputs;
# (3) Level 3: centred k-box ramp (Lemma 11_k) + Lemma 13, three zero inputs.
import sys, time, math
import numpy as np
from flint import arb, ctx
t0 = time.time()
ok_all = True

# ---------- (0) Chirre-Helfgott constants
ctx.prec = 192
PI = arb.pi()
T = 3*arb(10)**12 + 1 + PI/3
CT = ((T/(2*PI)).log()**2)/(2*PI) - (T/(2*PI)).log()/(6*PI)
c1 = CT < arb('113.67')
c2 = PI**2/(3*T**2) + PI/(T - 1) <= PI/(3*arb(10)**12)
c3 = T <= arb(3000175332800)
ok0 = bool(c1 and c2 and c3)
ok_all &= ok0
print(f'(0) {"PASS" if ok0 else "FAIL"}: C_T = {CT.str(10)} < 113.67: {c1}; pi^2/(3T^2) + pi/(T-1) <= pi/3e12: {c2}; '
      f'T = 3e12 + 1 + pi/3 <= 3000175332800 (Platt-Trudgian): {c3}', flush=True)

# ---------- (1) direct Nicolas check
ctx.prec = 128
NMAX = 10**7
sv = np.ones(NMAX + 1, dtype=bool); sv[:2] = False
for i in range(2, int(NMAX**0.5) + 1):
    if sv[i]: sv[i*i::i] = False
primes = np.nonzero(sv)[0].tolist()
EG = arb.const_euler(); Ssum = arb(0); TH = arb(0); worst = None; ok1 = True
for k, p in enumerate(primes, 1):
    Ssum += -(arb(-1)/p).log1p(); TH += arb(p).log()
    if k == 1: continue                      # N_1 = 2: log log 2 < 0
    m = Ssum - EG - TH.log().log()
    if not (m > 0):
        print('FAIL direct Nicolas at', p); ok1 = False; break
    v = float((m*math.sqrt(p)*math.log(p)).lower())
    if worst is None or v < worst[0]: worst = (v, p)
ok_all &= ok1
print(f'(1) {"PASS" if ok1 else "FAIL"}: Nicolas holds for all {len(primes)} primorials with p_k <= {primes[-1]} '
      f'(min of margin*sqrt(p)log p = {worst[0]:.4f} at p = {worst[1]}); {time.time()-t0:.0f}s', flush=True)

# ---------- (2), (3) analytic certificates
ctx.prec = 192
import certify3 as Q
LN10 = Q.LN10
e1 = Q.eta1_ch()
print(f'    eta_1 = eta_CH(1e19) = {e1.str(6)}; eta_CH(1e40) = {Q.eta_ch(40*LN10, e1).str(6)}', flush=True)
claims = [  # (part, level, label, b0, X_N, X_R)
    ('(2) Level 2-CH', 2, '0.87338196661320752', '7.17e41', '8.47e40'),
    ('(2) Level 2-CH', 2, '7/8',                 '5.99e41', '7.11e40'),
    ('(2) Level 2-CH', 2, '1',                   '3.57e36', '5.51e35'),
    ('(3) Level 3',    4, '0.87338196661320752', '1.78e48', '1.43e47'),
    ('(3) Level 3',    4, '7/8',                 '1.76e48', '1.42e47'),
    ('(3) Level 3',    4, '1',                   '6.30e47', '6.46e46'),
]
for part, level, label, XN, XR in claims:
    b0 = Q.parse_b0(label)
    okN, LN, nbN, mN, kN = Q.run3('N', b0, 7*LN10, arb(XN).log(), e1, level)
    okR, LR, nbR, mR, kR = Q.run3('R', b0, 10*LN10, arb(XR).log(), e1, level)
    ok_all &= okN and okR
    LX = arb(XR).log()
    eta = Q.best_eta3(LX, e1 if LX + arb('1e-3') > Q.LY1 else arb(0))     # >= |theta(X_R)/X_R - 1| (Lemma 12)
    log10n = float((arb(XR)*(1 - eta.upper())/LN10).lower())                      # log N_{eps_R} >= theta(X_R)
    print(f'{part}, b0 = {label}: Nicolas on [1e7, {XN}]: {"PASS" if okN else "FAIL"} ({nbN} boxes, min {mN:.3g}, k <= {kN}); '
          f'Robin on CA with largest prime in [1e10, {XR}]: {"PASS" if okR else "FAIL"} ({nbR} boxes, min {mR:.3g}, k <= {kR}) '
          f'=> Robin for 5041 <= n <= 10^(10^{math.floor(math.log10(log10n)*10000)/10000:.4f}); {time.time()-t0:.0f}s', flush=True)
print('ALL PASS' if ok_all else 'SOME CHECK FAILED', f'total {time.time()-t0:.0f}s')
sys.exit(0 if ok_all else 1)
