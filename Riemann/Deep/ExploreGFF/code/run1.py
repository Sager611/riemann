# run1.py (explore17, T1): JSW / Kahane "positive piece + remainder" certificates for Weil windows.
import os, sys, time, math, numpy as np, mpmath as mp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gw import *
t0 = time.time()
N = 40; nq = 600
# 1. scale-invariant LOG piece: lambda_min(LOG on W_L) + log L should be constant (exact scaling of log|t|)
print('# 1. LOG piece (multiplier log|t|) bottom: lam(L) + log L')
for c in (1.2, 2.0, 5.0, 20.0, 100.0):
    for sec in ('even', 'odd'):
        B = Basis(c, N, sec); n = len(B.idx)
        LOG = B.jump(lambda x: 1/x, nq) - (EG + math.log(B.L))*np.eye(n)
        lam = np.linalg.eigvalsh(LOG)[0]
        print('  c=%-6g L=%.4f %s lam_min(LOG)=%.6f  lam+logL=%.6f' % (c, B.L, sec, lam, lam + math.log(B.L)), flush=True)
# 2. r0 = inf_t (phi(t) - log|t|), phi = -log pi + Re psi(1/4 + i t/2)
f = lambda t: float(-mp.log(mp.pi) + mp.re(mp.digamma(mp.mpf(1)/4 + 1j*t/2)) - mp.log(abs(t)))
ts = np.concatenate([np.linspace(0.05, 50, 4000), np.geomspace(50, 1e6, 400)])
vals = [f(t) for t in ts]; i = int(np.argmin(vals))
print('# 2. r0 = inf (phi - log|t|) = %.6f at t = %.3f; limit -log(2 pi) = %.6f' % (vals[i], ts[i], -math.log(2*math.pi)))
r0 = min(vals[i], -math.log(2*math.pi))
# 3. margins on a grid of windows
print('# 3. margins (min over sectors); M0 analytic JSW, M1 LK+V exact/pole exact/primes by norm, M2 Weyl split (LK+V) + (POLE+PR), exact = lam_min Q')
print('  c      halfw   M0        M1        M2        exact      pole-free  Qinf       M3(w_n)   M3(Galerkin norm)')
rows = []
for c in [1.1, 1.15, 1.2, 1.25, 1.3, 1.4, 1.5, 1.8, 2.0, 2.01, 2.05, 2.1, 2.27, 2.3, 2.5, 3.0, 4.0]:
    out = {}
    for sec in ('even', 'odd'):
        P = weil_parts(c, N, sec, nq); B = P['B']
        lamLOG = np.linalg.eigvalsh(P['LOG'])[0]
        polemin = min(0.0, np.linalg.eigvalsh(P['POLE'])[0])
        sw = sum(mangoldt(n)/math.sqrt(n) for n in range(2, int(math.ceil(c))) if n < c)
        M0 = lamLOG + r0 + polemin - 2*sw
        prnorm = max(abs(np.linalg.eigvalsh(P['PR'])).max(), 0.0) if sw > 0 else 0.0
        M1 = np.linalg.eigvalsh(P['LK'])[0] + P['V'][0, 0] + polemin - prnorm
        M2 = np.linalg.eigvalsh(P['LK'] + P['V'])[0] + np.linalg.eigvalsh(P['POLE'] + P['PR'])[0]
        ex = np.linalg.eigvalsh(P['LK'] + P['V'] + P['POLE'] + P['PR'])[0]
        pf = np.linalg.eigvalsh(P['LK'] + P['V'] + P['PR'])[0]
        qinf = np.linalg.eigvalsh(P['LK'] + P['V'] + P['POLE'])[0]
        M3 = qinf - sum(mangoldt(n)/math.sqrt(n) for n in range(2, int(math.ceil(c))) if n < c)
        M3g = qinf - prnorm
        out[sec] = (M0, M1, M2, ex, pf, qinf, M3, M3g)
    m = [min(out['even'][k], out['odd'][k]) for k in range(8)]
    rows.append((c, m))
    print('  %-6g %.4f  %+.4f  %+.4f  %+.4f  %+.3e  %+.4f  %+.3e  %+.4f  %+.4f' % (c, math.log(c)/2, *m), flush=True)
print('time', round(time.time() - t0, 1))
