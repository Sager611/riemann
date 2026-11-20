#!/usr/bin/env python3
# model.py -- FLOAT model (not a certificate) of the Li-range constant c = n/T0^2 with KLN density inserted (lidens).
# Q(n) = Q_all(a*) + 2 sum_j [upper Stieltjes sum of f_j against count bound U_j] + tail at a_B beyond tK.
import sys, os, json, math, importlib.util
import numpy as np
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..', '..'))
spec = importlib.util.spec_from_file_location('li200k', os.path.join(ROOT, 'Riemann/Deep/Li200k/verify.py'))
L = importlib.util.module_from_spec(spec); spec.loader.exec_module(L)
rows = json.load(open(os.path.join(ROOT, 'Riemann/Deep/LiDens/logs/kln_rows.json')  # written by code/kln_rows.py))
T0 = 3000175332800.0; H = 3e12 - 1
PI2 = 2 * math.pi
def E(t): return 0.1038 * math.log(t) + 0.2573 * math.log(math.log(t)) + 9.4926 + 0.875
def Ntil(t): return t / PI2 * math.log(t / PI2) - t / PI2 + 0.875
def psi(n, t, a): return 4 * np.sinh(n * a / (2 * t * t)) ** 2 + 4 * n * a * a / t ** 3
def D(s, t, use=('0.85', '0.86', '0.87', '0.88'), logform=False, tag='main'):
    best = None
    for rw in rows[tag]:
        if rw['name'] not in use: continue
        A = math.exp(rw['c0'] + rw['c1'] * s); Bc = rw['q0'] - rw['q1'] * s; lt = np.log(t)
        if logform:
            x = 2 * math.pi * rw['d'] * A * lt ** (4 - 2 * s) * t ** (8 * (1 - s) / 3) / (t - H)
            v = (t - H) * lt / (2 * math.pi * rw['d']) * np.log1p(x) + Bc * lt ** 2
        else:
            v = A * lt ** (5 - 2 * s) * t ** (8 * (1 - s) / 3) + Bc * lt ** 2
        best = v if best is None else np.minimum(best, v)
    return best
def Qall(n, a, t1=None):
    # 2 psi(t1) E(t1) + int_{t1}^inf psi w dt, w = (1/2pi)log(t/2pi) + 0.13/t   (t1 = T0 default)
    t1 = T0 if t1 is None else t1
    u = np.concatenate([[0.0], np.geomspace(1e3, 200 * t1, 6000)]); t = t1 + u
    w = np.log(t / PI2) / PI2 + 0.13 / t
    y = psi(n, t, a) * w
    return 2 * psi(n, t1, a) * E(t1) + float(np.sum((y[1:] + y[:-1]) / 2 * np.diff(t)))
def Qdens(n, sig, aB, mintot=True, logform=False, use=('0.85', '0.86', '0.87', '0.88'), tag='main', tK=100 * T0, K=1500):
    u = np.concatenate([[0.0], np.geomspace(1e4, tK - T0, K)]); t = T0 + u
    Ntot = np.array([Ntil(x) for x in t]) - Ntil(T0) + np.array([E(x) for x in t]) + E(T0)
    a = np.array(sig) - 0.5
    P = np.array([psi(n, t, aa) for aa in a])            # (J+1) x K
    tot = 0.0
    for j in range(len(sig) - 1):
        f = P[j + 1] - P[j]
        U = 2 * D(sig[j], t, use, logform, tag)
        if mintot: U = np.minimum(U, Ntot)
        tot += float(np.sum(U[1:] * (f[:-1] - f[1:])) + f[-1] * U[-1])
    tail = 2 * psi(n, tK, aB) * E(tK) + Qall(n, aB, tK) - 2 * psi(n, tK, aB) * E(tK)
    return tot + tail
def Mlo(n): return float(L.M_lo(int(n), int(n)))
def cmax(Qfun, lo=100.0, hi=125.0, it=40):
    for _ in range(it):
        mid = (lo + hi) / 2; n = mid * T0 * T0
        if Qfun(n) < Mlo(n): lo = mid
        else: hi = mid
    return lo
if __name__ == '__main__':
    B = 0.87334585610671857680; aB = B - 0.5
    print('baseline Li200k model (all zeros at a_B): c = %.4f' % cmax(lambda n: Qall(n, aB)))
    for s0 in [float(x) for x in sys.argv[1:]] or [0.85]:
        for J in (40, 160):
            sig = list(np.linspace(s0, B, J + 1))
            for mt, lf in ((False, False), (True, False), (False, True), (True, True)):
                c = cmax(lambda n: Qall(n, s0 - 0.5) + Qdens(n, sig, aB, mt, lf))
                print('s*=%.3f J=%d mintot=%d logform=%d: c = %.4f' % (s0, J, mt, lf, c), flush=True)
