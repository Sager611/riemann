# explore47: Bérard–Ehrenfest (CCM 1503.00955) bound for S(t) tested on Odlyzko's zeros.
# Usage: python3 bt.py validate | window <name>
import sys, os, math, time
import numpy as np
import mpmath as mp
D = os.path.dirname(os.path.abspath(__file__))
Z = os.path.join(D, 'zeros')  # Odlyzko files live in .scratch/explore47 (symlinked as zeros/)

def sieve(N):
    s = np.ones(N + 1, dtype=bool); s[:2] = False
    for p in range(2, int(N ** 0.5) + 1):
        if s[p]: s[p * p::p] = False
    return np.nonzero(s)[0]

def prime_powers(N):
    ps = sieve(N).astype(np.int64); n = [ps]; base = [ps]; q = ps.copy(); k = 1
    while True:
        k += 1; m = ps[ps.astype(float) ** k <= N]
        if len(m) == 0: break
        n.append(m ** k); base.append(m)
    n = np.concatenate(n); base = np.concatenate(base); o = np.argsort(n)
    n = n[o]; base = base[o]
    return n, np.log(base.astype(float)), base

def weights(logn, lam, n, Dl):
    # columns: for each Delta, a_n, b_n with O^± = -(1/Δ)Σ Λ n^{-1/2}[a sin ± b cos]
    A = np.zeros((len(logn), len(Dl))); B = np.zeros_like(A)
    for j, Dt in enumerate(Dl):
        eta = logn / Dt; m = eta < 1
        e = eta[m]
        A[m, j] = (1 / math.pi + (1 - e) / np.tan(math.pi * e)) * lam[m] / np.sqrt(n[m]) / Dt
        B[m, j] = (1 - e) * lam[m] / np.sqrt(n[m]) / Dt
    return A, B

def TB(Dt, logn, lam, n):
    eta = logn / Dt; m = eta < 1; e = eta[m]
    a = 1 / math.pi + (1 - e) / np.tan(math.pi * e); b = 1 - e
    return float(np.sum(lam[m] / np.sqrt(n[m]) * np.sqrt(a * a + b * b)) / Dt)

def Hpm(x, sgn):
    # Beurling's H^± at complex/real x (mpmath)
    x = mp.mpmathify(x)
    s = mp.sin(mp.pi * x)
    core = (s / mp.pi) ** 2 * (mp.psi(1, 1 - x) - mp.psi(1, 1 + x) + 2 / x)
    return core + sgn * (s / (mp.pi * x)) ** 2

def Rpm(x, t, Dt, sgn):
    d = mp.mpf(Dt) / (2 * mp.pi)
    return (Hpm(d * (t + x), sgn) + Hpm(d * (t - x), sgn)) / 2

def load(name):
    mp.mp.dps = max(mp.mp.dps, 50)  # offsets have 22 digits: never round them at the default 15 digits
    p = os.path.join(Z, name); vals = []
    with open(p) as f:
        for line in f:
            s = line.strip()
            try: float(s); vals.append(s)
            except ValueError: pass
    offs = {'zeros1': ('0', 0), 'zeros3': ('267653395647', 10 ** 12), 'zeros4': ('144176897509546973000', 10 ** 21),
            'zeros5': ('1370919909931995300000', 10 ** 22)}
    o, k0 = offs[name]
    return mp.mpf(o), k0, vals
