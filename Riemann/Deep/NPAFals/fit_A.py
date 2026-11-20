#!/usr/bin/env python3
"""npafals (A): fits of the immune inside gain g(N) = minSimple(L4o) - minSimple(P) (continuous colgen, logs/A_N*_*.log).
Models: M1 g = g_inf - A/N; M2 g = g_inf - A/N + B/N^2; M3 g = A*N^-alpha (alpha > 0 means -> 0);
M4 g = A*(N-N0)*exp(-N/L) (rise then decay to 0). Prints parameters, residuals, and predictions."""
import re
import numpy as np
from scipy.optimize import curve_fit

def val(N, S):
    try: s = open(f'logs/A_N{N}_{S}.log').read()
    except FileNotFoundError: return None
    m = re.findall(r'FINAL N=\d+ S=\w+ minSimple=([0-9.]+) lb~([0-9.]+)', s)
    if m: return float(m[-1][0]), float(m[-1][1]), True
    m = re.findall(r'minSimple=([0-9.]+) slack=([0-9.e+-]+) maxred=[0-9.e+-]+ lb~(-?[0-9.]+)', s)
    if m and float(m[-1][1]) < 1e-9: return float(m[-1][0]), float(m[-1][2]), False
    return None

rows = []
for N in range(3, 13):
    p, o = val(N, 'P'), val(N, 'L4o')
    if p and o:
        rows.append((N, o[0] - p[0]))
        print(f"N={N} P={p[0]:.6f} L4o={o[0]:.6f} (lb {o[1]:.6f}{'' if o[2] else ', unfinished'}) inside gain={o[0]-p[0]:.6f} "
              f"[lb-based {o[1]-p[0]:.6f}]")
    full, oo = val(N, 'L4'), val(N, 'L4oo')
    if p and full: print(f"   full L4={full[0]:.6f} gain={full[0]-p[0]:.6f}")
    if p and oo: print(f"   L4oo (cutoff N-2, theta={(N-2)/N:.3f}) gain={oo[0]-p[0]:.6f}")
if len(rows) >= 3:
    N = np.array([r[0] for r in rows], float); g = np.array([r[1] for r in rows])
    models = {
        'M1 g_inf - A/N': (lambda n, gi, a: gi - a / n, [0.003, 0.01]),
        'M2 g_inf - A/N + B/N^2': (lambda n, gi, a, b: gi - a / n + b / n**2, [0.003, 0.01, 0.0]),
        'M3 A N^-alpha': (lambda n, a, al: a * n**(-al), [0.002, -1.0]),
        'M4 A (N-N0) exp(-N/L)': (lambda n, a, n0, L: a * (n - n0) * np.exp(-n / L), [0.001, 3.0, 10.0]),
    }
    for name, (f, p0) in models.items():
        if len(p0) >= len(N): print(f"{name}: too few points"); continue
        try:
            par, _ = curve_fit(f, N, g, p0=p0, maxfev=20000)
            res = g - f(N, *par)
            pred = {n: f(n, *par) for n in (10, 12, 16, 24, 50, 1000)}
            print(f"{name}: params={np.round(par, 6).tolist()} max|res|={np.abs(res).max():.2e} res={np.round(res, 6).tolist()} "
                  f"pred " + " ".join(f"N={n}:{v:.5f}" for n, v in pred.items()))
        except Exception as e:
            print(f"{name}: fit failed ({e})")
