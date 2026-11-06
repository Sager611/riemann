#!/usr/bin/env python3
"""odlyzko31 float layer: Poitou explicit-formula lower bounds for root discriminants.
Kernel F(x) = S(x/b) * sech(wp*x) on [0, 2b] (0 beyond), S in {G, H} (Odlyzko 1976 shapes).
log D >= r1*logA + 2*r2*logB - E  (prime terms >= 0 and zero terms >= 0 dropped), where
  Ic = int_0^inf (1-F)/(2 sinh(x/2)),  Ir = int_0^inf (1-F)/(2 cosh(x/2)),  E = 4 int_0^inf F cosh(x/2),
  logB = gamma + log(8 pi) - Ic,  logA = logB + pi/2 - Ir.
Admissible iff wp >= w, where all zeros satisfy |Re rho - 1/2| <= w (w = 1/2 unconditional, 15/32 Cureton-Thorner,
B - 1/2 abelian (Bebel), 0 GRH)."""
import math, sys
import numpy as np
EG = 0.57721566490153286
L8P = math.log(8 * math.pi)

def G(t):
    return (1 - t / 2) * np.cos(np.pi * t / 2) + np.sin(np.pi * t / 2) / np.pi

def H(t):
    return (2 - t) * (1 + np.cos(np.pi * t) / 2) / 3 + np.sin(np.pi * t) / (2 * np.pi)

SHAPES = {'G': G, 'H': H}

def sech(y):
    e = np.exp(-np.abs(y))
    return 2 * e / (1 + e * e)
_T, _W = np.polynomial.legendre.leggauss(40)

def nodes(X):
    br = [0.0]
    x = 0.125
    while x < X:
        br.append(x); x *= 2
    br.append(X)
    xs, ws = [], []
    for a, c in zip(br[:-1], br[1:]):
        xs.append((c - a) / 2 * _T + (c + a) / 2); ws.append((c - a) / 2 * _W)
    return np.concatenate(xs), np.concatenate(ws)

def consts(shape, b, wp):
    X = 2 * b
    x, wt = nodes(X)
    F = SHAPES[shape](x / b) * sech(wp * x)
    Ic = float(np.sum(wt * (1 - F) / (2 * np.sinh(x / 2)))) - math.log(math.tanh(X / 4))
    Ir = float(np.sum(wt * (1 - F) / (2 * np.cosh(x / 2)))) + math.pi / 2 - 2 * math.atan(math.tanh(X / 4))
    E = 4 * float(np.sum(wt * F * np.cosh(x / 2)))
    logB = EG + L8P - Ic
    return logB + math.pi / 2 - Ir, logB, E

def limit_consts(w, X=400.0):
    # F = sech(w x): the pointwise supremum of admissible kernels
    x, wt = nodes(X)
    F = 1 / np.cosh(w * x)
    Ic = float(np.sum(wt * (1 - F) / (2 * np.sinh(x / 2))))
    Ir = float(np.sum(wt * (1 - F) / (2 * np.cosh(x / 2))))
    logB = EG + L8P - Ic
    return math.exp(logB + math.pi / 2 - Ir), math.exp(logB), Ic, Ir

if __name__ == '__main__':
    # reproduce Odlyzko 1976 Tables 3 (GRH, G, wp=0) and 4 (unconditional, H, wp=1/2)
    print('Table 3 check (GRH): b, A, B, E')
    for b in (0.3, 0.42, 0.85, 1.5, 5.0, 20.0, 70.0):
        lA, lB, E = consts('G', b, 0.0)
        print(f'  {b:7.3f} {math.exp(lA):9.3f} {math.exp(lB):8.3f} {E:.4e}  E_formula={8*b*math.pi**2*((2*math.cosh(b/2))/(math.pi**2+b**2))**2:.4e}')
    print('Table 4 check (unconditional): b, A, B, E')
    for b in (0.3, 0.42, 0.85, 1.5, 5.0, 20.0, 75.0):
        lA, lB, E = consts('H', b, 0.5)
        print(f'  {b:7.3f} {math.exp(lA):9.3f} {math.exp(lB):8.3f} {E:.5e}  8b/3={8*b/3:.5e}')
    for w in (0.5, 15/32, 0.3733819666132076, 0.0):
        A, B, Ic, Ir = limit_consts(w) if w > 0 else (None, None, 0, 0)
        if w > 0:
            print(f'limit w={w:.10f}: A_r={A:.5f} A_c={B:.5f} Ic={Ic:.8f} Ir={Ir:.8f}')
    print('GRH limits', math.exp(EG + L8P + math.pi/2), math.exp(EG + L8P))
    print('classical', 4*math.pi*math.exp(1+EG), 4*math.pi*math.exp(EG))
