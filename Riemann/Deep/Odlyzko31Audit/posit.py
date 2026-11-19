#!/usr/bin/env python3
"""odlaudit control: Re Phi(1/2 + a + i g) = int F(x) cosh(ax) cos(gx) dx for the certified crossover kernel
F = H(x/20.412) sech(15x/32) on a gamma grid: must be >= 0 for |a| <= 15/32 and may fail for |a| > 15/32 (sharpness of the input).
Usage: uv run --no-project --with numpy python3 posit.py > logs/posit.log"""
import numpy as np
b, w = 20.412, 15 / 32
X = 2 * b
xg, wg = np.polynomial.legendre.leggauss(400)
# 64 panels on [0, X], 400 Gauss nodes each
edges = np.linspace(0, X, 65)
xs = np.concatenate([(e1 - e0) / 2 * xg + (e1 + e0) / 2 for e0, e1 in zip(edges[:-1], edges[1:])])
ws = np.concatenate([(e1 - e0) / 2 * wg for e0, e1 in zip(edges[:-1], edges[1:])])
t = xs / b
H = (2 - t) * (1 + np.cos(np.pi * t) / 2) / 3 + np.sin(np.pi * t) / (2 * np.pi)
F = H / np.cosh(w * xs)
gam = np.linspace(0, 3, 6001)
print(f'kernel H b={b} w=15/32; gamma grid [0, 3] step {gam[1]:.4f}; Gauss-Legendre 64x400 on [0, {X}]')
for a in (0.0, 0.3, 0.46, 15 / 32, 0.475, 0.49, 0.5):
    g = F * np.cosh(a * xs) * ws
    vals = 2 * np.array([np.dot(g, np.cos(gg * xs)) for gg in gam])
    i = int(np.argmin(vals))
    scale = vals.max()
    print(f'a = {a:.5f}: min Re Phi = {vals[i]: .3e} at gamma = {gam[i]:.4f} (max {scale:.3e}); '
          f"{'NEGATIVE' if vals[i] < -1e-12 * scale else 'nonnegative (to rounding)'}")
