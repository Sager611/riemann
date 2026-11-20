#!/usr/bin/env python3
"""explore58 smallest case on zeta: first M zeros.
(a) Gonek ratio at finite height: R = sum_{T1<g<T2} Z'(g)^2 / int_{T1}^{T2} Z'(t)^2 rho(t) dt, rho = log(t/2pi)/(2pi) (-> 1).
(b) model premise: periodised windows of 2n consecutive unfolded zeros -> real-rooted type-pi trig polynomial; mean d/m1 and spectrum flatness
    vs CUE (same n), and min gaps.
"""
import sys, math, json
import numpy as np
import mpmath as mp
sys.path.insert(0, 'Riemann/Deep/ExploreLK')
from lk_cue import stats, cue_angles
mp.mp.dps = 20
M = int(sys.argv[1]); n2 = int(sys.argv[2])
g = [float(mp.zetazero(k).imag) for k in range(1, M + 1)]
g = np.array(g)
zp = np.array([float(mp.siegelz(t, derivative=1)) for t in g])
T1, T2 = g[0], g[-1]
ts = np.linspace(T1, T2, 40 * M)
zpt = np.array([float(mp.siegelz(t, derivative=1)) for t in ts])
rho = np.log(ts / (2 * math.pi)) / (2 * math.pi)
f = zpt ** 2 * rho
integ = np.sum((f[1:] + f[:-1]) / 2 * np.diff(ts))
R = (np.sum(zp ** 2) - 0.5 * (zp[0] ** 2 + zp[-1] ** 2)) / integ
print(f'(a) zeta zeros 1..{M} (t in [{T1:.2f},{T2:.2f}]): Gonek ratio sum Z\'(g)^2 / int Z\'^2 rho = {R:.4f}  (asymptotically 1)', flush=True)
x = np.array([float(mp.siegeltheta(t)) / math.pi + 1 for t in g])  # unfolded
res = []
for i in range(0, M - n2, n2 // 2):
    xs = x[i:i + n2] - x[i]; span = x[i + n2] - x[i]
    xs = xs * n2 / span
    w, d, m1 = stats(xs, float(n2))
    gm = np.min(np.diff(np.concatenate([xs, [xs[0] + n2]])))
    res.append((d / m1, w, gm))
dm = np.mean([r[0] for r in res]); W = np.mean([r[1] for r in res], axis=0); gmin = np.median([r[2] for r in res])
rng = np.random.default_rng(3)
cue = [stats(cue_angles(n2, rng) * n2 / (2 * math.pi), float(n2)) for _ in range(400)]
cdm = np.mean([c[1] / c[2] for c in cue]); CW = np.mean([c[0] for c in cue], axis=0)
flat = np.array([1.0] + [2.0] * (n2 // 2)) / (n2 + 1)
print(f'(b) windows of {n2} zeta zeros: {len(res)} windows, mean d/m1 = {dm:.4f}, median min gap = {gmin:.3f}; CUE(N={n2}): mean d/m1 = {cdm:.4f}', flush=True)
print('    mean folded spectrum (zeta windows) / flat:', np.round(W / flat, 3).tolist(), flush=True)
print('    mean folded spectrum (CUE)          / flat:', np.round(CW / flat, 3).tolist(), flush=True)
json.dump(dict(M=M, n2=n2, gonek_ratio=R, zeta_dm=dm, cue_dm=cdm, zeta_spec=(W / flat).tolist(), cue_spec=(CW / flat).tolist()),
          open(f'Riemann/Deep/ExploreLK/logs/zeta_small_M{M}.json', 'w'), indent=1)
