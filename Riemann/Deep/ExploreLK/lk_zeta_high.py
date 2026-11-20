#!/usr/bin/env python3
"""explore58: model premise at height 1e8 and 1e11 (ExploreDPP zero windows, offsets from tc).
Periodised windows of 2n consecutive zeros -> real-rooted type-pi trig polynomial; mean d/m1 and spectrum vs flat, vs CUE."""
import sys, math, json, re
import numpy as np
sys.path.insert(0, 'Riemann/Deep/ExploreLK')
from lk_cue import stats, cue_angles
n2 = int(sys.argv[1]); out = {}
for tag in ['1e8', '1e11']:
    tc = {}; runs = []
    for line in open(f'Riemann/Deep/ExploreDPP/logs/zeros_{tag}.txt'):
        if line.startswith('ZG'):
            m = re.match(r'ZG n=(\d+) tc=([0-9.eE+-]+)', line); tc[m.group(1)] = float(m.group(2))
        elif line.startswith('ZZ'):
            parts = line.split(); key = parts[1][2:]
            offs = np.array([float(p) for p in parts[2:]])
            runs.append((tc.get(key), offs))
    res = []
    for t0, offs in runs:
        if t0 is None or len(offs) < n2 + 1: continue
        L = math.log(t0 / (2 * math.pi)); x = np.sort(offs) * L / (2 * math.pi)
        for i in range(0, len(x) - n2, n2 // 2):
            xs = x[i:i + n2] - x[i]; span = x[i + n2] - x[i]; xs = xs * n2 / span
            w, d, m1 = stats(xs, float(n2)); res.append((d / m1, w, np.min(np.diff(np.concatenate([xs, [n2]])))))
    flat = np.array([1.0] + [2.0] * (n2 // 2)) / (n2 + 1)
    dm = float(np.mean([r[0] for r in res])); W = np.mean([r[1] for r in res], axis=0)
    print(f'height {tag}: {len(res)} windows of {n2} zeros; mean d/m1 = {dm:.4f}; median min gap {np.median([r[2] for r in res]):.3f}; spectrum/flat = {np.round(W/flat,2).tolist()}', flush=True)
    out[tag] = dict(windows=len(res), dm=dm, spec=(W / flat).tolist())
rng = np.random.default_rng(5)
cue = [stats(cue_angles(n2, rng) * n2 / (2 * math.pi), float(n2)) for _ in range(1000)]
CW = np.mean([c[0] for c in cue], axis=0); flat = np.array([1.0] + [2.0] * (n2 // 2)) / (n2 + 1)
print(f'CUE N={n2}: mean d/m1 = {np.mean([c[1]/c[2] for c in cue]):.4f}; spectrum/flat = {np.round(CW/flat,2).tolist()}', flush=True)
out['cue'] = dict(dm=float(np.mean([c[1] / c[2] for c in cue])), spec=(CW / flat).tolist())
json.dump(out, open(f'Riemann/Deep/ExploreLK/logs/zeta_high_n{n2}.json', 'w'), indent=1)
