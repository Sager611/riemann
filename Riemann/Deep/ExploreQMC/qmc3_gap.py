# explore30 qmc3: needed gap at d = 0.37 for narrower bumps vs the asymptotic allowed gaps 4 pi A
# (A = 0.0946: arXiv 2606.31441 Cor. improvedzeta; A = 0.1038: HSW 2022; A = 0.112: Trudgian 2014), plus the S_1 check
# (max |int (S_mod - mean)| over the block, compared with Trudgian's |S_1(t)| <= 0.059 log t + 2.067).
# usage: uv run --no-project --with numpy python3 qmc3_gap.py
import numpy as np, math, sys, time
sys.path.insert(0, '.')
from qmc1_adversary import adversary, memcheck
t0 = time.time()
DT = math.log((3e12 + 1e4) / (2 * math.pi)) / (2 * math.pi)
lims = {A: 4 * math.pi * A for A in (0.0946, 0.1038, 0.112)}
print('# asymptotic allowed gaps 4 pi A:', {A: round(g, 4) for A, g in lims.items()})
for wf in (0.03, 0.06):
    for a in (1.7, 1.75, 1.8, 1.85, 1.925, 2.131, 2.303):
        r = adversary(a, wf * a, DT, 0.37, 1.0)
        tag = ' '.join(f'A={A}:{"ok" if (r["gap"] <= g and r["W"] < 0) else "no"}' for A, g in lims.items())
        print(f'  w = {wf}a, a = {a}: needed gap {r["gap"]:.4f}, W/(D2pi||f||^2) = {r["W"]:+.4f}, pairs {r["npairs"]}  [{tag}]', flush=True)
# narrow-bump limit (w -> 0): gap = (1 - arccos(sech(2 d a))/pi) pi / a
for a in (1.65, 1.7, 1.75, 1.77, 1.78, 1.8, 1.925, 2.131, 2.303):
    p = math.acos(1 / math.cosh(2 * 0.37 * a)) / math.pi
    print(f'  w -> 0: a = {a}: p = {p:.4f}, gap = {(1 - p) * math.pi / a:.4f}')
# S_1 check for the a = 1.925 block at T = 3e12 (rebuild S on a grid)
from qmc1_adversary import TwoBump
import qmc1_adversary as Q
a = 1.925; w = 0.06 * a; d = 0.37
r = adversary(a, w, DT, d, 1.0)
print(f'# a = 1.925 block at T = 3e12: max|S| {r["Sabs"]:.3f}, swing {r["swing"]:.3f}, max|int(S - mean)| {r["S1"]:.3f}; Trudgian S_1 bound at 3e12: {0.059*math.log(3e12)+2.067:.3f}')
print(f'# maxrss {memcheck("end")/1e6:.0f} MB, {time.time()-t0:.0f}s')
