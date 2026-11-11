# explore30 qmc7: balanced adversary at T = 3e12 for smaller delta (beta = 1/2 + delta) at the record windows.
# usage: uv run --no-project --with numpy python3 qmc7_delta.py
import numpy as np, math, sys, time
sys.path.insert(0, '.')
from qmc1_adversary import adversary, Sbounds, memcheck
t0 = time.time()
Sb, D = Sbounds(3e12 + 1e4); Sallow = Sb - 1.0; S1a = 0.059 * math.log(3e12) + 2.067 - 1.0
print(f'# T = 3e12: |S| allowed {Sallow:.3f}, |S_1| allowed {S1a:.3f}')
for dd in (0.2, 0.1, 0.05):
    for a in (1.925, 2.131, 2.303):
        best = None
        for wf in (0.03, 0.06):
            for q in (0.6, 0.8, 0.95):
                r = adversary(a, wf * a, D, dd, q, balance=True)
                ok = r['W'] < 0 and r['Sabs'] <= Sallow and r['S1'] <= S1a
                key = (not ok, max(r['Sabs'] - Sallow, r['S1'] - S1a))
                if best is None or key < best[0]:
                    best = (key, wf, q, r, ok)
        _, wf, q, r, ok = best
        print(f'  delta = {dd}, a = {a}: w = {wf}a, q = {q}: W {r["W"]:+.4f}, max|S| {r["Sabs"]:.2f}, max|int S| {r["S1"]:.2f} '
              f'-> {"ADMISSIBLE, W < 0" if ok else "not admissible"}', flush=True)
print(f'# maxrss {memcheck("end")/1e6:.0f} MB, {time.time()-t0:.0f}s')
