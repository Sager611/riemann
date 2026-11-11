# explore30 qmc5: balanced adversary: threshold at T = 3e12 (S and S_1 constraints) and large heights.
# usage: uv run --no-project --with numpy python3 qmc5_thresh.py
import numpy as np, math, sys, time
sys.path.insert(0, '.')
from qmc1_adversary import adversary, Sbounds, memcheck
t0 = time.time()
def run(label, t, alist, wfs=(0.03, 0.06), qs=(0.6, 0.8, 0.95)):
    Sb, D = Sbounds(t); Sallow = Sb - 1.0
    S1b = 0.059 * math.log(t) + 2.067
    A2 = 0.0946 * math.log(t) + 5.746 - 1.0                      # 2606.31441 count bound (minus margin)
    Sallow = min(Sallow, A2)
    print(f'# {label}: D = {D:.4f}, |S| allowed {Sallow:.3f}, |S_1| allowed {S1b - 1:.3f}; d = 0.37', flush=True)
    for a in alist:
        best = None
        for wf in wfs:
            for q in qs:
                r = adversary(a, wf * a, D, 0.37, q, balance=True)
                ok = r['W'] < 0 and r['Sabs'] <= Sallow and r['S1'] <= S1b - 1.0
                viol = max(r['Sabs'] - Sallow, r['S1'] - (S1b - 1.0))
                key = (not ok, viol)
                if best is None or key < best[0]:
                    best = (key, wf, q, r, ok)
        _, wf, q, r, ok = best
        print(f'  a = {a} (c = {math.exp(2*a):.2f}): w = {wf}a, q = {q}: W {r["W"]:+.4f}, max|S| {r["Sabs"]:.2f}, '
              f'max|int S| {r["S1"]:.2f}, gap {r["gap"]:.3f} -> {"ADMISSIBLE, W < 0" if ok else "not admissible"}', flush=True)
run('T = 3e12', 3e12 + 1e4, (1.4, 1.5, 1.6, 1.7, 1.8))
run('t = 1e100', 1e100, (1.85, 1.925, 2.131, 2.303), wfs=(0.03,), qs=(0.8, 0.95))
print(f'# maxrss {memcheck("end")/1e6:.0f} MB, {time.time()-t0:.0f}s')
