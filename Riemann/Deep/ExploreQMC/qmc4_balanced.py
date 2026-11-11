# explore30 qmc4: balanced adversary (second-order control of S and S_1 = int S) at T = 3e12 and the needed gaps.
# usage: uv run --no-project --with numpy python3 qmc4_balanced.py
import numpy as np, math, sys, time
sys.path.insert(0, '.')
from qmc1_adversary import adversary, Sbounds, memcheck
t0 = time.time()
Sb, DT = Sbounds(3e12 + 1e4); Sallow = Sb - 1.0
S1b = 0.059 * math.log(3e12) + 2.067
print(f'# T = 3e12: D = {DT:.4f}, |S| allowed {Sallow:.3f} (Trudgian {Sb:.3f} - 1), Trudgian |S_1| <= {S1b:.3f}; d = 0.37')
for (a, wf, q) in ((1.2, 0.06, 0.8), (1.4, 0.06, 0.8), (1.85, 0.03, 0.8), (1.925, 0.03, 0.8), (1.925, 0.06, 0.8), (2.131, 0.03, 0.8), (2.303, 0.03, 0.8)):
    r = adversary(a, wf * a, DT, 0.37, q, balance=True)
    ok = r['W'] < 0 and r['Sabs'] <= Sallow and r['S1'] <= S1b - 1.0
    print(f'  a = {a}, w = {wf}a, q = {q}: W/(D2pi||f||^2) = {r["W"]:+.4f}, max|S| {r["Sabs"]:.3f}, max|int S| {r["S1"]:.3f}, '
          f'needed gap {r["gap"]:.4f}, pairs {r["npairs"]} -> {"ADMISSIBLE" if ok else "check"}', flush=True)
print(f'# maxrss {memcheck("end")/1e6:.0f} MB, {time.time()-t0:.0f}s')
