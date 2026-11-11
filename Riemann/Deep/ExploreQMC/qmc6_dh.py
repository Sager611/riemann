# explore30 qmc6: localized two-bump test functions on the Davenport-Heilbronn zeros (ExploreHelio dh_zeros.json, T <= 420):
# is the zero-side form min over (theta, centre, w) negative near DH's off-line pairs?  (explore19: DH's window form fails
# at c ~ 30.7 through low-frequency near-null modes; localized f are expected to miss isolated pairs.)
# usage: uv run --no-project --with numpy python3 qmc6_dh.py
import numpy as np, math, sys, time, json
sys.path.insert(0, '.')
from qmc1_adversary import TwoBump, memcheck
t0 = time.time()
d = json.load(open('../ExploreHelio/results/dh_zeros.json'))
on = np.array(d['online']); off = np.array(d['offline'])          # off: (beta, gamma), beta > 1/2
on_all = np.concatenate([on, -on])
off_g = np.concatenate([off[:, 1], -off[:, 1]]); off_d = np.concatenate([off[:, 0] - 0.5] * 2)
print(f'# DH: {len(on)} on-line, {len(off)} off-line pairs (delta {", ".join(f"{x:.3f}" for x in off[:,0]-0.5)})')
for a in (1.5, 1.71, 1.925, 2.303):
    for wf in (0.15, 0.3, 0.45):
        f0 = TwoBump(a, wf * a); w = wf * a; Umax = 60.0 / w
        best = None
        for (b, gp) in off[:4]:
            for tc in gp + np.linspace(-3, 3, 25):
                if tc + Umax > 420: continue
                for th in np.linspace(0, 2 * math.pi, 37)[:-1]:
                    f = TwoBump(a, w, th)
                    uo = on_all - tc; uo = uo[np.abs(uo) <= Umax]
                    uf = off_g - tc; m = np.abs(uf) <= Umax
                    W = float(np.sum(np.abs(f.F(uo))**2)) + 2 * float(np.sum(f.h(uf[m], off_d[m]))) if m.any() else float(np.sum(np.abs(f.F(uo))**2))
                    W /= 2 * math.pi * f.norm2
                    if best is None or W < best[0]:
                        best = (W, gp, tc, th)
        if best is None:
            print(f'  a = {a}, w = {wf}a: window too wide for T <= 420'); continue
        print(f'  a = {a}, w = {wf}a: min W/(2pi||f||^2) = {best[0]:+.4f} (pair at {best[1]:.2f}, centre {best[2]:.2f}, theta {best[3]:.2f})', flush=True)
print(f'# maxrss {memcheck("end")/1e6:.0f} MB, {time.time()-t0:.0f}s')
