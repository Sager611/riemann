"""Build an explicit periodic dual configuration for the TEN window: blocks '.xx' (a copies) and '.xx.xx.xxx' (b copies)
in one period, relax positions, report e vs nu - H_W. Writes dual_TEN.json (rational positions)."""
import sys, json
from fractions import Fraction
import numpy as np
from scipy.optimize import minimize
from tenwin import TEN_C, MT_C, H_of
from periodic import setup, energy

def build(a, b, spread=True):
    blocks = []
    if spread and b > 0:
        # distribute the b long blocks evenly among the a short ones
        per = a // b
        rem = a - per * b
        for i in range(b):
            blocks += ['.xx'] * (per + (1 if i < rem else 0)) + ['.xx.xx.xxx']
    else:
        blocks = ['.xx'] * a + ['.xx.xx.xxx'] * b
    s = ''.join(blocks)
    x0 = np.array([i for i, ch in enumerate(s) if ch == 'x'], dtype=float)
    return s, x0

def run(c, a, b):
    H = H_of(c)
    s, x0 = build(a, b)
    S = len(s); p = len(x0)
    ms, kh, k0 = setup(c, S)
    # shift the integer seeds toward the relaxed '.xx' geometry (0.98, 2.02) before relaxing
    r = minimize(lambda x: energy(x, S, ms, kh, k0), x0, jac=True, method="L-BFGS-B",
                 options={"maxiter": 20000, "gtol": 1e-15, "ftol": 1e-17})
    x = np.sort(np.mod(r.x, S))
    e = energy(x, S, ms, kh, k0)[0]
    nu = p / S
    return dict(a=a, b=b, S=S, p=p, nu=nu, e=e, slack=nu - H - e, H=H, x=x.tolist(), pattern=s)

if __name__ == "__main__":
    c = TEN_C if (len(sys.argv) < 2 or sys.argv[1] == "TEN") else MT_C
    res = []
    for (a, b) in [(12, 1), (11, 1), (10, 1), (13, 1), (14, 1), (23, 2), (25, 2), (35, 3), (9, 1)]:
        d = run(c, a, b)
        res.append(d)
        print(f"a={a} b={b} S={d['S']} p={d['p']} nu={d['nu']:.7f} e={d['e']:.7e} nu-H-e={d['slack']:.3e}", flush=True)
    ok = [d for d in res if d['slack'] > 0]
    best = min(ok, key=lambda d: d['nu']) if ok else None
    if best:
        print("best dual: nu =", best['nu'], "slack", best['slack'])
        json.dump(best, open("dual_TEN.json", "w"))
