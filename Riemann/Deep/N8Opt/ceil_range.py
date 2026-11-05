# ceil_range.py -- upper bounds (ceilings) for certificates of the critical line that see pairs at index distance <= s,
# from periodic trial states, as in Prop. 7.3 of arXiv 2610.08965 (floating point; the paper certifies s=7 and all pairs)
import json, sys, itertools
import numpy as np
from scipy.optimize import linprog
d = json.load(open('src/anc/certificate_line/inputs.json'))
a = np.array(d['window_coefficients'], float) / 1e6
jj = np.arange(1, len(a) + 1)
C = 4 / 3 + 0.5 * (a**2).sum() + (((-1.0) ** jj * a / jj**2 - a**2 / (4 * jj**2)).sum()) / np.pi**2
def K(t):
    s_ = np.sinc(t)
    for j, aj in enumerate(a, 1): s_ = s_ + aj / 2 * (np.sinc(t - j) + np.sinc(t + j))
    return s_
def energy(g, s):  # pair energy per period, ordered pairs, index distance <= s (s=None: all, truncated)
    g = np.asarray(g, float); p = len(g); P = g.sum()
    S = s if s is not None else 4000
    reps = S // p + 2
    y = np.concatenate([[0.0], np.cumsum(np.tile(g, reps))])
    e = 0.0
    for k in range(1, S + 1):
        e += (K(y[k:k + p] - y[:p]) ** 2).sum()
    return 2 * e, p, P
def sbar(states, s):
    E = [energy(g, s) for g in states]
    # min sum w n  s.t. sum w P = 1, sum w (e - n) <= -(2 - C), w >= 0
    c = np.array([n for e, n, P in E]); Aeq = np.array([[P for e, n, P in E]]); Aub = np.array([[e - n for e, n, P in E]])
    r = linprog(c, A_ub=Aub, b_ub=[-(2 - C)], A_eq=Aeq, b_eq=[1.0], bounds=[(0, None)] * len(E), method='highs')
    return r.fun if r.status == 0 else np.inf
st7 = [[1.03627, 1.03627, 1.96582, 1.04317, 1.96582], [1.04588, 1.97681, 1.04588, 1.98634, 1.98634]]
stall = [[1.03846, 1.96154], [1.03009, 1.03009, 1.95116, 1.03749, 1.95117]]
print('C(f) = %.10f, 2 - C = %.7f' % (C, 2 - C))
print('paper states, s=7: %.7f (paper 0.6735660)' % sbar(st7, 7))
print('paper states, all pairs: %.7f (paper 0.6738216)' % sbar(stall, None))
rng = np.random.default_rng(1)
def improve(states, s, iters=400, step=0.01):
    best = sbar(states, s); cur = [list(g) for g in states]
    while step > 1e-5:
        moved = False
        for k in range(len(cur)):
            for i in range(len(cur[k])):
                for sg in (1, -1):
                    tr = [list(g) for g in cur]; tr[k][i] += sg * step
                    v = sbar(tr, s)
                    if v < best - 1e-12: best, cur, moved = v, tr, True
        if not moved: step /= 2
    return best, cur
res = {}
for s in [7, 8, 9, 10, 11]:
    cands = []
    for base in (st7, stall):
        cands.append(improve(base, s))
    # also mix one state of each family
    cands.append(improve([st7[0], st7[1], stall[0], stall[1]], s))
    b, stt = min(cands, key=lambda z: z[0])
    res[s] = b
    print('s=%d (frames of %d points): ceiling <= %.7f  states %s' % (s, s + 1, b, [np.round(g, 4).tolist() for g in stt]), flush=True)
