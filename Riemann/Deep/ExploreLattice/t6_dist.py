# t6_dist.py (explore48): where are the off-line zeros born? Exact off-line count N(T)-N0(T) at checkpoints every 0.05 T_s
# (T_s = 2 pi y) up to R T_s, with local refinement of N0 when the count is odd/negative/decreasing. Births = jumps of the count.
import math, sys, json, bisect, time
sys.path.insert(0, '.')
from epstein_cs import F_line
from count import NE

def dist(label, x, y, R=5.0, step=0.05, t_start=0.05):
    T0 = time.time(); Ts = 2 * math.pi * y; Tmax = R * Ts + 5
    t = t_start; Fp = F_line(t, x, y); pos = []
    while t < Tmax:
        h = math.pi / (8 * math.log(max(t * y / math.pi, 3.0)))
        tn = t + h; Fn = F_line(tn, x, y)
        if (Fn > 0) != (Fp > 0): pos.append(tn - 0.5 * h)
        t, Fp = tn, Fn
    def refine(lo, hi, fac=40):
        nonlocal pos
        t = lo; Fp = F_line(t, x, y); new = []
        while t < hi:
            h = math.pi / (8 * fac * math.log(max(t * y / math.pi, 3.0)))
            tn = min(t + h, hi); Fn = F_line(tn, x, y)
            if (Fn > 0) != (Fp > 0): new.append(tn - 0.5 * h)
            t, Fp = tn, Fn
        old = [p for p in pos if lo <= p <= hi]
        pos = sorted([p for p in pos if not (lo <= p <= hi)] + new)
        return len(new) - len(old)
    def snap(T):
        i = bisect.bisect_right(pos, T)
        return 0.5 * (pos[i - 1] + pos[i]) if 0 < i < len(pos) else T
    rows = []; prev = 0; prevT = 0.05
    r = 0.9
    while r <= R + 1e-9:
        T = snap(r * Ts); v, _ = NE(T, x, y); o = v - bisect.bisect_right(pos, T)
        if round(o) % 2 == 1 or round(o) < prev:
            ex = refine(max(0.05, prevT - 0.3), T + 0.3)
            T = snap(r * Ts); v, _ = NE(T, x, y); o = v - bisect.bisect_right(pos, T)
        rows.append((round(T / Ts, 4), round(T, 3), round(o, 3)))
        if round(o) != prev: print('JUMP %s y=%.3f  T/Ts in (%.3f, %.3f]  offline %d -> %d' % (label, y, prevT / Ts, T / Ts, prev, round(o)), flush=True)
        prev, prevT = round(o), T
        r = round(r + step, 6)
    print('DIST %s y=%.4f: final offline count %d up to %.1f (%.2f T_s); total zeros %.0f  [%.0fs]' % (label, y, prev, prevT, prevT / Ts, v, time.time() - T0), flush=True)
    return {'label': label, 'x': x, 'y': y, 'Ts': Ts, 'rows': rows}

if __name__ == '__main__':
    cases = {'rect10': ('rect', 0.0, 10.0), 'D427': ('D-427', -0.5, math.sqrt(427) / 2), 'D163': ('D-163 h1', -0.5, math.sqrt(163) / 2),
             'rect20': ('rect', 0.0, 20.0)}
    k = sys.argv[1]; R = float(sys.argv[2]) if len(sys.argv) > 2 else 5.0
    res = dist(*cases[k], R=R)
    json.dump(res, open('logs/t6_%s.json' % k, 'w'), indent=1)
