# t4_onset.py (explore48): first off-line zero T_off of E*(z,.) vs the smoothing height T_s = 2 pi y, with exact counting (count.py).
# For each lattice: N0 scan of F on (0, Tmax]; off-line counts N(T)-N0(T) at checkpoints r*T_s; bisection to the first off-line zero;
# complex Newton to locate it. Also T1 = first t with |beta| >= 1 (smoothing-tail criterion) and max|beta| below T_off.
import math, cmath, sys, time, json, bisect
sys.path.insert(0, '.')
from epstein_cs import F_line, Estar
from count import NE

def newton(f, s0, it=60):
    s = s0; h = 1e-7
    for _ in range(it):
        fs = f(s); d = (f(s + h) - fs) / h
        if d == 0: break
        ds = fs / d
        if abs(ds) > 0.5: ds = ds / abs(ds) * 0.25
        s -= ds
        if abs(ds) < 1e-11: break
    return s, abs(f(s))

def run(label, x, y, rmax=3.0, tcap=600.0):
    T0 = time.time(); Ts = 2 * math.pi * y; Tmax = rmax * Ts + 10
    # N0 scan + beta
    t = 0.05; Fp = F_line(t, x, y); pos = []; T1 = None; bmax_hist = []
    while t < Tmax and time.time() - T0 < tcap * 0.5:
        h = math.pi / (8 * math.log(max(t * y / math.pi, 3.0)))
        tn = t + h
        Fn, c, b = F_line(tn, x, y, return_parts=True)[:3]
        if (Fn > 0) != (Fp > 0): pos.append(tn - 0.5 * h)
        if T1 is None and abs(b) >= 1: T1 = tn
        bmax_hist.append((tn, abs(b)))
        t, Fp = tn, Fn
    Tscan = t
    def N0(T): return bisect.bisect_right(pos, T)
    def snap(T):  # move T to the midpoint between consecutive sign changes
        i = bisect.bisect_right(pos, T)
        if 0 < i < len(pos): return 0.5 * (pos[i - 1] + pos[i])
        return T
    def off(T):
        T = snap(T); v, _ = NE(T, x, y)
        return T, v - N0(T), v
    def refine(lo, hi, fac=40):
        # rescan F on [lo, hi] with step h/fac; replace sign-change positions there; return number of new sign changes
        nonlocal pos
        t = lo; Fp = F_line(t, x, y); new = []
        while t < hi:
            h = math.pi / (8 * fac * math.log(max(t * y / math.pi, 3.0)))
            tn = min(t + h, hi); Fn = F_line(tn, x, y)
            if (Fn > 0) != (Fp > 0): new.append(tn - 0.5 * h)
            t, Fp = tn, Fn
        old = [p_ for p_ in pos if lo <= p_ <= hi]
        pos = sorted([p_ for p_ in pos if not (lo <= p_ <= hi)] + new)
        return len(new) - len(old)
    cps = []
    prevT, prevO = None, 0
    first = None
    for r in (0.7, 0.85, 0.95, 1.0, 1.05, 1.1, 1.2, 1.35, 1.5, 1.75, 2.0, 2.5, 3.0):
        T = r * Ts
        if T > Tscan - 1 or time.time() - T0 > tcap: break
        Tq, o, v = off(T)
        cps.append((Tq, round(o, 3), round(v, 3)))
        if abs(o - round(o)) > 1e-3: print('WARN nonint', label, y, Tq, o, v, flush=True)
        if round(o) > 0 and first is None:
            lo, hi = (prevT if prevT else 0.5), Tq
            while hi - lo > 1.0 and time.time() - T0 < tcap:
                mid = snap(0.5 * (lo + hi))
                if not (lo < mid < hi): break
                _, om, _ = off(mid)
                if round(om) > 0: hi = mid
                else: lo = mid
            # first make sure the on-line count in [lo, hi] is not missing close pairs
            extra = refine(max(0.05, lo - 0.5), hi + 0.5)
            if extra > 0:
                _, o2, _ = off(hi)
                print('REFINE %s y=%.4f [%.2f,%.2f] +%d sign changes; offline now %.3f' % (label, y, lo, hi, extra, o2), flush=True)
                if round(o2) <= 0:
                    prevT, prevO = hi, 0
                    continue
            # Newton in [lo, hi]
            best = None
            for tq in [lo + (hi - lo) * j / 6 for j in range(7)]:
                for s0 in (0.502, 0.51, 0.52, 0.6, 0.75, 0.95, 1.3, 1.8):
                    try: s, res = newton(lambda s: Estar(s, x, y), complex(s0, tq))
                    except Exception: continue
                    if res < 1e-9 and s.real > 0.5 + 1e-6 and lo - 0.6 <= s.imag <= hi + 0.6:
                        if best is None or s.imag < best[0].imag: best = (s, res)
                if best: break
            first = {'lo': lo, 'hi': hi, 'zero': [best[0].real, best[0].imag, best[1]] if best else None}
            break
        prevT, prevO = Tq, round(o)
    tz = first['zero'][1] if first and first['zero'] else (first['hi'] if first else None)
    bmax_below = max([b for (tt, b) in bmax_hist if tz is None or tt < tz] or [0.0])
    out = {'label': label, 'x': x, 'y': y, 'Ts': Ts, 'T1': T1, 'T1r': T1 / Ts if T1 else None, 'Tscan': Tscan,
           'checkpoints': cps, 'first': first, 'Toff_r': (tz / Ts if tz else None), 'bmax_below_Toff': bmax_below,
           'secs': time.time() - T0}
    fz = first['zero'] if first else None
    print('SUM %-18s x=%5.2f y=%8.4f Ts=%8.3f T1/Ts=%s  Toff=%s  Toff/Ts=%s  max|beta|<Toff=%.3f  scanned %.2f Ts  [%.0fs]' % (
        label, x, y, Ts, ('%.4f' % (T1 / Ts)) if T1 else '-',
        ('%.6f%+.6fi' % (fz[0], fz[1])) if fz else ('bracket %.2f-%.2f' % (first['lo'], first['hi']) if first else 'none'),
        ('%.4f' % (tz / Ts)) if tz else '-', bmax_below, Tscan / Ts, out['secs']), flush=True)
    return out

if __name__ == '__main__':
    from t2_sweep import batches
    b = sys.argv[1]; rmax = float(sys.argv[2]) if len(sys.argv) > 2 else 3.0; tcap = float(sys.argv[3]) if len(sys.argv) > 3 else 600
    sel = batches[b]
    if len(sys.argv) > 4: sel = sel[int(sys.argv[4]):int(sys.argv[5])]
    res = [run(l, x, y, rmax, tcap) for (l, x, y) in sel]
    json.dump(res, open('logs/t4_%s%s.json' % (b, ('_' + sys.argv[4]) if len(sys.argv) > 4 else ''), 'w'), indent=1)
