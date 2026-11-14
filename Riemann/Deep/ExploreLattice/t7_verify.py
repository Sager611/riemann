# t7_verify.py (explore48): verify every birth (jump of the off-line count) found by t6_dist.py:
# exact zero count in (Ta, Tb] from NE, on-line count by a fine scan (step = base/fac), Newton localisation of off-line zeros.
import math, sys, json, time
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

def online(x, y, a, b, fac):
    t = a; Fp = F_line(t, x, y); n = 0; mins = []; Fpp = None
    while t < b:
        h = math.pi / (8 * fac * math.log(max(t * y / math.pi, 3.0)))
        tn = min(t + h, b); Fn = F_line(tn, x, y)
        if (Fn > 0) != (Fp > 0): n += 1
        if Fpp is not None and abs(Fp) < abs(Fpp) and abs(Fp) < abs(Fn): mins.append((abs(Fp), t))
        Fpp, Fp, t = Fp, Fn, tn
    return n, sorted(mins)[:6]

def verify(path, fac=50):
    d = json.load(open(path)); x, y, Ts = d['x'], d['y'], d['Ts']
    rows = d['rows']; out = []
    prev = rows[0]
    for r in rows[1:]:
        if round(r[2]) > round(prev[2]):
            Ta, Tb = prev[1], r[1]
            na, _ = NE(Ta, x, y); nb, _ = NE(Tb, x, y)
            n_on, mins = online(x, y, Ta, Tb, fac)
            n_off = round(nb - na) - n_on
            zeros = []
            if n_off > 0:
                for tq in [Ta + (Tb - Ta) * j / 8 for j in range(9)] + [m[1] for m in mins]:
                    for s0 in (0.502, 0.51, 0.53, 0.6, 0.75, 0.95, 1.3):
                        try: s, res = newton(lambda s: Estar(s, x, y), complex(s0, tq))
                        except Exception: continue
                        if res < 1e-9 and s.real > 0.5 + 1e-6 and Ta - 1e-6 <= s.imag <= Tb + 1e-6:
                            if all(abs(s - z) > 1e-6 for z in zeros): zeros.append(s)
                    if len(zeros) * 2 >= n_off: break
            print('%s  (%.3f, %.3f] T/Ts (%.4f, %.4f]: zeros %d, on-line %d, off-line %d, located %s' % (
                d['label'], Ta, Tb, Ta / Ts, Tb / Ts, round(nb - na), n_on, n_off,
                ', '.join('%.6f%+.6fi (%.4f T_s)' % (z.real, z.imag, z.imag / Ts) for z in zeros) or '-'), flush=True)
            out.append({'Ta': Ta, 'Tb': Tb, 'zeros': round(nb - na), 'online': n_on, 'offline': n_off,
                        'located': [[z.real, z.imag] for z in zeros]})
        prev = r
    return out

if __name__ == '__main__':
    k = sys.argv[1]; fac = int(sys.argv[2]) if len(sys.argv) > 2 else 50
    res = verify('logs/t6_%s.json' % k, fac)
    json.dump(res, open('logs/t7_%s.json' % k, 'w'), indent=1)
