# t2_sweep.py (explore48): for lattices L_z scan the critical line from just below the smoothing height 2 pi y;
# report T1 = first t with |beta| >= 1, and the first off-line zero (deficit of sign changes of F vs zeros of the constant term,
# then complex Newton on E*). Usage: python3 t2_sweep.py <batch> ; output lines 'RES ...' (also JSON at end).
import sys, math, time, json
sys.path.insert(0, '.')
from epstein_cs import F_line, Estar

def newton(f, s0, it=60):
    s = s0; h = 1e-7
    for _ in range(it):
        fs = f(s); d = (f(s + h) - fs) / h
        if d == 0: break
        ds = fs / d
        if abs(ds) > 1.0: ds = ds / abs(ds) * 0.5
        s -= ds
        if abs(ds) < 1e-11: break
    return s, abs(f(s))

def scan(label, x, y, mult_hi=3.0, tcap=150.0, stop_after=1):
    T0 = time.time()
    ts = 2 * math.pi * y
    t = max(0.5, ts - 8 * ts ** (1 / 3))
    t_hi = mult_hi * ts + 20
    Fp, cp, bp = F_line(t, x, y, return_parts=True)[:3]
    Na = NF = 0; T1 = None; maxb = abs(bp); found = []; lastF = []  # recent F sign-change positions
    deficit_hist = []
    while t < t_hi and time.time() - T0 < tcap:
        h = math.pi / (8 * math.log(max(t * y / math.pi, 3.0)))
        tn = t + h
        Fn, cn, bn = F_line(tn, x, y, return_parts=True)[:3]
        if (cn > 0) != (cp > 0): Na += 1
        if (Fn > 0) != (Fp > 0):
            NF += 1; lastF.append(tn); lastF = lastF[-6:]
        if abs(bn) > maxb: maxb = abs(bn)
        if T1 is None and abs(bn) >= 1: T1 = tn
        d = Na - NF
        deficit_hist.append((tn, d))
        if d >= 2:
            # candidate: off-line pair between the last F sign change and now (+ a bit); try Newton at several points
            g0 = lastF[-1] if lastF else tn - 2
            ok = None
            for tc in (0.5 * (g0 + tn), g0 + 0.25 * (tn - g0), g0 + 0.75 * (tn - g0), tn):
                for dl in (0.08, 0.2, 0.4, 0.7):
                    try:
                        s, r = newton(lambda s: Estar(s, x, y), complex(0.5 + dl, tc))
                    except Exception as e:
                        continue
                    if r < 1e-9 and s.real > 0.5 + 1e-7 and abs(s.imag - tc) < 4 and s.real < 3:
                        ok = (s, r); break
                if ok: break
            if ok:
                s, r = ok
                found.append({'re': s.real, 'im': s.imag, 'resid': r, 'ratio': s.imag / ts})
                print('RES %s x=%.4f y=%.6f offline %.10f%+.10fi  t/(2pi y)=%.4f  resid=%.1e' % (label, x, y, s.real, s.imag, s.imag / ts, r), flush=True)
                Na -= 2  # account for the pair
                if len(found) >= stop_after: break
            else:
                # possibly missed close sign changes: refine [g0, tn]
                k = 200; cnt = 0; Fq = F_line(g0, x, y)
                for j in range(1, k + 1):
                    tq = g0 + (tn - g0) * j / k
                    Fr = F_line(tq, x, y)
                    if (Fr > 0) != (Fq > 0): cnt += 1
                    Fq = Fr
                # cnt includes the sign change at g0? g0 is a sign-change point (F changes between grid pts around g0); recount conservatively
                NF += max(0, cnt - 1)
                if Na - NF >= 2:
                    print('WARN %s y=%.4f unresolved deficit near t=%.3f (refined %d)' % (label, y, tn, cnt), flush=True)
                    Na -= 2
        t, Fp, cp, bp = tn, Fn, cn, bn
    el = time.time() - T0
    out = {'label': label, 'x': x, 'y': y, 'ts': ts, 'T1': T1, 'T1ratio': (T1 / ts if T1 else None), 'maxbeta': maxb,
           'scanned_to': t, 'scan_ratio': t / ts, 'offline': found, 'secs': el}
    print('SUM %s y=%.4f  2piy=%.3f  T1=%s (%.4f)  max|beta|=%.3f  scanned to %.2f (%.3f x 2piy)  first offline: %s  [%.0fs]' % (
        label, y, ts, ('%.3f' % T1) if T1 else None, (T1 / ts) if T1 else float('nan'), maxb, t, t / ts,
        ('%.4f%+.4fi ratio %.4f' % (found[0]['re'], found[0]['im'], found[0]['ratio'])) if found else 'none', el), flush=True)
    return out

batches = {
 'val': [('D-20 x2+5y2', 0.0, math.sqrt(5)), ('D-15 x2+xy+4y2', -0.5, math.sqrt(15) / 2), ('D-20 2x2+2xy+3y2', -0.5, math.sqrt(5) / 2)],
 'h1': [('D-163 h1', -0.5, math.sqrt(163) / 2), ('D-67 h1', -0.5, math.sqrt(67) / 2), ('D-43 h1', -0.5, math.sqrt(43) / 2)],
 'h2': [('D-24 x2+6y2', 0.0, math.sqrt(6)), ('D-40 x2+10y2', 0.0, math.sqrt(10)), ('D-52 x2+13y2', 0.0, math.sqrt(13)),
        ('D-88 x2+22y2', 0.0, math.sqrt(22)), ('D-148 x2+37y2', 0.0, math.sqrt(37)), ('D-232 x2+58y2', 0.0, math.sqrt(58)),
        ('D-35 x2+xy+9y2', -0.5, math.sqrt(35) / 2), ('D-51', -0.5, math.sqrt(51) / 2), ('D-91', -0.5, math.sqrt(91) / 2),
        ('D-115', -0.5, math.sqrt(115) / 2), ('D-123', -0.5, math.sqrt(123) / 2), ('D-187', -0.5, math.sqrt(187) / 2),
        ('D-235', -0.5, math.sqrt(235) / 2), ('D-267', -0.5, math.sqrt(267) / 2), ('D-403', -0.5, math.sqrt(403) / 2),
        ('D-427', -0.5, math.sqrt(427) / 2)],
 'rect': [('rect', 0.0, yy) for yy in (2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 15.0, 20.0, 25.0, 30.0, 40.0)],
 'rhomb': [('rhomb', 0.5, yy) for yy in (2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 15.0, 20.0, 25.0, 30.0, 40.0)],
 'redo': [('D-403', -0.5, math.sqrt(403) / 2), ('rect', 0.0, 15.0), ('rect', 0.0, 30.0)],
 'big': [('rect', 0.0, 50.0), ('rect', 0.0, 60.0)],
}
if __name__ == '__main__':
    b = sys.argv[1]
    tcap = float(sys.argv[2]) if len(sys.argv) > 2 else 150.0
    mult = float(sys.argv[3]) if len(sys.argv) > 3 else 3.0
    res = [scan(l, x, y, mult_hi=mult, tcap=tcap) for (l, x, y) in batches[b]]
    json.dump(res, open('logs/t2_%s.json' % b, 'w'), indent=1)
