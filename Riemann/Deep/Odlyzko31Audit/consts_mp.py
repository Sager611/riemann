#!/usr/bin/env python3
"""odlaudit: independent mpmath (256-bit) recomputation of the odlyzko31 constants, thresholds and random table entries.
Own code (not imported from Riemann/Deep/Odlyzko31). Usage:
  uv run --no-project --with mpmath python3 consts_mp.py [nrandom] > logs/consts_mp.log
"""
import json, random, sys, time
import mpmath as mp
mp.mp.prec = 256
pi, EG = mp.pi, mp.euler
SRC = '/Users/adriansager/playground/riemann/Riemann/Deep/Odlyzko31/'

def SG(t): return (1 - t / 2) * mp.cos(pi * t / 2) + mp.sin(pi * t / 2) / pi
def SH(t): return (2 - t) * (1 + mp.cos(pi * t) / 2) / 3 + mp.sin(pi * t) / (2 * pi)
SHAPE = {'G': SG, 'H': SH}

def rat(s):
    if '/' in s:
        p, q = s.split('/'); return mp.mpf(int(p)) / int(q)
    return mp.mpf(s)

def kconsts(shape, b, w):
    """I_c, I_r, E for F(x) = S(x/b) sech(w x), support [-2b, 2b]; tails integrated numerically (independent of closed forms)."""
    S = SHAPE[shape]; X = 2 * b
    F = lambda x: S(x / b) / mp.cosh(w * x)
    N = int(mp.ceil(X / 3)) + 1
    pts = [X * k / N for k in range(N + 1)]
    Ic = mp.quad(lambda x: (1 - F(x)) / (2 * mp.sinh(x / 2)), pts) + mp.quad(lambda x: 1 / (2 * mp.sinh(x / 2)), [X, X + 20, mp.inf])
    Ir = mp.quad(lambda x: (1 - F(x)) / (2 * mp.cosh(x / 2)), pts) + mp.quad(lambda x: 1 / (2 * mp.cosh(x / 2)), [X, X + 20, mp.inf])
    E = 4 * mp.quad(lambda x: F(x) * mp.cosh(x / 2), pts)
    return Ic, Ir, E

def limit(w):
    pts = [0, 1, 2, 4, 8, 16, 32, 64, 128, 256, mp.inf]
    Ic = mp.quad(lambda x: (1 - 1 / mp.cosh(w * x)) / (2 * mp.sinh(x / 2)), pts)
    Ir = mp.quad(lambda x: (1 - 1 / mp.cosh(w * x)) / (2 * mp.cosh(x / 2)), pts)
    lB = EG + mp.log(8 * pi) - Ic
    return Ic, Ir, mp.exp(lB), mp.exp(lB + pi / 2 - Ir)

def bound(Ic, Ir, E, n, r1):
    lB = EG + mp.log(8 * pi) - Ic; lA = lB + pi / 2 - Ir
    return lA, lB, mp.exp((r1 * lA + (n - r1) * lB - E) / n)

def main():
    nr = int(sys.argv[1]) if len(sys.argv) > 1 else 40
    t0 = time.time()
    print('## limits (mpmath 256-bit, [0,inf) split at powers of 2)')
    for name, ws in (('classical', '1/2'), ('ct', '15/32'), ('ab', '3733819666132076/10000000000000000'), ('liu', '374886647040/1000000000000'), ('grh', '0')):
        Ic, Ir, Ac, Ar = limit(rat(ws))
        print(f'limit {name} w={ws}: Ic={mp.nstr(Ic, 15)} Ir={mp.nstr(Ir, 15)} A_c={mp.nstr(Ac, 12)} A_r={mp.nstr(Ar, 12)}')
    print('closed forms w=1/2: log2 =', mp.nstr(mp.log(2), 15), ' pi/2-1 =', mp.nstr(pi / 2 - 1, 15),
          ' 4pi e^g =', mp.nstr(4 * pi * mp.exp(EG), 12), ' 4pi e^(1+g) =', mp.nstr(4 * pi * mp.exp(1 + EG), 12),
          ' GRH 8pi e^g =', mp.nstr(8 * pi * mp.exp(EG), 12), ' 8pi e^(g+pi/2) =', mp.nstr(8 * pi * mp.exp(EG + pi / 2), 12))
    print('## Lemma A: 2 int_0^inf cosh(ax)cos(tx)/cosh(wx) dx vs (2pi/w)cos(pi a/2w)cosh(pi t/2w)/(cosh(pi t/w)+cos(pi a/w))')
    for a, w, t in ((0.3, 15 / 32, 0), (0.3, 15 / 32, 2.5), (0.46, 15 / 32, 0.7), (-0.2, 0.3733819666, 1.3), (0.0, 0.25, 4.0)):
        a, w, t = mp.mpf(a), mp.mpf(w), mp.mpf(t)
        f = lambda x: 2 * mp.cosh(a * x) * mp.cos(t * x) / mp.cosh(w * x)
        num = mp.quadosc(f, [0, mp.inf], omega=t) if t > 0 else mp.quad(f, [0, 1, 4, 16, 64, 256, 1024, 4096, mp.inf])
        cf = (2 * pi / w) * mp.cos(pi * a / (2 * w)) * mp.cosh(pi * t / (2 * w)) / (mp.cosh(pi * t / w) + mp.cos(pi * a / w))
        print(f'  a={mp.nstr(a,4)} w={mp.nstr(w,6)} t={mp.nstr(t,3)}: quad {mp.nstr(num, 15)} closed {mp.nstr(cf, 15)} diff {mp.nstr(num - cf, 3)}')
    print('## Lemma B: int_{-2}^{2} S(x) cos(tx) dx vs |hhat(t)|^2/||h||^2')
    for sh, h, nh in (('G', lambda u: mp.cos(pi * u / 2), 1), ('H', lambda u: 1 + mp.cos(pi * u), 3)):
        for t in (0.0, 1.7, 5.3, 12.0):
            t = mp.mpf(t)
            ft = 2 * mp.quad(lambda x: SHAPE[sh](x) * mp.cos(t * x), [0, 0.5, 1, 1.5, 2])
            hh = 2 * mp.quad(lambda u: h(u) * mp.cos(t * u), [0, 0.5, 1])
            print(f'  {sh} t={mp.nstr(t,3)}: FT {mp.nstr(ft, 15)} |hhat|^2/|h|^2 {mp.nstr(hh * hh / nh, 15)} diff {mp.nstr(ft - hh * hh / nh, 3)}')
    cert = json.load(open(SRC + 'certified.json'))['results']
    special = [e for e in cert if e['kind'] != 'table']
    table = [e for e in cert if e['kind'] == 'table']
    random.seed(20261119)
    pick = special + random.sample(table, nr)
    print(f'## entries: {len(special)} special + {nr} random table entries (seed 20261119)')
    worst = 0; fails = 0
    for e in pick:
        b = mp.mpf(e['b']); w = rat(e['wq'])
        Ic, Ir, E = kconsts(e['shape'], b, w)
        lA, lB, rd = bound(Ic, Ir, E, e['n'], e['r1'])
        rel = max(abs(lA / e['logA'] - 1), abs(lB / e['logB'] - 1), abs(E / e['E'] - 1))
        worst = max(worst, rel)
        ok = rd > e['cert_rd']
        tgt = ''
        if 'target' in e: ok = ok and rd > e['target']; tgt = f" target {e['target']}"
        if e['kind'].startswith('cross_'):
            lim = 4 * pi * mp.exp(EG) if e['kind'] == 'cross_tc' else 4 * pi * mp.exp(1 + EG)
            ok = ok and rd > lim; tgt = f' limit {mp.nstr(lim, 12)}'
        fails += 0 if ok else 1
        print(f"{e['w']:9s} {e['kind']:8s} n={e['n']:7d} r1={e['r1']:7d} {e['shape']} b={e['b']:8s}: rd>={mp.nstr(rd, 12)} cert {e['cert_rd']}{tgt} rel(logA,logB,E)={float(rel):.1e} {'OK' if ok else 'FAIL'}", flush=True)
    print(f'entries {len(pick)} fails {fails} worst rel diff vs certified.json floats {float(worst):.2e} time {time.time() - t0:.0f}s')

if __name__ == '__main__':
    main()
