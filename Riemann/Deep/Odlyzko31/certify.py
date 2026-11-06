#!/usr/bin/env python3
"""odlyzko31 CERTIFICATE (python-flint / Arb ball arithmetic).
Usage: uv run --no-project --with python-flint python3 certify.py kernels.json
For every kernel F(x) = S(x/b) sech(w x) (S = G or H of Odlyzko 1976, support |x| <= 2b) in kernels.json it encloses
  Ic = int_0^inf (1-F)/(2 sinh(x/2)) dx,  Ir = int_0^inf (1-F)/(2 cosh(x/2)) dx,  E = 4 int_0^inf F(x) cosh(x/2) dx
rigorously, and certifies, for the stated (n, r1),
  log rd_K >= [r1 logA + (n - r1) logB - E]/n,  logB = gamma + log(8 pi) - Ic,  logA = logB + pi/2 - Ir
(Poitou's explicit formula; zero terms >= 0 because every zero has |Re rho - 1/2| <= w, prime terms >= 0 because F >= 0;
see DERIVATION.md). Each printed bound is a 3-decimal number c with exp(lower bound) > c checked in Arb.
Also encloses the limit constants A_c(w) = exp(gamma + log 8 pi - Ic(sech(w x))), A_r(w) = A_c(w) exp(pi/2 - Ir(sech(w x))).
"""
import sys, json, math, time
from flint import arb, acb, fmpq, ctx
ctx.prec = 128
PI = arb.pi()
EULER = arb.const_euler()
HALF = arb(1) / 2
DELTA = arb(1) / 2 ** 24          # [0, DELTA] handled by a derivative bound

def wq(s):
    p, q = s.split('/')
    return arb(fmpq(int(p), int(q)))

def S_G(t):   # Odlyzko GRH shape on [0, 2]
    return (1 - t / 2) * (PI * t / 2).cos() + (PI * t / 2).sin() / PI

def dS_G(t):
    return -(PI / 2) * (1 - t / 2) * (PI * t / 2).sin()

def S_H(t):   # Odlyzko unconditional shape on [0, 2]
    return (2 - t) * (1 + (PI * t).cos() / 2) / 3 + (PI * t).sin() / (2 * PI)

def dS_H(t):
    return -(1 - (PI * t).cos()) / 3 - (2 - t) * PI * (PI * t).sin() / 6

SH = {'G': (S_G, dS_G), 'H': (S_H, dS_H)}

def kernel_consts(shape, b, w):
    S, dS = SH[shape]
    X = 2 * b
    F = lambda x: S(x / b) / (w * x).cosh()
    f_c = lambda x, a: (1 - F(x)) / (2 * (x / 2).sinh())
    f_r = lambda x, a: (1 - F(x)) / (2 * (x / 2).cosh())
    f_e = lambda x, a: F(x) * (x / 2).cosh()
    # [0, DELTA]: |1 - F(x)| <= x * sup|F'| and 2 sinh(x/2) >= x, so the piece lies in [-M*DELTA, M*DELTA]
    xb = arb(DELTA / 2, DELTA / 2)                       # ball containing [0, DELTA]
    dF = dS(xb / b) / b / (w * xb).cosh() - S(xb / b) * w * (w * xb).tanh() / (w * xb).cosh()
    M = abs(dF).upper()
    small = arb(0, M * DELTA)
    Ic = acb.integral(f_c, DELTA, X).real + small - (b / 2).tanh().log()
    Ir = acb.integral(f_r, 0, X).real + PI / 2 - 2 * (b / 2).tanh().atan()
    E = 4 * acb.integral(f_e, 0, X).real
    logB = EULER + (8 * PI).log() - Ic
    logA = logB + PI / 2 - Ir
    return Ic, Ir, E, logA, logB

def limit_consts(w, X=arb(200)):
    # kernel sech(w x) (the supremum of admissible kernels); tail [X, inf): 0 <= integrand <= 1/(2 sinh(x/2)) resp. 1/(2 cosh(x/2))
    f_c = lambda x, a: (1 - 1 / (w * x).cosh()) / (2 * (x / 2).sinh())
    f_r = lambda x, a: (1 - 1 / (w * x).cosh()) / (2 * (x / 2).cosh())
    small = arb(0, w * w * DELTA * DELTA / 4)            # 0 <= 1 - sech(wx) <= (wx)^2/2 on [0, DELTA]
    tc = -(X / 4).tanh().log(); tr = PI / 2 - 2 * (X / 4).tanh().atan()
    Ic = acb.integral(f_c, DELTA, X).real + small + arb(tc / 2, tc / 2)
    Ir = acb.integral(f_r, 0, X).real + arb(tr / 2, tr / 2)
    lB = EULER + (8 * PI).log() - Ic
    return (lB + PI / 2 - Ir).exp(), lB.exp(), Ic, Ir

def floor3(x):
    """largest c = k/1000 with x > c certified (x an arb ball)."""
    k = math.floor(float(x.mid()) * 1000) + 1
    while not (x > arb(fmpq(k, 1000))):
        k -= 1
    return k / 1000

def main(path):
    t0 = time.time()
    entries = json.load(open(path))
    cache, results, bad = {}, [], 0
    for e in entries:
        key = (e['wq'], e['shape'], e['b'])
        if key not in cache:
            b = arb(fmpq(*[int(v) for v in (e['b'].replace('.', '') if '.' in e['b'] else e['b'], 10 ** (len(e['b'].split('.')[1]) if '.' in e['b'] else 0))]))
            cache[key] = kernel_consts(e['shape'], b, wq(e['wq']))
        Ic, Ir, E, lA, lB = cache[key]
        n, r1 = e['n'], e['r1']
        v = (r1 * lA + (n - r1) * lB - E) / n
        rd = v.exp()
        if not rd.is_finite() or rd.rad() > 1e-9:
            bad += 1; print('BAD', e); continue
        c = floor3(rd)
        rec = dict(e, cert_rd=c, logA=float(lA.mid()), logB=float(lB.mid()), E=float(E.mid()))
        if 'target' in e:
            tgt = arb(fmpq(int(round(e['target'] * 1000)), 1000))
            ok = bool(rd > tgt)
            rec['beats_target'] = ok
            print(f"THRESHOLD {e['w']} {e['kind']} n={n} r1={r1} kernel {e['shape']} b={e['b']}: rd > {rd.str(10, radius=False)} > {e['target']}: {'PASS' if ok else 'FAIL'}")
            bad += 0 if ok else 1
        if e['kind'].startswith('cross_'):
            target = 4 * PI * (EULER.exp() if e['kind'] == 'cross_tc' else (1 + EULER).exp())
            ok = bool(rd > target)
            rec['beats_classical_limit'] = ok
            print(f"CROSSOVER {e['w']} {e['kind']} n={n} r1={r1} kernel {e['shape']} b={e['b']}: rd > {rd.str(10, radius=False)} > {target.str(10, radius=False)}: {'PASS' if ok else 'FAIL'}")
            bad += 0 if ok else 1
        results.append(rec)
    lim = {}
    for wk, ws in (('classical', '1/2'), ('ct', '15/32'), ('ab', '3733819666132076/10000000000000000'), ('liu_arxiv', '374886647040/1000000000000')):
        Ar, Ac, Ic, Ir = limit_consts(wq(ws))
        lim[wk] = dict(A_r=Ar.str(8, radius=False), A_c=Ac.str(8, radius=False), A_r_lo=math.floor(float(Ar.lower()) * 1e5) / 1e5,
                       A_c_lo=math.floor(float(Ac.lower()) * 1e5) / 1e5, A_r_hi=math.ceil(float(Ar.upper()) * 1e5) / 1e5,
                       A_c_hi=math.ceil(float(Ac.upper()) * 1e5) / 1e5, Ic=Ic.str(10), Ir=Ir.str(10))
    json.dump({'results': results, 'limits': lim}, open('certified.json', 'w'), indent=0)
    print(f'kernels {len(cache)}  entries {len(results)}  bad {bad}  time {time.time() - t0:.1f}s')
    for wk, d in lim.items():
        print('limit', wk, d)
    print('ALL CERTIFIED' if bad == 0 else 'FAILURES')

if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else 'kernels.json')
