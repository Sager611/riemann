#!/usr/bin/env python3
# dump.py (robinchaudit) -- run RobinCH's own Level-3 search (certify3.run3, level 4) and record every box
# (exact endpoints as (man, exp), k, lam, margin lower bound, pass/fail), then probe the endpoint.
#   uv run --no-project --with python-flint --with mpmath python3 dump.py OUT.json   (run from the RobinCH copy dir)
import sys, json, time
from flint import arb
import certify3 as Q
rec = []
orig = Q.margins4
def me(a): return [int(v) for v in a.man_exp()]
def wrapped(kind, La, Lb, b0, eta1):
    m, arg = orig(kind, La, Lb, b0, eta1)
    rec.append(dict(kind=kind, b0=str(b0.mid()), La=me(La.lower()), Lb=me(Lb.upper()), La_rad=float(La.rad()), Lb_rad=float(Lb.rad()),
                    k=(arg[0] if arg else None), lam=(arg[1].hex() if arg else None), m=float(m.lower()), ok=bool(m > 0)))
    return m, arg
Q.margins4 = wrapped
e1 = Q.eta1_ch(); LN10 = Q.LN10; t0 = time.time()
runs = [('N', '0.87338196661320752', 7*LN10, '1.78e48'), ('R', '0.87338196661320752', 10*LN10, '1.43e47'),
        ('N', '1', 7*LN10, '6.30e47'), ('R', '1', 10*LN10, '6.46e46')]
summ = []
for kind, b0s, L0, X in runs:
    b0 = Q.parse_b0(b0s); n0 = len(rec)
    ok, Lend, nb, minm, km = Q.run3(kind, b0, L0, arb(X).log(), e1, 4)
    summ.append(dict(kind=kind, b0=b0s, X=X, ok=ok, nbox=nb, minm=minm, kmax=km, first=n0, last=len(rec)))
    print(kind, b0s, X, ok, nb, minm, km, f'{time.time()-t0:.0f}s', flush=True)
    # endpoint probe: continue the search from the claimed end until the certificate fails
    n1 = len(rec)
    ok2, Lf, nb2, minm2, km2 = Q.run3(kind, b0, arb(X).log(), arb(170), e1, 4)
    summ.append(dict(kind=kind, b0=b0s, X='probe', ok=ok2, Lfail_log10=float((Lf/LN10).mid()), nbox=nb2, first=n1, last=len(rec)))
    print('  probe: fails at x = 10^%.6f after %d boxes' % (float((Lf/LN10).mid()), nb2), f'{time.time()-t0:.0f}s', flush=True)
json.dump(dict(summary=summ, boxes=rec), open(sys.argv[1], 'w'))
print('dumped', len(rec), 'boxes')
