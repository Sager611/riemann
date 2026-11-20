#!/usr/bin/env python3
# kln_rows.py -- KLN Lemma 4.14 rows at H0 = 3e12 from Ingham4/cert.py build_kln_H0 (read-only import); prints A(s), B(s) data.
import sys, os, json, importlib.util
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..', '..'))
spec = importlib.util.spec_from_file_location('ing4cert', os.path.join(ROOT, 'Riemann/Deep/Ingham4/cert.py'))
C = importlib.util.module_from_spec(spec); spec.loader.exec_module(C)
from flint import arb, ctx
ctx.prec = 160
out = {}
for tag, de in (('main', None), ('delta1', '1')):
    P = C.build_kln_H0('3e12', de)
    rows = []
    for p in P:
        rows.append(dict(name=p.name, k=float(p.k), d=float(p.d), c0=float(p.c0), c1=float(p.c1), q0=float(p.q0), q1=float(p.q1),
                         smin=float(p.smin), A_s0=float(p.logA(p.s0).exp()), B_s0=float(p.Bc(p.s0)), ok=bool(p.ok)))
    out[tag] = rows
    for rw in rows:
        if rw['name'] in ('0.80', '0.85', '0.86', '0.87', '0.88', '0.89'):
            print(tag, rw['name'], 'k=%.3f d=%.3f c1=%.4f A(s0)=%.5f B(s0)=%.5f smin=%.4f ok=%s' % (rw['k'], rw['d'], rw['c1'], rw['A_s0'], rw['B_s0'], rw['smin'], rw['ok']))
json.dump(out, open(os.path.join(ROOT, 'Riemann/Deep/LiDens/logs/kln_rows.json'), 'w'), indent=1)
