#!/usr/bin/env python3
# k_of_L.py -- K as a function of the range: E*(L) = majorant(L) (cert.py, same envelopes) for several L, and the least L
# (to 0.01) giving K+1, for KLN (RH 3e12) and CS (b4 corrected). Informational (not part of the certificate).
#   uv run --no-project --with python-flint python3 Riemann/Deep/Robin2Adic/k_of_L.py
import os, sys, importlib.util
sys.dont_write_bytecode = True                             # no __pycache__ next to the imported files
HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location('r2a_cert', os.path.join(HERE, 'cert.py'))
C = importlib.util.module_from_spec(spec); spec.loader.exec_module(C)
from flint import arb
C.KLN_P = C.I4.build_kln_H0('3e12')
C.CS_R = []
for row in C.I4.CS_ROWS:
    r = C.I4.CSRow(row); C.CS_R.append((r.s1, r.s2, arb(row[2]), r.B1, arb(row[3]), arb(row[4])))
G = C.grid()
for name, Wf in [('KLN_RH3e12', C.W_kln), ('CS_b4_corrected', C.W_cs)]:
    c, Wt = C.envelope(G, Wf)
    row = []
    for Lv in ['80', '90', '100', '108.57', '110', '115', '120', '130', '150']:
        E = C.majorant(arb(Lv), G, c); lg, K = C.kval(E); row.append('L=%s: E*<=%s K=%d' % (Lv, E.str(4, radius=False), K))
    print(name + ': ' + '; '.join(row), flush=True)
    E1 = C.majorant(arb(10857)/100, G, c); lg, K1 = C.kval(E1)
    lo_, hi_ = 108.57, 200.0
    for _ in range(40):
        mid = (lo_ + hi_)/2
        if C.kval(C.majorant(arb(str(round(mid, 6))), G, c))[1] >= K1 + 1: hi_ = mid
        else: lo_ = mid
    Lk = round(hi_ + 0.005, 2)
    Ek = C.majorant(arb(str(Lk)), G, c)
    print('  %s: K = %d at L = 108.57; K = %d needs L >= %s (E* <= %s there), i.e. a Robin range to about 10^(10^%.2f)' % (
        name, K1, C.kval(Ek)[1], Lk, Ek.str(5, radius=False), (Lk - 0.834032)/2.302585), flush=True)
