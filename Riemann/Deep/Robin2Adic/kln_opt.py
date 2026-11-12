#!/usr/bin/env python3
# kln_opt.py -- exploratory: re-optimise KLN Lemma 4.14 parameters (k, mu, alpha, delta, d) for W(s) near s = 0.85-0.873 with RH(3e12),
# using Ingham4's KP class (its admissibility test = KLN Lemma 4.14 hypotheses as coded there). Prints the best rows; not a certificate.
#   uv run --no-project --with python-flint python3 Riemann/Deep/Robin2Adic/kln_opt.py
import os, sys, itertools, importlib.util, time
sys.dont_write_bytecode = True                             # no __pycache__ next to the imported files
HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location('r2a_cert', os.path.join(HERE, 'cert.py'))
C = importlib.util.module_from_spec(spec); spec.loader.exec_module(C)
from flint import arb
I4 = C.I4
def make_rows(rows):
    old = (I4.H0, I4.HH, I4.LH0, I4.llH0)
    I4.H0 = arb('3e12'); I4.HH = I4.H0 - 1; I4.LH0 = I4.H0.log(); I4.llH0 = I4.LH0.log()
    try: P = [I4.KP(r) for r in rows]
    finally: I4.H0, I4.HH, I4.LH0, I4.llH0 = old
    return P
def W_row(p, s):
    if not (s > p.smin) or not p.ok: return None
    a = 8*(1 - s)/3; q = 2 - a; e = 5 - 2*s
    return 4*(p.logA(s).exp()*C.Ig(q, e) + p.Bc(s).max(arb(0))*C.I2)
S = [arb(int(v*2**16))/2**16 for v in (0.85, 0.86, 0.87, 0.873)]
base = make_rows([r for r in I4.KLN_ROWS if r[0] in ('0.85', '0.86', '0.87')])
bw = [min(float(W_row(p, s).mid()) for p in base if W_row(p, s) is not None) for s in S]
print('table rows: W(s) at s = 0.85, 0.86, 0.87, 0.873:', ['%.4e' % v for v in bw], flush=True)
t0 = time.time(); best = [(None, 1e9)]*len(S); n_ok = 0; n = 0
mus = ['1.245', '1.2', '1.15', '1.1']; als = ['%.3f' % (0.06 + 0.02*i) for i in range(10)]; des = ['0.15', '0.2', '0.25', '0.303', '0.4', '0.5', '0.7']
ds = ['%.3f' % (0.20 + 0.04*i) for i in range(9)]
for mu, al, de, d in itertools.product(mus, als, des, ds):
    row = ('0.86', '1.0', mu, al, de, d, '0', '0'); n += 1
    try: p = make_rows([row])[0]
    except Exception as ex: continue
    if not p.ok: continue
    n_ok += 1
    for i, s in enumerate(S):
        w = W_row(p, s)
        if w is not None and float(w.mid()) < best[i][1]: best[i] = (row, float(w.mid()))
    if time.time() - t0 > 900: print('time cap'); break
print('tried %d, admissible %d, %.0f s' % (n, n_ok, time.time() - t0), flush=True)
for s, (row, w), b in zip(S, best, bw): print('s=%.4f best W %.4e (table %.4e, ratio %.3f) row %s' % (float(s.mid()), w, b, w/b, row), flush=True)
