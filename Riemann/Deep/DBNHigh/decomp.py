#!/usr/bin/env python3
# dbnhigh decomp.py: endpoint decomposition of the theta141 chain into its source pieces (float, scheme 1).
# Each variant removes (pays back) one debit in every row, or the own pair, and reruns the whole chain.
import json, math, sys
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from fractions import Fraction as F
rows = json.load(open(sys.argv[1]))
import importlib.util
spec = importlib.util.spec_from_file_location('m', __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), 'model.py'))
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
R = m.load_rows(sys.argv[1])
direct = {r['row']: r for r in json.load(open(sys.argv[2]))['rows']}
for i, r in enumerate(R):
    hc = direct[i]['head_construction']
    if 'raw_head_bounds' in hc:
        r['Hm'], r['Hp'] = [float(F(x)) for x in hc['raw_head_bounds']]
        r['cred'] = float(F(hc['composite_credit'])) + float(F(hc['old_angular_credit'])) + float(F(hc['additional_clock_credit']))
    else:
        r['Hm'], r['Hp'] = [float(F(x)) for x in hc['raw_cap_sums']]
        r['cred'] = float(F(hc['fixed_credit'])) + float(F(hc['coordinate_credit'])) + float(F(hc['orbit_credit']))
    r['a'] = r['S0']; r['b'] = float(F(rows[i]['source_lower_slope']))
def own(q, row): return float(m.R_row(row, __import__('numpy').zeros(1), math.sqrt(q))[0])
B0 = 0.8733458561067187
def run(dS=lambda i, th: 0.0, Ofac=1.0, stop=1e-4):
    dt = 1/2048000; q = (2*B0 - 1)**2 + 1e-12; k = 0; i = 0
    while True:
        t = k*dt
        while t + dt > R[i]['t1'] + 1e-15: i += 1
        row = R[i]; th = (t - row['t0'])/(row['t1'] - row['t0'])
        S = row['a'] + row['b']*(t - row['t0']) + dS(i, th)
        r = 2*(S - Ofac*ownf(q, row)) - 1e-12
        q = (q + 2/r)*math.exp(-r*dt) - 2/r; k += 1
        if q <= stop: return k*dt + q/2
def ownf(q, row):
    s = 0.0
    for h, a, w in row['nodes']:
        A = h*h + a*a - q; s += 2*w*h*A/(A*A + 4*a*a*q)
    if row['ov'] is None: return s
    al, cs = row['ov']; s *= al; y = math.sqrt(q)
    for p, bs in cs:
        for mm, b in enumerate(bs): s += 2*b*p**(mm-1)*(-1)**mm*((p - y)**-(mm+1) + (p + y)**-(mm+1))
    return s
E0 = run(); out = dict(base=E0)
print('base', repr(E0), flush=True)
V = dict(
  native_to_0=lambda i, th: R[i]['native'],
  credits_removed=lambda i, th: -R[i]['cred'],
  raw_head_to_0=lambda i, th: (1 - th)*R[i]['Hm'] + th*R[i]['Hp'] - R[i]['cred'],
  head_minus_10pct=lambda i, th: 0.1*((1 - th)*R[i]['Hm'] + th*R[i]['Hp']),
  main_plus_1pct=lambda i, th: 0.01*R[i]['main'],
  native_x2=lambda i, th: -R[i]['native'])
for k, f in V.items():
    E = run(dS=f); out[k] = E; print(k, repr(E), 'dE %.3e' % (E - E0), flush=True)
for fac in (0.99, 0.9):
    E = run(Ofac=fac); out['own_pair_x%.2f' % fac] = E; print('own_pair_x%.2f' % fac, repr(E), 'dE %.3e' % (E - E0), flush=True)
json.dump(out, open(sys.argv[3], 'w'), indent=1)
