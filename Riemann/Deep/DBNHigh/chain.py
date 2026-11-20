#!/usr/bin/env python3
# dbnhigh chain.py: float model of Raymond's 67-row chain restarted at theta141 (scheme 1 of HP141/DBNB2/float_check.py).
# Outputs per-row visited (t,q) ranges, S, O, rate; endpoint; per-row sensitivity dE/dS_i (finite difference).
import json, math, sys, time
from fractions import Fraction as F
rows = json.load(open(sys.argv[1]))
R = []
for r in rows:
    nodes = [(float(F(n['height'])), float(F(n['shift'])), float(F(n['weight']))) for n in r['nodes']]
    o = r['joint_jet_overlay']
    ov = None if o is None else (float(F(o['alpha'])), [(float(F(c['p'])), [float(F(b)) for b in c['b']]) for c in o['centres']])
    R.append(dict(start=float(F(r['time'][0])), end=float(F(r['time'][1])), cap=float(F(r['q_interval'][1])),
                  a=float(F(r['source_lower_intercept_at_row_start'])), b=float(F(r['source_lower_slope'])), nodes=nodes, ov=ov))
def own(q, row):
    s = 0.0
    for h, v, w in row['nodes']:
        A = h*h + v*v - q
        s += 2*w*h*A/(A*A + 4*v*v*q)
    if row['ov'] is None: return s
    alpha, cs = row['ov']; s *= alpha; y = math.sqrt(q)
    for p, bs in cs:
        for m, b in enumerate(bs):
            s += 2*b*p**(m-1)*(-1)**m*((p - y)**-(m+1) + (p + y)**-(m+1))
    return s
B0 = 0.8733458561067187
def run(dS=None, stop=1e-4, log=None):
    dt = 1/2048000; q = (2*B0 - 1)**2 + 1e-12; k = 0; i = 0
    while True:
        t = k*dt
        while t + dt > R[i]['end'] + 1e-15: i += 1
        row = R[i]
        S = row['a'] + row['b']*(t - row['start']) + (dS[i] if dS else 0.0)
        O = own(q, row)
        r = 2*(S - O) - 1e-12
        if log is not None:
            g = log.setdefault(i, dict(t0=t, q0=q, S0=S, O0=O, r0=r))
            g.update(t1=t + dt, q1=None, S1=S, O1=O, r1=r)
        qn = (q + 2/r)*math.exp(-r*dt) - 2/r
        if log is not None: log[i]['q1'] = qn
        q = qn; k += 1
        if q <= stop: return k*dt + q/2, k, i
t0 = time.time()
log = {}
E, k, ilast = run(log=log)
print('endpoint', repr(E), 'arcs', k, 'last row', ilast, 's', round(time.time() - t0, 1), flush=True)
out = dict(endpoint=E, arcs=k, rows={})
for i, g in sorted(log.items()):
    out['rows'][i] = g
# sensitivities: dE/dS_i by finite difference with dS = 1e-3 on row i
h = 1e-3
sens = {}
for i in sorted(log):
    dS = [0.0]*67; dS[i] = h
    Ei = run(dS=dS)[0]
    sens[i] = (Ei - E)/h
print('sum dE/dS_i', sum(sens.values()), flush=True)
out['sens'] = sens
json.dump(out, open(sys.argv[2], 'w'), indent=0)
for i in sorted(log):
    g = log[i]
    print('%2d t %.5f-%.5f q %.5f-%.5f S %.4f O %.4f-%.4f r %.3f-%.3f dE/dS %.3e cap %.4f' % (i, g['t0'], g['t1'], g['q0'], g['q1'], g['S0'], g['O0'], g['O1'], g['r0'], g['r1'], sens[i], R[i]['cap']))
