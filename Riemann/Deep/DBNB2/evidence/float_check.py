#!/usr/bin/env python3
"""dbnb2 independent float re-integration (written from the TeX formulas, not the producer).

Rows: verification/arcs/ACTUAL_67_ROW_ARC_INPUTS.json (intercept, slope, signed nodes, jets).
Own pair: 2wh(h^2+v^2-q)/((h^2+v^2-q)^2+4v^2 q); jets (-1)^m[(p-y)^-(m+1)+(p+y)^-(m+1)], y=sqrt q.
Scheme 1 (certificate model): rate frozen per arc at (t0,q0), exact exponential, no rounding.
Scheme 2 (ideal): RK4 of Q'=-2-(2(S(t)-O(Q))-1e-12)Q, same step.
"""
import json, math, sys
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
            Sm1 = (p - y)**-(m+1) + (p + y)**-(m+1)
            s += 2*b*p**(m-1)*(-1)**m*Sm1
    return s

def src(t, row): return row['a'] + row['b']*(t - row['start'])

def run(b0, scheme, stop=1e-4):
    dt = 1/2048000; q = (2*b0 - 1)**2 + 1e-12; k = 0; i = 0
    while True:
        t = k*dt
        while t + dt > R[i]['end'] + 1e-15: i += 1
        row = R[i]
        assert q <= row['cap']
        if scheme == 1:
            r = 2*(src(t, row) - own(q, row)) - 1e-12
            q = (q + 2/r)*math.exp(-r*dt) - 2/r
        else:
            f = lambda tt, qq: -2 - (2*(src(tt, row) - own(qq, row)) - 1e-12)*qq
            k1 = f(t, q); k2 = f(t+dt/2, q+dt/2*k1); k3 = f(t+dt/2, q+dt/2*k2); k4 = f(t+dt, q+dt*k3)
            q = q + dt/6*(k1 + 2*k2 + 2*k3 + k4)
        k += 1
        if q <= stop: return k, k*dt, q, k*dt + q/2

for name, b0 in (('7/8', 7/8), ('b0', 0.87338196661320752)):
    for sc in (1, 2):
        k, T, q, e = run(b0, sc)
        print(name, 'scheme', sc, 'arcs', k, 'T', T, 'qT', q, 'endpoint', repr(e))
# local slope dLambda/db0 (certificate model)
e1 = run(0.87338196661320752, 1)[3]; e2 = run(0.87338196661320752 + 1e-4, 1)[3]
print('dLambda/db0 near b0 (model):', (e2 - e1)/1e-4)
# author's research number at 0.8748867 (TeX l.210: 0.107936602217)
print('model at 0.8748867:', repr(run(0.8748867, 1)[3]))
# stop-threshold option
print('stop 2e-6 (model):', repr(run(0.87338196661320752, 1, stop=2e-6)[3]))
