#!/usr/bin/env python3
# dbnhigh model.py: float model of one source row of Raymond's chain (Zenodo 23752103), from the TeX formulas.
# S(v) = c0 v(0) - (1/2) sum_{n=2}^{N} d_t(n) n^{-1/2-t c0} |v(log n)|  (head n<=4096, native n>4096, no credits, no RS errors)
# O(q) = R(0,sqrt q);  R(d,b) = sum_j w_j (1/2)[K_h(d+a,b)+K_h(d-a,b)],  K_h(d,b) = Re[1/(h-b+id) + 1/(h+b+id)]
import json, math, sys
import numpy as np
from fractions import Fraction as F
X = 6000000185827.0
L0 = math.log(X/(4*math.pi)); c0 = L0/4
NMAX = 690988
def dt_coeffs(t, N=NMAX):
    n = np.arange(N + 1, dtype=float); ln = np.zeros(N + 1); ln[1:] = np.log(n[1:])
    b = np.exp(t*ln*ln/4)
    d = ln*b  # (log n) b_t(n)
    d[0] = 0; d[1] = 0
    for r in range(2, N//2 + 1):
        if d[r] != 0.0:
            m = np.arange(2, N//r + 1)
            d[r*m] -= d[r]*b[m]
    return d, ln
def load_rows(path):
    rows = json.load(open(path)); out = []
    for r in rows:
        nodes = [(float(F(x['height'])), float(F(x['shift'])), float(F(x['weight']))) for x in r['nodes']]
        o = r['joint_jet_overlay']
        ov = None if o is None else (float(F(o['alpha'])), [(float(F(c['p'])), [float(F(bb)) for bb in c['b']]) for c in o['centres']])
        out.append(dict(t0=float(F(r['time'][0])), t1=float(F(r['time'][1])), cap=float(F(r['q_interval'][1])), nodes=nodes, ov=ov,
                        main=float(F(r['main_lo'])), native=float(F(r['native_debit'])), A=float(F(r['first_hi'])), B=float(F(r['last_hi'])),
                        S0=float(F(r['source_lower_intercept_at_row_start']))))
    return out
def v_nodes(nodes, u):
    v = np.zeros_like(u)
    for h, a, w in nodes: v += w*np.exp(-h*u/2)*np.cos(a*u/2)
    return v
def v_jet(ov, u):
    # joint jet centre p, coefficients b_m: zero-side 2 b_m p^(m-1) (-1)^m [(p-y)^-(m+1)+(p+y)^-(m+1)] at d=0
    # = b_m * (2 p^(m-1)) * (1/m!)(d/dp)^m [1/(p-y)+1/(p+y)] -> Dirichlet profile b_m 2 p^(m-1) (-u/2)^m e^{-pu/2}/m!
    v = np.zeros_like(u)
    for p, bs in ov[1]:
        for m, bm in enumerate(bs):
            v += bm*2*p**(m - 1)*(-u/2)**m*np.exp(-p*u/2)/math.factorial(m)
    return v
def v_row(row, u):
    v = v_nodes(row['nodes'], u)
    if row['ov'] is not None: v = row['ov'][0]*v + v_jet(row['ov'], u)
    return v
def R_nodes(nodes, d, b):
    s = 0.0
    for h, a, w in nodes:
        for dd in (d + a, d - a):
            s = s + w*0.5*((h - b)/(dd*dd + (h - b)**2) + (h + b)/(dd*dd + (h + b)**2))
    return s
def R_jet(ov, d, b):
    s = 0.0
    for p, bs in ov[1]:
        for m, bm in enumerate(bs):
            z1 = (p - b) + 1j*d; z2 = (p + b) + 1j*d
            s = s + bm*2*p**(m - 1)*(-1)**m*np.real(z1**(-m - 1) + z2**(-m - 1))
    return s
def R_row(row, d, b):
    r = R_nodes(row['nodes'], d, b)
    if row['ov'] is not None: r = row['ov'][0]*r + R_jet(row['ov'], d, b)
    return r
def band_weights(t, d, ln, nlo=2, nhi=NMAX):
    n = np.arange(nlo, nhi + 1)
    return n, 0.5*d[n]*np.exp(-(0.5 + t*c0)*ln[n])
if __name__ == '__main__':
    rows = load_rows(sys.argv[1])
    for i in [int(x) for x in sys.argv[2].split(',')]:
        row = rows[i]
        for t in (row['t0'], row['t1']):
            d, ln = dt_coeffs(row['t1'])  # TeX: d_{t+} with n^{-1/2 - t c0}
            n, beta = band_weights(t, d, ln)
            v = v_row(row, ln[n])
            head = float(np.sum((beta*np.abs(v))[n <= 4096])); nat = float(np.sum((beta*np.abs(v))[n > 4096]))
            print('row %d t %.5f main %.6f (file %.6f) head %.6f (file A %.6f B %.6f) tail(4096,N] %.6f (file native %.6f)' % (i, t, c0*float(v_row(row, np.zeros(1))[0]), row['main'], head, row['A'], row['B'], nat, row['native']), flush=True)
        dd = np.concatenate([np.linspace(0, 5, 501), np.linspace(5.05, 30, 500), np.linspace(30.5, 300, 540)])
        worst = []
        for qq in np.linspace(0, row['cap'], 9):
            y = math.sqrt(qq); Rv = R_row(row, dd, y); bound = 4/(dd*dd + 4*qq) if qq > 0 else np.where(dd > 0, 4/np.maximum(dd*dd, 1e-300), np.inf)
            rel = (bound - Rv)/bound; k = int(np.argmin(rel)); worst.append((round(qq, 4), float(rel[k]), float(dd[k])))
        print('  kernel min rel slack by q (q, slack, d):', worst)
        print('  O(cap)=%.6f O(0)=%.6f far-field d^2 R(1e4,0)=%.7f' % (float(R_row(row, np.zeros(1), math.sqrt(row['cap']))[0]), float(R_row(row, np.zeros(1), 0.0)[0]), float(1e8*R_row(row, np.array([1e4]), 0.0)[0])))
