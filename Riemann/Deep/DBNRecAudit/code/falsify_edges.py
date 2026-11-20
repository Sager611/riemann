#!/usr/bin/env python3
# dbnrecaudit falsify_edges.py: (L2) on the actual signed head at the domain edges (falsification only): the band edge x = X (u = 0), the
# shell-0/1 boundary u = 1/100, and u = u_F (start of the far chord), at t = t-, t+; 10^6 points on x' in [c, c + 2] per window (C helper 'dense').
# usage: falsify_edges.py ROW [ROW ...]   (cwd = repo root; writes Riemann/Deep/DBNRecAudit/logs/edges_r<ROW>.json)
import sys, json, subprocess, os
sys.argv, args = sys.argv[:1], sys.argv[1:]
sys.path.insert(0, 'Riemann/Deep/DBNRecAudit/code')
import falsify_dense as fd
mp = fd.mp
from fractions import Fraction as Fr
K_ = fd.K; ln = fd.ln; X = fd.X; PI = fd.PI
q = lambda r: mp.mpf(r.numerator)/r.denominator
for row in map(int, args):
    p = fd.plan[row]; ar = fd.aff[row]
    tm, tp = (q(Fr(v)) for v in p['time'])
    nodes = [(Fr(n['height']), Fr(n['shift']), Fr(n['weight'])) for n in p['nodes']]
    M = sum(w for _, _, w in nodes); s = (1 - Fr(ar['beta']))*M/4; sf = q(s)
    cdir = 'Riemann/Deep/DBNRec/logs/' if row <= 3 else 'Riemann/Deep/DBNRec2/logs/'
    cert = json.load(open(cdir + 'cert_r%d.json' % row)); Ap, Bp = q(Fr(cert['A_new'])), q(Fr(cert['B_new'])); uF = q(Fr(cert['far']['uF']))
    F = [mp.mpf(0)]*(K_ + 1)
    for n in range(2, K_ + 1): F[n] = sum(q(w)*mp.exp(-q(h)*ln[n]/2)*mp.cos(q(a)*ln[n]/2) for h, a, w in nodes)
    ds = {0: fd.dcoef(tm), 1: fd.dcoef(tp)}
    res = dict(row=row, uF=float(uF), runs=[])
    for name, u in (('band edge u=0', mp.mpf(0)), ('shell boundary u=1/100', mp.mpf(1)/100), ('far-chord start u=uF', uF)):
        for th, t in ((0, tm), (1, tp)):
            x0 = X*mp.exp(u); c = x0 + PI*t/4          # x' window [c, c + 2] <-> x in [x0, x0 + 2]
            L = mp.log(x0/(4*PI)); dt = ds[th]
            f = fd.TMP + 'e_r%d.txt' % row
            with open(f, 'w') as fo:
                fo.write('%d\n' % (K_ - 1))
                for n in range(2, K_ + 1):
                    amp = dt[n]*mp.exp(-(mp.mpf(1)/2 + t*L/4)*ln[n])*F[n]/2
                    fo.write('%.17g %.17g 0 %.17g\n' % (float(ln[n]/2), float(amp), float(mp.fmod(c*ln[n]/2, 2*PI))))
            out = subprocess.run([fd.DENSE, f, 'a', '0', '2', '1000000'], capture_output=True, text=True, check=True).stdout.split()
            os.remove(f)
            G = float(out[0]) - float(sf*u); rhs = float((1 - th)*Ap + th*Bp)
            res['runs'].append(dict(window=name, theta=th, Qmax=float(out[0]), G=G, rhs=rhs, margin=rhs - G))
    json.dump(res, open('Riemann/Deep/DBNRecAudit/logs/edges_r%d.json' % row, 'w'), indent=1)
    print(json.dumps(dict(row=row, minmargin=min(r['margin'] for r in res['runs']), runs=[(r['window'], r['theta'], round(r['Qmax'], 6), round(r['margin'], 4)) for r in res['runs']])), flush=True)
