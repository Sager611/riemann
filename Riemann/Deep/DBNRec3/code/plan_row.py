#!/usr/bin/env python3
# dbnrec plan_row.py (float helper for thresholds and shell lists only; nothing certified here).
#  plan_row.py t0 ROW PREP.json        -> shell-0 target T0 = tri+ - 0.965 cr*(kappa 1.32)  (row 1: 0.7715 gave S = 0.7740)
#  plan_row.py T  ROW PREP.json UNRES   -> T = max(thr, max over dumped cells of ub - s log(x/X))  (location-aware, float)
#  plan_row.py shells PREP.json T       -> shell boundaries 0.01 ... uF, uF >= max((tri+ - T)/s, (tri- - T - D)/s) + 0.002
import sys, json, math
X = 6000000185827.0; D = 4*math.pi/math.log(2)
cmd = sys.argv[1]
if cmd == 't0':
    row = int(sys.argv[2]); pj = json.loads(open(sys.argv[3]).readline())
    cr = [r for r in json.load(open('Riemann/Deep/CeilAudit5/logs/selfcons.json')) if r['row'] == row][0]['cr_1.32']
    print('%.6f' % (pj['tri_plus'] - 0.965*cr))
elif cmd == 'T':
    pj = json.loads(open(sys.argv[3]).readline()); s = pj['s_float']; best = pj['thr']
    for line in open(sys.argv[4]):
        w = line.split()
        off = float(w[0])*D + float(w[1]) if len(w) == 4 else float(w[0]); ub = float(w[-1])
        best = max(best, ub - s*math.log1p((off - 1.0)/X))
    print('%.6f' % (best + 2e-6))
else:
    pj = json.loads(open(sys.argv[2]).readline()); T = float(sys.argv[3]); s = pj['s_float']
    uF = max((pj['tri_plus'] - T)/s, (pj['tri_minus'] - T - pj['D_up'])/s) + 0.002
    us = [0.01, 0.02, 0.035, 0.05]; u = 0.05
    while u + 0.02 < uF: u = round(u + 0.02, 6); us.append(u)
    us.append(round(uF, 4)); print(' '.join('%g' % v for v in us))
