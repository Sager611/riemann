#!/usr/bin/env python3
# dbnrec check.py: exact/Arb checks of the lemma's numeric premises (I6, I7) for the given rows.
# (I6) X_* = X - 4 > 4 pi N0^2 (so Delta = c - c0 >= log(x/X)/4 for x >= X); (I7) main_lo <= M c0 with M = sum_j w_j = AFFINE mass, s = retained_main_mass/4.
# usage: check.py row ...
import sys, json
from flint import arb, fmpq, ctx
ctx.prec = 200
def Q(v): return fmpq(*map(int, v.split('/'))) if '/' in v else fmpq(int(v))
X = arb(6000000185827); N0 = 690988; c0 = arb(N0).log()/2
r6 = (arb(X - 4) - 4*arb.pi()*N0**2)
print('I6: X-4 - 4 pi N0^2 =', r6.str(10), 'OK' if r6 > 0 else 'FAIL')
src = json.load(open('.scratch/dbnrec2/ws/lambda_round15/flow_main_v2/MAIN_CLOCK_SOURCE.json'))['records']
plan = json.load(open('.scratch/dbnrec2/ws/lambda_round14/material/m8_late/M8_LATE_PLAN.json'))['rows']
aff = {int(r['row']): r for r in json.load(open('.scratch/dbnrec2/dl/Proof_and_Replay__Lambda_Peer_Review__verification__referee__ALL_67_ORIGINAL_AFFINE_NATIVE.json'))['rows']}
for i in map(int, sys.argv[1:]):
    s, p, a = src[i], plan[i], aff[i]
    M = sum((Q(n['weight']) for n in p['nodes']), fmpq(0))
    ok_mass = M == Q(a['mass']); beta = Q(a['beta']); ok_ret = (1 - beta)*M == Q(a['retained_main_mass'])
    gap = arb(M)*c0 - arb(Q(s['main_lo']))
    print('row %d: mode %s jet %s | M = sum w = AFFINE mass: %s | (1-beta)M = retained_main_mass: %s | beta %s | M c0 - main_lo = %s %s'
          % (i, s['mode'], p.get('joint_jet_overlay') is not None, ok_mass, ok_ret, beta, gap.str(8), 'OK' if gap > 0 else 'FAIL'))
