#!/usr/bin/env python3
# dbnrecaudit check_inputs.py: rows 1-8: plan (M8_LATE_PLAN) = arc inputs (prep_end's source) in time, nodes, mode, jet; source t_lo/t_hi = plan;
# plan node mass = AFFINE mass; native_debit (final) <= AFFINE N0 (refund >= 0); main_lo <= M c0 (Arb); X - 4 - 4 pi N0^2 > 0 (Arb).
import json
from fractions import Fraction as Fr
from flint import arb, fmpq, ctx
ctx.prec = 200
DL = '.scratch/dbnrecaudit/dl/Proof_and_Replay__Lambda_Peer_Review__'
plan = {int(r['row']): r for r in json.load(open(DL + 'inputs__workspace__lambda_round14__material__m8_late__M8_LATE_PLAN.json'))['rows']}
arcs = {int(r['row']): r for r in json.load(open(DL + 'verification__arcs__ACTUAL_67_ROW_ARC_INPUTS.json'))}
A = json.load(open(DL + 'verification__referee__ALL_67_ORIGINAL_AFFINE_NATIVE.json')); aff = {int(r['row']): r for r in A['rows']}
src = json.load(open(DL + 'inputs__workspace__lambda_round15__flow_main_v2__MAIN_CLOCK_SOURCE.json'))['records']
print('AFFINE original_native_nodes', A['original_native_nodes'], '| input_sha256', json.dumps(A['input_sha256'])[:700])
print('AFFINE claim', str(A['mathematical_claim'])[:500])
def nodes(lst): return [(Fr(n['height']), Fr(n['shift']), Fr(n['weight'])) for n in lst]
N0 = 690988; c0 = arb(N0).log()/2
print('X-4-4piN0^2 =', (arb(6000000185827 - 4) - 4*arb.pi()*N0**2).str(10))
for i in range(1, 9):
    p, a, s, f = plan[i], arcs[i], src[i], aff[i]
    ok = [p['time'] == a['time'], nodes(p['nodes']) == nodes(a['nodes']), p['mode'] == a['mode'] == s['mode'] == 'retained',
          p.get('joint_jet_overlay') is None and a['joint_jet_overlay'] is None, [Fr(s['t_lo']), Fr(s['t_hi'])] == [Fr(x) for x in p['time']]]
    M = sum(w for _, _, w in nodes(p['nodes']))
    mlo = Fr(s['main_lo']); Mc0 = arb(fmpq(M.numerator, M.denominator))*c0
    nd, n0 = Fr(s['native_debit']), Fr(f['N0'])
    print(i, all(ok), ok, 'mass=AFF', M == Fr(f['mass']), 'Mc0-main_lo', (Mc0 - arb(fmpq(mlo.numerator, mlo.denominator))).str(3),
          'native_debit %.9f AFF_N0 %.9f refund %.3e beta %s nodes %d' % (float(nd), float(n0), float(n0 - nd), f['beta'], len(p['nodes'])))
