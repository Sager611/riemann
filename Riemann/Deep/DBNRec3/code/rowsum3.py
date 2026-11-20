#!/usr/bin/env python3
# dbnrec3 rowsum3.py (log summary only, no computation): one markdown bullet per certified row from logs/cert_r<row>.json, prep.json and progress_r<row>.txt.
import sys, json, re, os
L = 'Riemann/Deep/DBNRec3/logs/'
def fq(v): n, d = (v.split('/') + ['1'])[:2]; return int(n)/int(d)
for r in map(int, sys.argv[1:]):
    c = json.load(open(L + 'cert_r%d.json' % r)); sh = c['shells']; s0 = sh[0]
    pj = json.loads(open(L + 'r%d_00.prep.json' % r).readline())
    T = fq(json.loads(open(L + 'r%d_01.prep.json' % r).readline())['T']) if os.path.exists(L + 'r%d_01.prep.json' % r) else float('nan')
    cpu = [int(m) for m in re.findall(r'maxUB3 (\d+) s', open(L + 'progress_r%d.txt' % r).read())]
    w = s0.get('worst') or {}; cell = w.get('cell', '').split()
    hot = [s['tag'] for s in sh[1:] if s['cells']]
    print('- Row %d (t ∈ [%s, %s], s = %.8f): shell 0 at T₀ = %.6f: %d unresolved cells; binding cell anchor %s (u = %.6f), Arb UB %.8f (float %.9f), S₀ = %.9f; T = %.6f; shells 1–%d (to u_F = %.4f): %s. bb CPU %d s (%s); certified (logs/cert_r%d.json): A′ = %.9f, B′ = %.9f; far a_F = %.6f, b_F = %.6f. Raymond\'s A = %.9f, B = %.9f: **dA = %.6f, dB = %.6f**.'
          % (r, pj['t'][0], pj['t'][1], pj['s_float'], s0['thr'], s0['cells'], cell[0] if cell else '-', w.get('u', float('nan')), w.get('ub_arb', float('nan')), w.get('ub_float', float('nan')), s0['S'],
             T, len(sh) - 1, fq(c['far']['uF']), 'cells above threshold in ' + ', '.join(hot) if hot else 'no cell above threshold', sum(cpu), '+'.join(map(str, cpu)), r,
             c['A_new_float'], c['B_new_float'], c['far']['a'], c['far']['b'], pj['A'], pj['B'], c['dA'], c['dB']))
