# cell_detail.py -- decomposition of the d4 bound and of the density bound for a few cells (for NOTE.md)
import sys, os, json
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from flint import arb
import d4_mine as D4
P = json.load(open(os.path.join(HERE, 'out', 'cell_parts.json')))
F = json.load(open(os.path.join(HERE, 'out', 'cells_final.json')))
for key in ['80,31', '100,31', '150,31', '100,30', '200,32', '120,29']:
    p = P.get(key)
    if p is None: print(key, 'no cell'); continue
    u, t, x = [arb(float(p[n])) for n in ('u', 't', 'x')]; L = arb(p['l'])
    d4, det = D4.d4_bound(u, t, x, L, detail=True)
    rem = {k: float(v.upper()) for k, v in det.items() if k.startswith('R')}
    print(key, 'u t w x', [round(float(p[n]), 4) for n in ('u', 't', 'w', 'x')], 'nt', p['nt'], 'kp', p['kp'])
    print('   d4 <=', d4.str(8), ' D4 main', det['D4m'].str(8), ' REM', det['rem'].str(4), ' S#', det['St'].str(4))
    print('   R_i', {k: '%.2e' % v for k, v in sorted(rem.items())})
    print('   C', arb(p['C']).str(6), ' kap2', arb(p['kap2']).str(6), ' d5', arb(p['d5']).str(6), ' c1', arb(p['c1']).str(8), ' c5/d4', arb(p['c5r']).str(6), ' logN <=', F[key]['logN_up'])

