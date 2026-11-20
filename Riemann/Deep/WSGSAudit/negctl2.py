# negctl2.py (wsgsaudit): boundary sweep just outside the first chain step 2 -> 9/4. Odd certificate at eps = 2.670380e-03 (= 2.00001 U(2))
# for c in (9/4, 3]; consistency demanded: a pass at c implies eps <= lambda_1^odd(E_60, c) (a truncated eigenvalue is an upper bound).
import sys, os, json
sys.path.insert(0, '/tmp/wsgs_copy')
from flint import acb_mat, ctx, fmpq
from weil_matrix import WeilMat
from run_windows import run_cert
env = dict(os.environ, WZ_RSCAN='/tmp/ws_rscan', WZ_SCAN='/tmp/wsgs_scan', WS_M='1', WS_NG='30', OMP_NUM_THREADS='1')
eps = '2.670380e-03'
for cs in ('9/4', '7/3', '19/8', '29/12', '5/2', '21/8', '11/4', '3'):
    p, q = (cs.split('/') + ['1'])[:2]; ctx.prec = 192
    T, _ = WeilMat(fmpq(int(p), int(q))).block(60, 'odd')
    o1 = sorted([e.real.mid() for e in acb_mat(T).eig(algorithm='approx')], key=lambda z: z.mid())[0]
    r = run_cert([cs, 'odd', eps, '192', '401'], env, 600)
    o1f = float(o1.str(20, radius=False)); ok = (not r.get('rigorous')) or float(eps) <= o1f
    print(json.dumps(dict(c=cs, lambda1_odd_E60=o1.str(6, radius=False), eps=eps, rigorous=r.get('rigorous'), margin=r.get('margin'), consistent=ok)), flush=True)
