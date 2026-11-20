# negctl.py (wsgsaudit): negative controls / falsification of the WeilSimple certificates (run from the /tmp copy, one process at a time).
# A sound certificate must FAIL whenever eps exceeds an upper bound of the eigenvalue it claims to bound from below:
# truncated eigenvalues are upper bounds (min-max on E_N), lambda_k(N) >= lambda_k.
import sys, os, json, math
sys.path.insert(0, '/tmp/wsgs_copy')
from flint import arb, acb_mat, ctx, fmpq
from weil_matrix import WeilMat
from run_windows import run_cert
env = dict(os.environ, WZ_RSCAN='/tmp/ws_rscan', WZ_SCAN='/tmp/wsgs_scan', WS_M='1', OMP_NUM_THREADS='1')
def eigs(cs, N, prec):
    ctx.prec = prec; W = WeilMat(fmpq(*[int(v) for v in cs.split('/')]) if '/' in cs else fmpq(int(cs))); out = {}
    for sec in ('even', 'odd'):
        T, _ = W.block(N, sec)
        out[sec] = sorted([e.real.mid() for e in acb_mat(T).eig(algorithm='approx')], key=lambda z: z.mid())[:2]
    return out
def run(tag, cs, kind, eps, bits, K2, NG, expect):
    e = dict(env, WS_NG=str(NG)); r = run_cert([cs, kind, eps, str(bits), str(K2)], e, 1200)
    got = r.get('rigorous'); line = dict(tag=tag, c=cs, kind=kind, eps=eps, K2=K2, NG=NG, rigorous=got, margin=r.get('margin'), U_up=r.get('U_up'),
                                       rss=r.get('rss_mb'), wall=r.get('wall'), expect=expect, PASS=(got == expect))
    print(json.dumps(line), flush=True)
for cs, N, bits, NG in (('5', 60, 256, 30), ('31/4', 80, 320, 60)):
    ev = eigs(cs, N, bits)
    o1 = ev['odd'][0]; e2 = ev['even'][1]
    print(json.dumps(dict(c=cs, N=N, e1=ev['even'][0].str(8, radius=False), e2=e2.str(8, radius=False), o1=o1.str(8, radius=False))), flush=True)
    run('odd eps=1.02*o1(N)', cs, 'odd', '%.6e' % (1.02*float(o1.str(20, radius=False))), bits, 401, NG, False)
    run('defl eps=1.02*e2(N)', cs, 'defl', '%.6e' % (1.02*float(e2.str(20, radius=False))), bits, 401, NG, False)
    run('odd eps=0.5*o1(N)', cs, 'odd', '%.6e' % (0.5*float(o1.str(20, radius=False))), bits, 401, NG, None)
run('stretched step 2->3 (eps = 2.00001 U(2))', '3', 'odd', '2.670380e-03', 192, 401, 30, False)
run('K2=201 at c=7 odd (logged failure)', '7', 'odd', '4.779037e-26', 256, 201, 40, False)
run('K2=401 at c=7 odd (logged success)', '7', 'odd', '4.779037e-26', 256, 401, 40, True)
