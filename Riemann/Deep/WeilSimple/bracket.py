# bracket.py: two-sided bounds L_1 <= mu_lambda <= U for the ground eigenvalue at the integer windows: the plain even certificate
# of cert_simple.py (= cert_zeta.py even sector) at eps = 0.5 U, 0.2 U, 0.05 U (first that succeeds). Usage: bracket.py logs/int.jsonl logs/bracket.jsonl
import sys, os, json, math
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from run_windows import run_cert
rows = [json.loads(l) for l in open(sys.argv[1])]
for r in sorted(rows, key=lambda r: int(r['c'])):
    c = r['c']; cf = float(c); U = float(r['U_up']); K2 = max(401, math.ceil(5*cf*cf/3)); K2 += (K2 % 2 == 0)
    env = dict(os.environ, OMP_NUM_THREADS='1')
    out = None
    for f in (0.5, 0.2, 0.05):
        eps = '%.4e' % (f*U)
        res = run_cert([c, 'even', eps, str(r['bits']), str(K2)], env, 1800)
        if res.get('rigorous'): out = dict(c=c, L1=eps, U=r['U_up'], margin=res.get('margin'), wall=res.get('wall'), rss_mb=res.get('rss_mb')); break
    if out is None: out = dict(c=c, L1=None, U=r['U_up'])
    with open(sys.argv[2], 'a') as fh: fh.write(json.dumps(out) + '\n')
    print(json.dumps(out), flush=True)
