# abracket.py: lower bounds kappa(c) >= eps by the plain even certificate (cert_simple.py even sector) at chosen eps, for the
# A401172 decade decisions. Usage: abracket.py out.jsonl c:eps:bits:K2 [...]
import sys, os, json
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from run_windows import run_cert
out = sys.argv[1]
for spec in sys.argv[2:]:
    c, eps, bits, K2 = spec.split(':')
    r = run_cert([c, 'even', eps, bits, K2], dict(os.environ, OMP_NUM_THREADS='1'), 1500)
    rec = dict(c=c, eps=eps, bits=int(bits), K2=int(K2), rigorous=r.get('rigorous'), margin=r.get('margin'), order=r.get('order'), wall=r.get('wall'), rss_mb=r.get('rss_mb'))
    with open(out, 'a') as f: f.write(json.dumps(rec) + '\n')
    print(json.dumps(rec), flush=True)
