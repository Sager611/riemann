# 2-level CEGAR on a sample of the extreme type-3 body: is there S (all pairs <= 1, containing the contact p1) whose
# graph of pairs with distance in [t, 1] is not 5-colourable?  CANDIDATE = a diameter-<=1 set with beta5 >= t.
import numpy as np, sys, json, time
sys.path.insert(0, 'Borsuk/code')
from p6_pd import body, interior
from p6_survey import boundary, d2
from cegar import cegar
spec = json.loads(open(sys.argv[1]).read()); M = int(sys.argv[2]); t = float(sys.argv[3]); seed = int(sys.argv[4]); tl = int(sys.argv[5])
rng = np.random.default_rng(seed); rho, P, cuts = body(spec)
C = np.concatenate([P, boundary(rho, P, cuts, M, rng), interior(rho, P, cuts, M // 5, rng)]); m = len(C)
D = np.zeros((m, m), np.float32)
for i0 in range(0, m, 500): D[i0:i0 + 500] = np.sqrt(np.maximum(d2(C[i0:i0 + 500], C), 0))
BAD = D > 1.0
U, V = np.nonzero(np.triu((D >= t) & (D <= 1.0), 1)); U = U.astype(np.int32); V = V.astype(np.int32)
del D
t0 = time.time(); res = cegar(m, U, V, BAD, tl=tl, root=0)
out = {'m': m, 't': t, 'conflict_edges': int(len(U)), 'res': res[0], 'info': [int(x) if not isinstance(x, (list, np.ndarray, str)) else (len(x) if not isinstance(x, str) else x) for x in res[1:]], 'sec': round(time.time() - t0)}
print(json.dumps(out), flush=True)
if res[0] == 'CANDIDATE':
    S = np.asarray(res[1]); np.save('Borsuk/.scratch/fuzzy_S_%s.npy' % t, C[S])
