# End-to-end test: random all-contact configuration -> optimised rotated simplex cones -> cert2 on each of the 5 cells.
import numpy as np, sys, json, subprocess
sys.path.insert(0, 'Borsuk/code')
from p6_survey import random_config, boundary
from p6_cones import optimise, rot, SIMP, pieces_diam
typ, seed, target = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3])
rng = np.random.default_rng(seed); rho, P, cuts, lam = random_config(typ, rng)
X = boundary(rho, P, cuts, 3000, rng); f, z = optimise(X, rng, iters=100)
Q = rot(z[:6]); ap = z[6:]; D = SIMP @ Q.T
Xd = boundary(rho, P, cuts, 30000, rng); fd = pieces_diam(Xd, Q, ap)
R = np.sqrt(rho)
txt = '%d\n' % 5
for j in range(5):
    hs = [(D[k] - D[j], (D[k] - D[j]) @ ap) for k in range(5) if k != j] + [(n, c) for n, c in cuts]
    txt += '%d %d\n' % (1 + len(P), len(hs))
    txt += '0 0 0 0 %.18f\n' % rho
    for p in P: txt += ' '.join('%.18f' % x for x in p) + ' 1\n'
    for a, b in hs: txt += ' '.join('%.18f' % x for x in a) + ' %.18f\n' % b
txt += '5\n' + ''.join('%d %d %.6f %.6f 3000000\n' % (j, j, target, R + 1e-9) for j in range(5))
open('/tmp/e2e.txt', 'w').write(txt)
out = subprocess.run(['bash', '-c', 'Borsuk/.scratch/cert2 < /tmp/e2e.txt'], capture_output=True, text=True).stdout
print(json.dumps({'type': typ, 'rho': round(rho, 4), 'phi_opt': round(f, 4), 'phi_dense': round(fd, 4), 'cert2': out.split('\n')[:5]}), flush=True)
