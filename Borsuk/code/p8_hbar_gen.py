# Generate a barrier point set for Tolmachev's H(u), u = (0.8,0,0,0): CEGAR UNSAT set at threshold thr (float), then
# rational points (denominator 2^40) shrunk by (1 - 1e-9) towards 0 (keeps membership with margin).  Writes POINTS file.
import numpy as np, sys, json
from fractions import Fraction as Fr
sys.path.insert(0, 'Borsuk/code')
from p6_survey import boundary, cegar
from p6_pd import body
u = json.loads(sys.argv[1]); M, seed, thr, out = int(sys.argv[2]), int(sys.argv[3]), float(sys.argv[4]), sys.argv[5]
rng = np.random.default_rng(seed); rho, P, cuts = body({'u': u}); C = boundary(rho, P, cuts, M, rng)
r = cegar(C, thr, TL=3600, rng=rng, want=True); print(r[0], r[1], flush=True)
if r[0] != 'UNSAT': sys.exit(1)
X = r[2] * (1 - 1e-9); den = 1 << 40
with open(out, 'w') as f:
    f.write('# H(u) barrier points, u = %s, float threshold %s; rational coordinates (denominator 2^40)\n' % (u, thr))
    for x in X: f.write(' '.join(str(Fr(int(round(c * den)), den)) for c in x) + '\n')
print('wrote', len(X), 'points to', out)
