# beta5 brackets of Tolmachev's H(u) on a list of u (CEGAR with kissat fallback): UNSAT = certified lower bound for the
# sample (hence for H(u)), COL = 5-colouring of the whole sample (upper bound for the sample only).
import numpy as np, sys, json
sys.path.insert(0, 'Borsuk/code')
from p6_survey import boundary, cegar
from p6_pd import body
us = json.loads(sys.argv[1]); M, seed = int(sys.argv[2]), int(sys.argv[3]); thrs = [float(x) for x in sys.argv[4].split(',')]
for u in us:
    rng = np.random.default_rng(seed); rho, P, cuts = body({'u': u})
    C = boundary(rho, P, cuts, M, rng); res = {}
    for t in thrs:
        r = cegar(C, t, TL=1800, rng=rng); res[t] = r[0] + ':' + str(r[1])
        if r[0] == 'COL': break
    print(json.dumps({'u': u, 'rho': round(float(rho), 4), 'res': res}), flush=True)
