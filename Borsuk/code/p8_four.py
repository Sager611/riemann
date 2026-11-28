# 4-contact test: H(u) ∩ B(P4, 1) with a random admissible 4th contact P4 below the triangle (x3 < 0, |P4| = sqrt(rho),
# |P4 - P_i| <= 1, 0 in conv{P1..P4} projected appropriately).  Strict CEGAR (tabu) at several thresholds.
import numpy as np, sys, json
sys.path.insert(0, 'Borsuk/code')
from p5_T import H_params
from p6_survey import boundary, cegar
u = json.loads(sys.argv[1]); n = int(sys.argv[2]); seed = int(sys.argv[3]); thrs = [float(x) for x in sys.argv[4].split(',')]
rng = np.random.default_rng(seed); rho, P = H_params(u); R = np.sqrt(rho); cuts = [(np.array([0, 0, 0, 1.0]), 0.5)]
done = 0
while done < n:
    v = rng.normal(size=4); v /= np.linalg.norm(v); p4 = R * v
    if p4[2] >= 0 or np.any(((P - p4) ** 2).sum(1) > 1) or p4[3] > 0.5: continue
    Q = np.vstack([P, p4]); A = np.vstack([Q.T, np.ones(4)]); lam, *_ = np.linalg.lstsq(A, np.r_[np.zeros(4), 1], rcond=None)
    inside = np.allclose(A @ lam, np.r_[np.zeros(4), 1], atol=1e-9) and np.all(lam > 0)
    C = boundary(rho, Q, cuts, 8000, rng); res = {}
    for t in thrs:
        r = cegar(C, t, rng=rng); res[t] = r[0]
        if r[0] == 'COL': break
    print(json.dumps({'u': u, 'p4': np.round(p4, 3).tolist(), '0 in conv4': bool(inside), 'res': res}), flush=True); done += 1
