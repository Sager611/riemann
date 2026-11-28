# Minimal-ish UNSAT core of the extreme type-3 body at threshold thr: CEGAR UNSAT set, then block-deletion minimisation
# with kissat; report the core's distance profile (how many edges are > 1, i.e. impossible in a diameter-1 set).
import numpy as np, sys, json, time
sys.path.insert(0, 'Borsuk/code')
from p6_survey import boundary, cegar, d2, kissat_colour
rho, M, seed, thr = float(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), float(sys.argv[4])
R = np.sqrt(rho); rng = np.random.default_rng(seed)
th = 2 * np.pi * np.arange(3) / 3; P = R * np.stack([np.cos(th), np.sin(th), 0 * th, 0 * th], 1)
cuts = [(np.array([0, 0, 1.0, 0]), 0.5), (np.array([0, 0, -1.0, 0]), 0.5), (np.array([0, 0, 0, 1.0]), 0.5)]
C = boundary(rho, P, cuts, M, rng)
r = cegar(C, thr, TL=1800, rng=rng, want=True)
print(r[0], r[1], flush=True)
if r[0] != 'UNSAT': sys.exit()
X = r[2]
def edges(Y):
    D = d2(Y, Y); U, V = np.nonzero(np.triu(D >= thr * thr, 1)); return U, V
def unsat(Y):
    U, V = edges(Y); st, _ = kissat_colour(len(Y), U, V, tl=600); return st == 'UNSAT'
keep = np.ones(len(X), bool); t0 = time.time()
for blk in [400, 200, 100, 50, 25, 12, 6, 3, 1]:
    idx = np.flatnonzero(keep); rng.shuffle(idx)
    for s in range(0, len(idx), blk):
        if time.time() - t0 > 6000: break
        trial = keep.copy(); trial[idx[s:s + blk]] = False
        if trial.sum() >= 6 and unsat(X[trial]): keep = trial
    print(json.dumps({'blk': blk, 'core': int(keep.sum()), 'sec': round(time.time() - t0)}), flush=True)
Y = X[keep]; np.save('Borsuk/.scratch/ext3_mincore_%s.npy' % thr, Y)
D = np.sqrt(np.maximum(d2(Y, Y), 0)); iu = np.triu_indices(len(Y), 1); dd = D[iu]; e = dd >= thr
print(json.dumps({'core': len(Y), 'edges': int(e.sum()), 'edges_gt1': int((dd > 1).sum()), 'edges_in_[thr,1]': int((e & (dd <= 1)).sum()),
                  'max_dist': round(float(dd.max()), 4), 'norms': np.round(np.sort(np.linalg.norm(Y, axis=1)), 3).tolist()[:40]}), flush=True)
