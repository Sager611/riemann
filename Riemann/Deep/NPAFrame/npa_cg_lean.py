#!/usr/bin/env python3
# npaframe: copied unmodified from NPAFals/npa_cg_lean.py; imports the local npa_cg.py (with T<n0> sets).
# npafals CHANGE vs ExploreNPA/npa_cg.py: memory-lean column generation (random initial on-line columns + big-M slacks, rays priced
# as zero-mass columns instead of all appended, phase 1 on-line only, phase 2 all configurations), so N = 7 fits under 290 MB;
# candidate columns copied (jobs after 23:33 UTC; before, views kept whole pricing blocks alive: N = 7 L4 off peaked at 287 MB).
# usage: npa_cg_lean.py N M s1,s2,... SETS world [-rays]   (SETS as npa_cg.py rowspec: P P3 L4 Q4 P3o L4o L4c)
import sys, time, resource, os
import numpy as np
from scipy.optimize import linprog
from npa_cg import rowspec, feats, blocks, rays

def solve(A, mass, c, b, BIG=50.0):
    R, n = A.shape
    Aeq = np.vstack([np.hstack([A, np.eye(R), -np.eye(R)]), np.concatenate([mass, np.zeros(2 * R)])[None, :]])
    cc = np.concatenate([c, BIG * np.ones(2 * R)]); beq = np.concatenate([b, [1.0]])
    cs = np.maximum(1.0, np.abs(Aeq).max(axis=0))
    r = linprog(cc / cs, A_eq=Aeq / cs, b_eq=beq, bounds=(0, None), method='highs-ds')
    x = r.x / cs
    return r.fun, x[:n], r.eqlin.marginals, float(x[n:].sum())

def main():
    N = int(sys.argv[1]); M = int(sys.argv[2]); depths = [float(t) for t in sys.argv[3].split(',')]
    sets = sys.argv[4]; world = sys.argv[5]; userays = '-rays' in sys.argv and world == 'off'
    spec = rowspec(N, sets); b = np.array([t for (_, _, t, _) in spec]); R = len(spec); t0 = time.time()
    RAYS = rays(N, spec) if userays else np.zeros((R, 0))
    print(f"N={N} M={M} sets={sets} world={world} rows={R} rays={RAYS.shape[1]}", flush=True)
    P0, ns0, lab0 = next(blocks(N, M, depths, 'on'))
    rng = np.random.default_rng(1); pick = rng.choice(len(ns0), size=min(3000, len(ns0)), replace=False)
    F0 = feats(P0[pick], spec)
    A = F0; c = ns0[pick] / N; mass = np.ones(len(pick)); kinds = ['on'] * len(pick)
    for phase in (['on'] if world == 'on' else ['on', world]):
        it = 0
        while True:
            it += 1
            val, x, y, sl = solve(A, mass, c, b)
            yr, ym = y[:R], y[R]
            cand = []; nconf = 0; minrc = 0.0
            for P, ns, lab in blocks(N, M, depths, phase):
                F = feats(P, spec); nconf += len(ns)
                sc = np.maximum(1.0, np.abs(F).max(axis=0))
                rc = (ns / N - yr @ F - ym) / sc
                idx = np.where(rc < -1e-10)[0]
                if len(idx):
                    minrc = min(minrc, rc[idx].min()); o = idx[np.argsort(rc[idx])[:300]]
                    cand += [(rc[i], F[:, i].copy(), ns[i] / N, 1.0, 'off' if lab(i)[1] else 'on') for i in o]   # copy: a view kept every block's F alive (287 MB at N = 7)
                del F, P
            if phase != 'on' and RAYS.shape[1]:
                sc = np.maximum(1.0, np.abs(RAYS).max(axis=0)); rc = (-(yr @ RAYS)) / sc
                idx = np.where(rc < -1e-10)[0]
                if len(idx):
                    minrc = min(minrc, rc[idx].min()); o = idx[np.argsort(rc[idx])[:200]]
                    cand += [(rc[i], RAYS[:, i], 0.0, 0.0, 'ray') for i in o]
            print(f"{phase} it={it} val={val:.6f} slack={sl:.1e} cols={A.shape[1]} priced={nconf} minrc={minrc:.2e} t={time.time()-t0:.0f}s", flush=True)
            if not cand or it > 80: break
            cand.sort(key=lambda z: z[0]); cand = cand[:600]
            A = np.hstack([A, np.array([z[1] for z in cand]).T]); c = np.concatenate([c, [z[2] for z in cand]])
            mass = np.concatenate([mass, [z[3] for z in cand]]); kinds += [z[4] for z in cand]
            if A.shape[1] > 9000:
                kp = np.unique(np.concatenate([np.where(x > 1e-13)[0], np.arange(len(x), A.shape[1])]))
                A = A[:, kp]; c = c[kp]; mass = mass[kp]; kinds = [kinds[i] for i in kp]
    supp = np.where(x > 1e-12)[0]
    sh = {k: float(sum(x[i] for i in supp if kinds[i] == k)) for k in ('on', 'off', 'ray')}
    print(f"RESULT N={N} M={M} sets={sets} world={world} rays={userays}: minSimpleOn={val:.6f} slack={sl:.1e} supp={len(supp)} "
          f"mass on/off={sh['on']:.4f}/{sh['off']:.2e} raycols={sum(1 for i in supp if kinds[i]=='ray')} "
          f"maxrss={resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/2**20:.0f}MB t={time.time()-t0:.0f}s", flush=True)

if __name__ == '__main__':
    main()
