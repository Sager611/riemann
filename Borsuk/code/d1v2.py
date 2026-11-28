#!/usr/bin/env python3
"""D1 (v2): subsets of uniform 4-polytope vertex sets and root systems, local check per threshold.

P is an orbit of a reflection group (vertex-transitive under isometries).  A subset S with diam(S)^2 = s
can be moved by an isometry of P to contain the base point x0; then S lies in L_s = {p : d^2(x0, p) <= s} and its
diameter graph is an induced subgraph of H_s[L_s] (pairs of L_s at squared distance exactly s).  So if every
H_s[L_s] is 5-colourable, no subset of P is a counterexample.  Colourings come from max degree < 5, tabu search
or SAT and are checked edge by edge; arithmetic is exact (integer coordinates in Z[w]).  A threshold whose local
graph is not 5-colourable goes to the CEGAR of d1_local.py.
Usage: python d1v2.py [group ...]   groups: ROOTS A4 B4 F4 H4
"""
import sys, json, time
import numpy as np
sys.path.insert(0, 'Borsuk/code'); sys.path.insert(0, 'Borsuk/code')
import d1
from tabu import colour_tabu

def local_graph(X, key, fv):
    n = len(X); DA, DB = d1.sqd_block(X, 0, 1); kk = d1.keys_of(DA[0], DB[0])
    sv = fv(key); vals = np.array([fv(int(t)) for t in kk])
    tol = 1e-9 * max(1.0, abs(sv))
    L = np.flatnonzero(vals <= sv + tol)
    assert all(int(kk[j]) == key or abs(vals[j] - sv) > tol for j in L), 'float tie'
    Y = X[L]; m = len(L); U = []; V = []
    for i0 in range(0, m, 128):
        DA, DB = d1.sqd_block(Y, i0, min(m, i0 + 128)); K2 = d1.keys_of(DA, DB)
        ii, jj = np.nonzero(K2 == key); ii = ii + i0; s = jj > ii; U.append(ii[s]); V.append(jj[s])
    return m, np.concatenate(U).astype(np.int32), np.concatenate(V).astype(np.int32)

def check(name, field, pts, out):
    d1.set_field(field); X = d1.to_int(pts); n = len(X); t0 = time.time()
    DA, DB = d1.sqd_block(X, 0, 1); kk = d1.keys_of(DA[0], DB[0])[1:]
    vals, cnt = np.unique(kk, return_counts=True)
    wf = d1.WF[field]; fv = lambda key: (lambda ab: ab[0] + ab[1] * wf)(d1.unkey(key))
    st = {}; bad = []
    for key, c in zip(vals.tolist(), cnt.tolist()):
        if c < 5: st['deg<5'] = st.get('deg<5', 0) + 1; continue
        m, U, V = local_graph(X, key, fv)
        deg = np.bincount(np.concatenate([U, V]), minlength=m)
        if deg.max() < 5: res = 'deg<5'
        else:
            res = 'hard'; col = None
            for sd in range(1, 5):
                c_, conf = colour_tabu(m, U, V, 5, iters=500000 * sd, seed=sd)
                if conf == 0: res = 'tabu'; col = c_; break
            if col is not None:
                c_ = np.asarray(col); assert ((c_ >= 0) & (c_ < 5)).all() and (c_[U] != c_[V]).all(), 'bad colouring'
            else: bad.append({'key': key, 's': round(fv(key), 6), 'L': m, 'edges': len(U), 'st': res})
        st[res] = st.get(res, 0) + 1
    rec = {'set': name, 'n': n, 'thresholds': len(vals), **st, 'bad': bad, 'sec': round(time.time() - t0, 1)}
    print(json.dumps(rec), flush=True); out.write(json.dumps(rec) + '\n'); out.flush()

if __name__ == '__main__':
    groups = sys.argv[1:] or ['ROOTS', 'A4', 'B4', 'F4', 'H4']
    with open('Borsuk/.scratch/d1v2_results.jsonl', 'a') as out:
        for g in groups:
            for name, field, pts, vt in d1.sets_for(g):
                if not vt:   # unions of two orbits: run the global check of d1.py
                    d1.check_set(name, field, pts, vt, out); continue
                check(name, field, pts, out)
