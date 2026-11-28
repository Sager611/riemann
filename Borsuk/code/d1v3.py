#!/usr/bin/env python3
"""D1 (v3): the local check of d1v2.py for large vertex-transitive orbits, with one streaming pass over all pairs.

For each threshold s with multiplicity >= 5 at the base point x0, the local graph H_s[L_s] has the edges {p, q} with
d^2(p, q) = s and d^2(x0, p), d^2(x0, q) <= s.  Pairs are streamed in row blocks; thresholds are processed in groups
that fit in memory.  Exact integer keys; colourings (tabu) checked edge by edge; uncoloured thresholds are listed as
'hard' (input for d1_cegar.py).  Usage: python d1v3.py NAME [NAME ...]   (names as in d1.sets_for, e.g. H4:0123)
"""
import sys, json, time
import numpy as np
sys.path.insert(0, 'Borsuk/code'); sys.path.insert(0, 'Borsuk/code')
import d1
from tabu import colour_tabu

def run(name, field, pts, out):
    d1.set_field(field); X = d1.to_int(pts); n = len(X); t0 = time.time()
    wf = d1.WF[field]; fv = lambda key: (lambda ab: ab[0] + ab[1] * wf)(d1.unkey(key))
    DA, DB = d1.sqd_block(X, 0, 1); k0 = d1.keys_of(DA[0], DB[0]); k0[0] = -1
    vals, cnt = np.unique(k0[1:], return_counts=True)
    order = np.argsort([fv(int(v)) for v in vals]); fl = np.array([fv(int(v)) for v in vals])[order]
    assert np.all(np.diff(fl) > 1e-9 * np.maximum(1, np.abs(fl[1:]))), 'float ordering not safe'
    sk = vals[order]; cnt = cnt[order]            # keys sorted by value; rank = index
    pos = np.argsort(sk); skey_sorted = sk[pos]    # integer-sorted keys -> rank
    def rank(K):
        i = np.searchsorted(skey_sorted, K); i = np.clip(i, 0, len(sk) - 1)
        ok = skey_sorted[i] == K; r = np.where(ok, pos[i], -1); return r
    rad = rank(k0); rad[0] = -1                     # x0 itself: radius rank -1 (inside every ball)
    need = np.flatnonzero(cnt >= 5)
    st = {'deg<5': int(len(sk) - len(need))}; bad = []
    est = {int(t): 0 for t in need}
    # group thresholds so that each group holds <= ~6M local edges (estimate by n * multiplicity / 2)
    groups = []; cur = []; tot = 0
    for t in need:
        e = n * int(cnt[t]) // 2
        if cur and tot + e > 6_000_000: groups.append(cur); cur = []; tot = 0
        cur.append(int(t)); tot += e
    if cur: groups.append(cur)
    for G in groups:
        Gs = np.array(G); inG = np.zeros(len(sk), bool); inG[Gs] = True
        E = {t: [] for t in G}
        for i0 in range(0, n, 64):
            i1 = min(n, i0 + 64); DA, DB = d1.sqd_block(X, i0, i1); K = d1.keys_of(DA, DB)
            R = rank(K.ravel()).reshape(K.shape)
            ii, jj = np.nonzero((R >= 0) & inG[np.maximum(R, 0)])
            ii = ii + i0; m = jj > ii; ii, jj = ii[m], jj[m]; r = R[ii - i0, jj]
            keep = (rad[ii] <= r) & (rad[jj] <= r); ii, jj, r = ii[keep], jj[keep], r[keep]
            for t in np.unique(r): s_ = r == t; E[int(t)].append(np.stack([ii[s_], jj[s_]], 1).astype(np.int32))
        for t in G:
            if not E[t]: st['deg<5'] += 1; continue
            e = np.concatenate(E[t]); del E[t]
            verts = np.unique(e); loc = -np.ones(n, np.int64); loc[verts] = np.arange(len(verts))
            U = loc[e[:, 0]].astype(np.int32); V = loc[e[:, 1]].astype(np.int32); m_ = len(verts)
            if np.bincount(np.concatenate([U, V]), minlength=m_).max() < 5: st['deg<5'] += 1; continue
            res = 'hard'
            for sd in range(1, 5):
                col, conf = colour_tabu(m_, U, V, 5, iters=500000 * sd, seed=sd)
                if conf == 0:
                    assert (col[U] != col[V]).all(); res = 'tabu'; break
            if res == 'hard': bad.append({'key': int(sk[t]), 's': round(float(fl[t]), 6), 'edges': len(U), 'st': 'hard'})
            st[res] = st.get(res, 0) + 1
    rec = {'set': name, 'n': n, 'thresholds': len(sk), **st, 'bad': bad, 'sec': round(time.time() - t0, 1)}
    print(json.dumps(rec), flush=True); out.write(json.dumps(rec) + '\n'); out.flush()

if __name__ == '__main__':
    want = set(sys.argv[1:])
    with open('Borsuk/.scratch/d1v2_results.jsonl', 'a') as out:
        for g in sorted({w.split(':')[0] for w in want}):
            for name, field, pts, vt in d1.sets_for(g):
                if name in want: run(name, field, pts, out)
