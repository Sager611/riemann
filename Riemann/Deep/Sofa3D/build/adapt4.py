"""adapt4.py -- lean quadtree-adaptive certification (per-level arrays, vectorized selection).
Lattice point i -> (i + off)/N, N = nc0 * 2^L; cell (level, j, k) spans lattice [j w, (j+1) w] x [k w, (k+1) w], w = 2^(L-level)."""
import math, time
import numpy as np
import sofa2 as S
from adapt2 import build_match, build_cut

def certify(C, nc0, Lmax, tau_rel, tau_abs, off=0.38196601125, gmerge=0.02, log=print, jr=None, kr=None):
    N = nc0 << Lmax; N1 = N + 1
    pos = lambda i: (np.asarray(i, float) + off) / N
    lines = {}
    st = dict(lines_s=0.0, pieces_s=0.0, build_s=0.0, nlines=0, npieces=0)
    def ensure(keys):
        miss = [k for k in keys if k not in lines]
        if not miss: return
        arr = np.array(miss, dtype=np.int64)
        P = np.zeros((len(arr), 10)); P[:, 0] = pos(arr // N1); P[:, 1] = pos(arr % N1); P[:, 4] = -4.0
        cnt, iv, dt, info = S.run(C, P)
        st['lines_s'] += dt; st['nlines'] += len(arr)
        o = np.zeros(len(cnt) + 1, np.int64); np.cumsum(cnt, out=o[1:])
        X = -4.0 + 4.0 * iv
        for i, k in enumerate(miss): lines[k] = X[o[i]:o[i + 1]]
    jr = jr or (0, nc0); kr = kr or (0, nc0)
    J0, K0 = np.meshgrid(np.arange(*jr), np.arange(*kr), indexing='ij')
    lev = []          # per level: dict(J, K, cert, refined, P, own, cnt, iv)
    J = J0.ravel(); K = K0.ravel()
    for level in range(Lmax + 1):
        if len(J) == 0: break
        w = 1 << (Lmax - level)
        keys = set()
        for a, b in ((0, 0), (1, 0), (0, 1), (1, 1)):
            keys.update(((J + a) * w * N1 + (K + b) * w).tolist())
        ensure(keys)
        t0 = time.time()
        hy = hz = w / (2.0 * N); area = 4 * hy * hz
        Pl = []; own = []; cand = []; ideal = np.zeros(len(J))
        for ci in range(len(J)):
            j = int(J[ci]); k = int(K[ci])
            cl = [lines[(j * w) * N1 + k * w], lines[((j + 1) * w) * N1 + k * w], lines[(j * w) * N1 + (k + 1) * w], lines[((j + 1) * w) * N1 + (k + 1) * w]]
            ideal[ci] = area * sum(float((L[:, 1] - L[:, 0]).sum()) for L in cl) / 4
            lists = [S.merge_gaps(L, gmerge) for L in cl]
            if any(len(L) == 0 for L in lists): continue
            yc = float(pos((j + 0.5) * w)); zc = float(pos((k + 0.5) * w))
            ps = build_match(lists, yc, zc, hy, hz)
            if len(ps): Pl.append(ps); own.extend([ci] * len(ps)); cand.extend([0] * len(ps))
            if any(len(L) > 1 for L in lists):
                pc = build_cut(lists, yc, zc, hy, hz)
                if len(pc): Pl.append(pc); own.extend([ci] * len(pc)); cand.extend([1] * len(pc))
        st['build_s'] += time.time() - t0
        cv = np.zeros(len(J)); rec = dict(J=J, K=K, w=w)
        if Pl:
            PS = np.vstack(Pl); del Pl
            own = np.array(own, np.int64); cand = np.array(cand, np.int8)
            cnt, iv, dt, info = S.run(C, PS)
            st['pieces_s'] += dt; st['npieces'] += len(PS)
            s, o = S.piece_len(cnt, iv)
            pv = s * (PS[:, 7] - PS[:, 4]) * area
            cvc = np.zeros((len(J), 2)); np.add.at(cvc, (own, cand), pv)
            win = np.argmax(cvc, 1); cv = cvc.max(1)
            keep = np.nonzero((cnt > 0) & (cand == win[own]))[0]
            kc = cnt[keep]; ko = np.zeros(len(keep) + 1, np.int64); np.cumsum(kc, out=ko[1:])
            sel = np.repeat(o[keep], kc) + (np.arange(int(kc.sum())) - np.repeat(ko[:-1], kc))
            rec.update(P=PS[keep], own=own[keep], cnt=kc, iv=iv[sel])
            del PS, iv, cnt
        loss = ideal - cv
        refined = (loss > tau_abs) & (loss > tau_rel * area) if level < Lmax else np.zeros(len(J), bool)
        rec.update(cert=cv, refined=refined)
        lev.append(rec)
        log('level %d cells %d cert %.6f ideal %.6f refine %d [lines %d %.0fs pieces %d %.0fs build %.0fs]' % (
            level, len(J), cv.sum(), ideal.sum(), int(refined.sum()), st['nlines'], st['lines_s'], st['npieces'], st['pieces_s'], st['build_s']))
        r = np.nonzero(refined)[0]
        J = (2 * J[r][:, None] + np.array([0, 0, 1, 1])[None, :]).ravel()
        K = (2 * K[r][:, None] + np.array([0, 1, 0, 1])[None, :]).ravel()
    lines.clear()
    # bottom-up best
    nlev = len(lev)
    best = [None] * nlev; usekids = [None] * nlev
    for l in range(nlev - 1, -1, -1):
        b = lev[l]['cert'].copy(); uk = np.zeros(len(b), bool)
        r = np.nonzero(lev[l]['refined'])[0]
        if len(r) and l + 1 < nlev:
            ks = best[l + 1].reshape(-1, 4).sum(1)
            better = ks > b[r]
            b[r[better]] = ks[better]; uk[r[better]] = True
        best[l] = b; usekids[l] = uk
    total = math.fsum(best[0].tolist())
    # top-down selection and sub-piece expansion
    sel = np.ones(len(lev[0]['J']), bool)
    subs = []; cells = []
    for l in range(nlev):
        R = lev[l]
        ownsel = sel & ~usekids[l]
        if 'P' in R and len(R['P']):
            m = ownsel[R['own']]
            P = R['P'][m]; c = R['cnt'][m]
            ko = np.zeros(len(R['cnt']) + 1, np.int64); np.cumsum(R['cnt'], out=ko[1:])
            idx = np.nonzero(m)[0]
            ivs = np.concatenate([R['iv'][ko[i]:ko[i + 1]] for i in idx]) if len(idx) else np.zeros((0, 2))
            Pr = np.repeat(P, c, axis=0)
            A = Pr[:, 4:7]; B = Pr[:, 7:10]
            s0 = ivs[:, :1]; s1 = ivs[:, 1:]
            sp = np.concatenate([Pr[:, :4], A + s0 * (B - A), A + s1 * (B - A)], 1)
            ok = ivs[:, 1] > ivs[:, 0]
            ow = np.repeat(R['own'][m], c)
            subs.append(sp[ok]); cells.append(np.stack([np.full(int(ok.sum()), l), R['J'][ow[ok]], R['K'][ow[ok]]], 1))
        if l + 1 < nlev:
            r = np.nonzero(R['refined'])[0]
            sel = np.repeat(sel[r] & usekids[l][r], 4)
    return total, np.vstack(subs), np.vstack(cells).astype(np.int64), st
