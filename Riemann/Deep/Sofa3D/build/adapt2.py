"""adapt2.py -- quadtree-adaptive certification with sheared prisms; two candidate
constructions per cell ('match' and 'cut'), best one kept (all verified by carve2)."""
import sys, math, time, json
import numpy as np
import sofa2 as S

def plane_fit(a, hy, hz):
    m = a.mean(1)
    u = (a[:, 1] + a[:, 3] - a[:, 0] - a[:, 2]) / (4 * hy)
    v = (a[:, 2] + a[:, 3] - a[:, 0] - a[:, 1]) / (4 * hz)
    t = (a[:, 0] - a[:, 1] - a[:, 2] + a[:, 3]) / 4
    return m, u, v, t

def pieces_from(A, B, yc, zc, hy, hz):
    """A, B: (n,3) planes (c0, slope_u, slope_v); clamp into [-4,0]; keep B>A at corners."""
    A = A.copy(); B = B.copy()
    A[:, 0] = np.maximum(A[:, 0], -4.0 + np.abs(A[:, 1]) * hy + np.abs(A[:, 2]) * hz)
    B[:, 0] = np.minimum(B[:, 0], -(np.abs(B[:, 1]) * hy + np.abs(B[:, 2]) * hz))
    n = len(A)
    P = np.concatenate([np.tile([yc, zc, hy, hz], (n, 1)), A, B], 1)
    ok = np.ones(n, bool)
    for sy in (-1, 1):
        for sz in (-1, 1):
            ok &= (B[:, 0] - A[:, 0]) + (B[:, 1] - A[:, 1]) * sy * hy + (B[:, 2] - A[:, 2]) * sz * hz > 1e-9
    return P[ok]

def fit_lo(a, hy, hz):
    m, u, v, t = plane_fit(a, hy, hz); return np.stack([m + np.abs(t), u, v], 1)

def fit_hi(b, hy, hz):
    m, u, v, t = plane_fit(b, hy, hz); return np.stack([m - np.abs(t), u, v], 1)

def build_match(lists, yc, zc, hy, hz):
    if all(len(L) == 1 for L in lists):
        a = np.array([[L[0, 0] for L in lists]]); b = np.array([[L[0, 1] for L in lists]])
        if b.min() - a.max() <= 0: return np.zeros((0, 10))
        return pieces_from(fit_lo(a, hy, hz), fit_hi(b, hy, hz), yc, zc, hy, hz)
    tups = set()
    for ci in range(4):
        for I in range(len(lists[ci])):
            lo, hi = lists[ci][I]; tup = []
            for cj in range(4):
                if cj == ci: tup.append(I); continue
                L = lists[cj]
                o = np.minimum(L[:, 1], hi) - np.maximum(L[:, 0], lo)
                bi = int(np.argmax(o))
                if o[bi] <= 0: tup = None; break
                tup.append(bi)
            if tup is not None: tups.add(tuple(tup))
    if not tups: return np.zeros((0, 10))
    tups = sorted(tups)
    a = np.array([[lists[c][t[c], 0] for c in range(4)] for t in tups]); b = np.array([[lists[c][t[c], 1] for c in range(4)] for t in tups])
    P = pieces_from(fit_lo(a, hy, hz), fit_hi(b, hy, hz), yc, zc, hy, hz)
    if len(P) <= 1: return P
    xlo = P[:, 4] - np.abs(P[:, 5]) * hy - np.abs(P[:, 6]) * hz; xhi = P[:, 7] + np.abs(P[:, 8]) * hy + np.abs(P[:, 9]) * hz
    order = np.argsort(-(P[:, 7] - P[:, 4])); keep = []
    for i in order:
        if all(xhi[i] < xlo[j] or xhi[j] < xlo[i] for j in keep): keep.append(i)
    return P[sorted(keep)]

def build_cut(lists, yc, zc, hy, hz):
    a = np.array([[L[0, 0] for L in lists]]); b = np.array([[L[-1, 1] for L in lists]])
    gaps = []
    for L in lists:
        for i in range(len(L) - 1): gaps.append((L[i, 1], L[i + 1, 0]))
    gaps.sort()
    G = []
    for g in gaps:
        if G and g[0] <= G[-1][1]: G[-1][1] = max(G[-1][1], g[1])
        else: G.append([g[0], g[1]])
    Al = fit_lo(a, hy, hz)[0]; Bh = fit_hi(b, hy, hz)[0]
    bounds_lo = [Al] + [np.array([g[1], 0.0, 0.0]) for g in G]
    bounds_hi = [np.array([g[0], 0.0, 0.0]) for g in G] + [Bh]
    return pieces_from(np.array(bounds_lo), np.array(bounds_hi), yc, zc, hy, hz)

def certify(C, nc0, Lmax, tau_rel, tau_abs, ts=(), gmerge=0.02, log=print, use_cut=True):
    N = nc0 << Lmax
    lines = {}
    st = dict(lines_s=0.0, pieces_s=0.0, build_s=0.0, nlines=0, npieces=0)
    def ensure(keys):
        miss = [k for k in keys if k not in lines]
        if not miss: return
        arr = np.array(miss, dtype=np.int64)
        P = np.zeros((len(arr), 10)); P[:, 0] = arr[:, 0] / N; P[:, 1] = arr[:, 1] / N; P[:, 4] = -4.0
        cnt, iv, dt, info = S.run(C, P, ts=ts)
        st['lines_s'] += dt; st['nlines'] += len(arr)
        off = np.zeros(len(cnt) + 1, np.int64); np.cumsum(cnt, out=off[1:])
        X = -4.0 + 4.0 * iv
        for i, k in enumerate(miss): lines[k] = X[off[i]:off[i + 1]]
    cert = {}; kids = {}; store = {}
    active = [(0, j, k) for j in range(nc0) for k in range(nc0)]
    for level in range(Lmax + 1):
        if not active: break
        w = 1 << (Lmax - level)
        keys = set()
        for (_, j, k) in active:
            keys.update([(j * w, k * w), ((j + 1) * w, k * w), (j * w, (k + 1) * w), ((j + 1) * w, (k + 1) * w)])
        ensure(keys)
        t0 = time.time()
        hy = hz = w / (2.0 * N); area = 4 * hy * hz
        Pl = []; own = []; cand = []; ideal = np.zeros(len(active))
        for ci, (_, j, k) in enumerate(active):
            cl = [lines[(j * w, k * w)], lines[((j + 1) * w, k * w)], lines[(j * w, (k + 1) * w)], lines[((j + 1) * w, (k + 1) * w)]]
            ideal[ci] = area * sum(float((L[:, 1] - L[:, 0]).sum()) for L in cl) / 4
            lists = [S.merge_gaps(L, gmerge) for L in cl]
            if any(len(L) == 0 for L in lists): continue
            yc = (j + 0.5) * w / N; zc = (k + 0.5) * w / N
            ps = build_match(lists, yc, zc, hy, hz)
            if len(ps): Pl.append(ps); own.extend([ci] * len(ps)); cand.extend([0] * len(ps))
            if use_cut and any(len(L) > 1 for L in lists):
                pc = build_cut(lists, yc, zc, hy, hz)
                if len(pc): Pl.append(pc); own.extend([ci] * len(pc)); cand.extend([1] * len(pc))
        st['build_s'] += time.time() - t0
        cv = np.zeros(len(active))
        if Pl:
            PS = np.vstack(Pl); own = np.array(own); cand = np.array(cand)
            cnt, iv, dt, info = S.run(C, PS, ts=ts)
            st['pieces_s'] += dt; st['npieces'] += len(PS)
            s, off = S.piece_len(cnt, iv)
            pv = s * (PS[:, 7] - PS[:, 4]) * area
            cvc = np.zeros((len(active), 2)); np.add.at(cvc, (own, cand), pv)
            win = np.argmax(cvc, 1); cv = cvc.max(1)
            for i in range(len(PS)):
                if cnt[i] == 0 or cand[i] != win[own[i]]: continue
                store.setdefault(active[own[i]], []).append((PS[i], iv[off[i]:off[i + 1]].copy()))
        nxt = []; nref = 0
        for ci, cell in enumerate(active):
            cert[cell] = cv[ci]
            loss = ideal[ci] - cv[ci]
            if level < Lmax and loss > tau_abs and loss > tau_rel * area:
                nref += 1; _, j, k = cell
                ch = [(level + 1, 2 * j + dj, 2 * k + dk) for dj in (0, 1) for dk in (0, 1)]
                kids[cell] = ch; nxt.extend(ch)
        log('level %d cells %d cert %.6f ideal %.6f refine %d  [lines %d %.0fs pieces %d %.0fs build %.0fs]' % (
            level, len(active), cv.sum(), ideal.sum(), nref, st['nlines'], st['lines_s'], st['npieces'], st['pieces_s'], st['build_s']))
        active = nxt
    best = {}; choice = {}
    def solve(cell):
        own_v = cert.get(cell, 0.0)
        if cell in kids:
            sv = sum(solve(c) for c in kids[cell])
            if sv > own_v:
                choice[cell] = 'kids'; return sv
        choice[cell] = 'own'; return own_v
    total = math.fsum(solve((0, j, k)) for j in range(nc0) for k in range(nc0))
    final = []
    def collect(cell):
        if choice.get(cell) == 'kids':
            for c in kids[cell]: collect(c)
        else: final.extend(store.get(cell, []))
    for j in range(nc0):
        for k in range(nc0): collect((0, j, k))
    return total, final, st

def subpieces(final, shrink=0.0):
    out = []
    for P, iv in final:
        A = P[[4, 5, 6]]; B = P[[7, 8, 9]]
        for s0, s1 in iv:
            if s1 - s0 <= 2 * shrink: continue
            s0 += shrink; s1 -= shrink
            out.append(np.concatenate([P[:4], A + s0 * (B - A), A + s1 * (B - A)]))
    return np.array(out)
