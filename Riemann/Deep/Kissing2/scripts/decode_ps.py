#!/usr/bin/env python3
"""decode_ps.py -- structure of a PackingStar-family kissing configuration in dimension 24+k (k = 3..7).

Input: an .npz with integer 'points' (PackingStar / KissingNumbers convention, squared norm 32 or 128).
Prints, with exact integer arithmetic:
  * layers by (body norm, tail norm) in the norm-32 scale (body = first 24 coordinates);
  * equator = body norm 32, tail 0: checks it is a subset of the Leech minimal shell L_min (Golay code read
    from the configuration's own (2^8) vectors), deleted set D = L_min minus equator;
  * caps = body norm 24, tail norm 8: every body is the midpoint (a+a')/2 of a Leech EDGE (a.a' = 16); the two
    endpoints are found from the class of 2b mod 2L (fast) and checked to lie in D;
  * incidences r_a = number of caps containing a (Lemma: r_a <= 8, see REPORT.md);
  * connected components of the edge graph (squares = 4-cycles) and square GROUPS sharing a direction pair (u,v).
Usage: python3 decode_ps.py config.npz      (NumPy; ~10-30 s)
"""
import sys, collections, random
import numpy as np

def golay_from_config(E):
    oct_ = set()
    for v in E:
        a = np.abs(v)
        if (a == 2).sum() == 8: oct_.add(int(sum(1 << i for i in range(24) if a[i] == 2)))
    piv = {}
    for o in oct_:
        x = o
        for p in sorted(piv, reverse=True):
            if x >> p & 1: x ^= piv[p]
        if x: piv[x.bit_length() - 1] = x
    words = [0]
    for r in piv.values(): words += [w ^ r for w in words]
    assert len(words) == 4096
    return words

def leech_min(words):
    mins = []
    for w in words:
        if bin(w).count('1') == 8:
            p = [i for i in range(24) if w >> i & 1]
            for m in range(256):
                if bin(m).count('1') % 2: continue
                v = [0] * 24
                for k, i in enumerate(p): v[i] = -2 if m >> k & 1 else 2
                mins.append(v)
    for a in range(24):
        for b in range(a + 1, 24):
            for s in range(4):
                v = [0] * 24; v[a] = -4 if s & 1 else 4; v[b] = -4 if s & 2 else 4; mins.append(v)
    for w in words:
        for a in range(24):
            v = [(3 if w >> i & 1 else -3) if i == a else (-1 if w >> i & 1 else 1) for i in range(24)]
            if sum(v) % 8 == 4: mins.append(v)
    U = np.array(mins, dtype=np.int64); assert U.shape == (196560, 24)
    return U

def hnf(vecs):
    basis = [None] * 24
    for r in vecs:
        r = list(map(int, r))
        for i in range(24):
            if r[i] == 0: continue
            if basis[i] is None:
                if r[i] < 0: r = [-x for x in r]
                basis[i] = r; break
            b = basis[i]
            while r[i] != 0:
                q = b[i] // r[i]; b = [x - q * y for x, y in zip(b, r)]; b, r = r, b
            if b[i] < 0: b = [-x for x in b]
            basis[i] = b
            if not any(r): break
    for i in range(24):
        for k in range(i):
            q = basis[k][i] // basis[i][i]
            if q: basis[k] = [x - q * y for x, y in zip(basis[k], basis[i])]
    return np.array(basis, dtype=np.int64)

def canon(X, H):
    X = np.array(X, dtype=np.int64).copy()
    for i in range(24):
        X -= np.outer(np.floor_divide(X[:, i], H[i, i]), H[i])
    return X

def load_points(path):
    """integer points of a KissingNumbers / PackingStar archive ('points', or 'rational' with zero 'radical')"""
    z = np.load(path)
    if 'points' in z.files: return z['points'].astype(np.int64)
    assert 'rational' in z.files and not np.any(z['radical']), 'only integer certificates are supported'
    return z['rational'].astype(np.int64)

def cl(v):
    v = tuple(int(x) for x in v); nz = next(x for x in v if x); return v if nz > 0 else tuple(-x for x in v)

def main(path):
    P = load_points(path)
    if (P * P).sum(1)[0] == 128:
        assert not np.any(P % 2); P //= 2
    n, d = P.shape; k = d - 24
    B, T = P[:, :24], P[:, 24:]
    bn, tn = (B * B).sum(1), (T * T).sum(1)
    assert np.all(bn + tn == 32)
    print(path, 'points', n, 'dimension', d)
    for key, c in sorted(collections.Counter(zip(bn.tolist(), tn.tolist())).items()):
        print('  layer body %2d tail %2d : %d' % (key[0], key[1], c))
    E = B[bn == 32]; assert np.all(tn[bn == 32] == 0)
    U = leech_min(golay_from_config(E))
    key = {tuple(x): i for i, x in enumerate(U.tolist())}
    inE = np.zeros(len(U), dtype=bool)
    for x in E.tolist(): inE[key[tuple(x)]] = True          # KeyError if an equator point is not in L_min
    print('  equator subset of L_min: yes;  deleted |D| =', int((~inE).sum()))
    random.seed(1); H = hnf(U[random.sample(range(len(U)), 400)])
    assert np.prod(np.diag(H)) == 2 ** 36
    H2 = 2 * H
    cls = {tuple(r): i for i, r in enumerate(canon(U, H2).tolist())}
    ci = np.flatnonzero(bn == 24); Bc = B[ci]
    w = [cls[tuple(r)] for r in canon(2 * Bc, H2).tolist()]    # class of 2b mod 2L is a minimal vector w
    ends = np.array([(key[tuple(((2 * b + U[j]) // 2).tolist())], key[tuple(((2 * b - U[j]) // 2).tolist())])
                     for b, j in zip(Bc, w)])
    assert not inE[ends].any()
    assert np.all((U[ends[:, 0]] * U[ends[:, 1]]).sum(1) == 16)
    print('  caps', len(ci), ': every body is an edge midpoint (a+a\')/2, a.a\'=16, both endpoints deleted')
    inc = collections.Counter(ends.ravel().tolist())
    print('  incidences r_a per deleted vertex:', sorted(collections.Counter(inc.values()).items()),
          ' caps/|D| = %.3f' % (len(ci) / max(1, int((~inE).sum()))))
    edges = {}
    for r, (x, y) in enumerate(ends): edges.setdefault((min(x, y), max(x, y)), set()).add(cl(T[ci[r]]))
    adj = collections.defaultdict(list)
    for (x, y) in edges: adj[x].append(y); adj[y].append(x)
    seen = set(); comp = collections.Counter(); groups = collections.Counter()
    for s in list(adj):
        if s in seen: continue
        st = [s]; seen.add(s); vs = []
        while st:
            x = st.pop(); vs.append(x)
            for y in adj[x]:
                if y not in seen: seen.add(y); st.append(y)
        ne = sum(len(adj[x]) for x in vs) // 2
        comp[(len(vs), ne)] += 1
        if len(vs) == 4 and ne == 4:
            dirs = sorted(set(cl(U[a] - U[b]) for a in vs for b in adj[a]))
            groups[tuple(dirs)] += 1
    print('  edges', len(edges), ' components (vertices,edges):count', sorted(comp.items(), key=lambda t: -t[1])[:8])
    print('  squares', sum(groups.values()), 'in', len(groups), 'direction-pair groups; group sizes',
          sorted(collections.Counter(groups.values()).items()))

if __name__ == '__main__':
    main(sys.argv[1])
