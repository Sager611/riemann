#!/usr/bin/env python3
"""discover.py -- how the sublattice M was found (NumPy; about one minute).  Not needed for verification.

1. Build the Leech lattice L (x sqrt8 scale) from golay_basis.txt and its 196560 minimal vectors.
2. Conway's elements xi_T = eps_{T1} eta_T (eta_T = reflection I - J/2 on each tetrad of a sextet) together with
   sign changes eps_c generate (part of) Co_0; a random word of order 78 has a power omega of order 3 with trace -12,
   i.e. a fixed-point-free automorphism with omega^2 + omega + 1 = 0 (an Eisenstein structure on L).
3. theta = omega - omega^2 satisfies theta^2 = -3; L/theta L = F_3^12.  The minimal vectors fall into 65520 classes T
   (three vectors u, omega u, omega^2 u per class).  The Fourier transform of the indicator of T takes only the values
   -90 (465920 times) and 639 (65520 times) off 0, so a subspace C of F_3^12 avoiding T has |C^perp| >= 1 + 65520/90
   = 729, i.e. dim C <= 6, with equality iff all nonzero vectors of C^perp are in the -90 set.
4. A greedy random search finds a 6-dimensional C avoiding T at once.  M = preimage of C: index 3^6, minimum 6
   (std scale), 26208 vectors of norm 6 -> 13104 lines.  M = {x in L : F x = 0 mod 3} with F written to F3.txt.
"""
import numpy as np, random, os, sys, itertools
HERE = os.path.dirname(os.path.abspath(__file__))
rows = [l.strip() for l in open(os.path.join(HERE, 'golay_basis.txt')) if l.strip()]
words = [0]
for r in rows:
    b = sum(1 << i for i, c in enumerate(r) if c == '1'); words += [w ^ b for w in words]
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
key = set(map(bytes, (U + 8).astype(np.uint8)))
octads = [w for w in words if bin(w).count('1') == 8]
def xi(T1, k):
    tet = [T1] + [o & ~T1 for o in octads if o & T1 == T1]
    E = np.zeros((24, 24), dtype=np.int64)
    for t in tet:
        idx = [i for i in range(24) if t >> i & 1]
        for i in idx:
            for j in idx: E[i, j] = (8 if i == j else 0) - 4
    Dg = np.eye(24, dtype=np.int64)
    for i in range(24):
        if tet[k] >> i & 1: Dg[i, i] = -1
    return Dg @ E                       # 8*xi, acting on column vectors
def is_aut(G8):
    W = U @ G8.T
    return not np.any(W % 8) and all(bytes(w) in key for w in (W // 8 + 8).astype(np.uint8))
def mul(a, b):
    p = a @ b; assert not np.any(p % 8); return p // 8
random.seed(1)
gens = []
for _ in range(6):
    T1 = sum(1 << p for p in random.sample(range(24), 4)); g = xi(T1, random.randrange(6)); assert is_aut(g); gens.append(g)
for c in random.sample(words, 4):
    gens.append(np.diag([-8 if c >> i & 1 else 8 for i in range(24)]).astype(np.int64))
I8 = 8 * np.eye(24, dtype=np.int64); omega = None
while omega is None:
    g = I8.copy()
    for _ in range(random.randrange(8, 30)): g = mul(g, random.choice(gens))
    h = g.copy(); o = None
    for k in range(1, 131):
        if np.array_equal(h, I8): o = k; break
        h = mul(h, g)
    if o and o % 3 == 0:
        h = I8.copy()
        for _ in range(o // 3): h = mul(h, g)
        if np.trace(h) == -96: omega = h; print('omega = (word of order %d)^%d' % (o, o // 3))
assert is_aut(omega) and np.array_equal(mul(omega, omega) + omega + I8, 0 * I8)
Wr = omega.T                                           # rows: x -> x Wr / 8
def rmul(X, G8):
    P = X @ G8; assert not np.any(P % 8); return P // 8
# Z-basis of L (echelon), coordinates, theta mod 3
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
random.seed(5); B = hnf(U[random.sample(range(len(U)), 400)])
assert np.prod([int(B[i, i]) for i in range(24)]) == 2 ** 36
def coords(X):
    X = X.copy(); c = np.zeros((len(X), 24), dtype=np.int64)
    for i in range(24):
        assert not np.any(X[:, i] % B[i, i]); c[:, i] = X[:, i] // B[i, i]; X -= np.outer(c[:, i], B[i])
    return c
def rref3(A):
    A = A.copy() % 3; r = 0; piv = []
    for c in range(A.shape[1]):
        p = next((i for i in range(r, A.shape[0]) if A[i, c] % 3), None)
        if p is None: continue
        A[[r, p]] = A[[p, r]]; A[r] = (A[r] * (1 if A[r, c] == 1 else 2)) % 3
        for i in range(A.shape[0]):
            if i != r and A[i, c]: A[i] = (A[i] - A[i, c] * A[r]) % 3
        piv.append(c); r += 1
        if r == A.shape[0]: break
    return A[:r], piv
def nullspace3(A):                                     # right null space of A (mod 3), as columns
    R, piv = rref3(A); free = [c for c in range(A.shape[1]) if c not in piv]
    Q = np.zeros((A.shape[1], len(free)), dtype=np.int64)
    for k, f in enumerate(free):
        Q[f, k] = 1
        for i, pc in enumerate(piv): Q[pc, k] = (-R[i, f]) % 3
    return Q
Th = coords(rmul(B, Wr) - rmul(rmul(B, Wr), Wr))     # theta in L-coordinates
Q = nullspace3(Th)                                     # 24x12: class of c is c Q mod 3
cls = (coords(U) @ Q) % 3
pw = 3 ** np.arange(11, -1, -1)
T = np.zeros(3 ** 12, dtype=bool); T[cls @ pw] = True
assert T.sum() == 65520
Fh = np.rint(np.real(np.fft.fftn(T.reshape([3] * 12).astype(float)))).astype(np.int64).ravel()
print('Fourier values of T:', dict(zip(*[x.tolist() for x in np.unique(Fh, return_counts=True)])))
# greedy search for a 6-dim T-free subspace
Tg = T.reshape([3] * 12); digits = np.array(np.unravel_index(np.arange(3 ** 12), [3] * 12)).T
def shift(A, x):
    for ax in range(12):
        if x[ax]: A = np.roll(A, int(x[ax]), axis=ax)
    return A
rng = np.random.default_rng(1); basisC = None
while basisC is None:
    Fb = Tg.copy(); Cl = [np.zeros(12, dtype=np.int64)]; bs = []
    for d in range(6):
        cand = np.flatnonzero((~Fb).ravel()); seen = set(map(tuple, Cl)); best = (-1, None, None)
        for s in rng.choice(cand, size=min(len(cand), 60 if d < 5 else len(cand)), replace=False):
            x = digits[s]
            if tuple(x) in seen: continue
            Fn = Fb | shift(Fb, x) | shift(Fb, (2 * x) % 3); a = int((~Fn).sum())
            if a > best[0]: best = (a, x, Fn)
            if d == 5 and a >= 729: break
        if best[0] < 729: break
        Fb = best[2]; x = best[1]; Cl = Cl + [(c + x) % 3 for c in Cl] + [(c + 2 * x) % 3 for c in Cl]; bs.append(x)
    if len(bs) == 6: basisC = np.array(bs)
Cl = np.array(Cl); assert len(set(map(tuple, Cl))) == 729 and not T[Cl @ pw].any()
# F: x in M iff (class of x) in C iff x B^{-1} Q Z^T = 0 mod 3, Z = basis of C^perp
def inv3(A):
    n = A.shape[0]; Mx = np.concatenate([A % 3, np.eye(n, dtype=np.int64)], axis=1)
    for c in range(n):
        p = next(i for i in range(c, n) if Mx[i, c] % 3); Mx[[c, p]] = Mx[[p, c]]
        Mx[c] = (Mx[c] * (1 if Mx[c, c] == 1 else 2)) % 3
        for i in range(n):
            if i != c and Mx[i, c]: Mx[i] = (Mx[i] - Mx[i, c] * Mx[c]) % 3
    return Mx[:, n:]
Z = nullspace3(basisC).T                               # 6x12, C^perp
F = (inv3(B) @ Q @ Z.T) % 3                            # 24x6
Ft, _ = rref3(F.T.copy())
inM = np.all((U @ Ft.T) % 3 == 0, axis=1)
print('norm-32 vectors in M:', int(inM.sum()))
out = os.path.join(HERE, 'F3_discovered.txt')
open(out, 'w').write(''.join(' '.join(map(str, r)) + '\n' for r in Ft.tolist()))
print('wrote', out, '(compare with F3.txt; a different random C gives an equivalent M)')

