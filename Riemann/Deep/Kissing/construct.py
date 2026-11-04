#!/usr/bin/env python3
"""construct.py -- regenerate equator_lines_13104.txt from golay_basis.txt and F3.txt (NumPy, ~1 minute).
Enumerates the norm-48 vectors of the Leech lattice L (shapes 2^12 0^12, 4 2^8 0^15, 3^3 1^21, 5 1^23), keeps one
vector per line (first nonzero entry positive) and those in M = {x in L : F x = 0 mod 3}; compares with the shipped
file."""
import numpy as np, os, itertools
HERE = os.path.dirname(os.path.abspath(__file__))
rows = [l.strip() for l in open(os.path.join(HERE, 'golay_basis.txt')) if l.strip()]
words = [0]
for r in rows:
    b = sum(1 << i for i, c in enumerate(r) if c == '1'); words += [w ^ b for w in words]
F = np.array([list(map(int, l.split())) for l in open(os.path.join(HERE, 'F3.txt')) if l.strip()]).T  # 24x6
out = []
def keep(X):
    X = np.asarray(X, dtype=np.int64)
    if len(X) == 0: return
    first = X[np.arange(len(X)), np.argmax(X != 0, axis=1)]
    X = X[(first > 0) & np.all((X @ F) % 3 == 0, axis=1) & (X.sum(axis=1) % 8 == np.where(X[:, 0] % 2 == 0, 0, 4))]
    out.extend(map(tuple, X.tolist()))
signs8 = np.array([[-1 if m >> j & 1 else 1 for j in range(8)] for m in range(256)])
signs12 = np.array([[-1 if m >> j & 1 else 1 for j in range(12)] for m in range(4096)])
trip = np.array(list(itertools.combinations(range(24), 3)))
for w in words:
    k = bin(w).count('1'); pos = [i for i in range(24) if w >> i & 1]
    if k == 12:
        X = np.zeros((4096, 24), dtype=np.int64); X[:, pos] = 2 * signs12; keep(X)
    if k == 8:
        for a in range(24):
            if w >> a & 1: continue
            for s4 in (4, -4):
                X = np.zeros((256, 24), dtype=np.int64); X[:, pos] = 2 * signs8; X[:, a] = s4; keep(X)
    base = np.array([-1 if w >> i & 1 else 1 for i in range(24)], dtype=np.int64)
    X = np.tile(base, (24, 1)); X[np.arange(24), np.arange(24)] *= 5; keep(X)
    X = np.tile(base, (len(trip), 1))
    for t in range(3): X[np.arange(len(trip)), trip[:, t]] *= -3
    keep(X)
out = sorted(set(out))
print(len(out), 'lines')
shipped = sorted(tuple(map(int, l.split())) for l in open(os.path.join(HERE, 'equator_lines_13104.txt')) if l.strip())
print('identical to shipped file:', out == shipped)

