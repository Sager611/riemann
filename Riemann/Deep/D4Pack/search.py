#!/usr/bin/env python3
"""search.py -- look for a Voronoi cell of volume < 8 (unit balls, centre 0).

Variables: neighbour centres Y (unit frame).  Packing: |y_i| >= 2, |y_i - y_j| >= 2.
Objective: volume of {x : x.y_i <= |y_i|^2/2} (own polytope code, voro4.float_cell),
gradient d vol / d y_i = A_i (y_i - c_i)/|y_i| (A_i facet 3-volume, c_i facet centroid).
A far box (halfspaces at distance 4, not centres) keeps every iterate bounded; a final
volume below 8 is re-checked without the box and in exact arithmetic.
Projected descent: gradient step, then repair of the packing constraints, accept if lower.
usage: search.py seed nstarts maxsteps kind   (kind: d4 | d4hole | c600 | shell | all)
"""
import sys, time, itertools
import numpy as np
from voro4 import float_cell, d4_roots

seed, nstarts, maxsteps = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
kind = sys.argv[4] if len(sys.argv) > 4 else "all"
rng = np.random.default_rng(seed)
ROOTS = np.array(d4_roots(), float) / np.sqrt(2)
BOX = np.vstack([np.eye(4), -np.eye(4)])   # far cube |x_k| <= 4
BOXB = np.full(len(BOX), 4.0)

def cell(Y):
    A = np.vstack([Y, BOX]); b = np.concatenate([0.5 * (Y * Y).sum(1), BOXB])
    vol, fd = float_cell(A, b)
    if vol is None:
        return None, None
    G = np.zeros_like(Y)
    for i, (area, cen) in fd.items():
        if i < len(Y):
            G[i] = area * (Y[i] - cen) / np.linalg.norm(Y[i])
    return vol, G

def repair(Y, margin=1e-9, iters=400):
    Y = Y.copy()
    for _ in range(iters):
        r = np.linalg.norm(Y, axis=1)
        low = r < 2 + margin
        Y[low] *= ((2 + margin) / r[low])[:, None]
        D = Y[:, None, :] - Y[None, :, :]
        d = np.linalg.norm(D, axis=2) + np.eye(len(Y)) * 9
        bad = np.argwhere(np.triu(d < 2 + margin, 1))
        if len(bad) == 0 and not low.any():
            return Y
        for i, j in bad:
            v = Y[i] - Y[j]; n = np.linalg.norm(v)
            push = (2 + margin - n) / 2 + 1e-12
            Y[i] += push * v / n; Y[j] -= push * v / n
    # fallback: inflate (always restores validity)
    D = np.linalg.norm(Y[:, None] - Y[None], axis=2) + np.eye(len(Y)) * 9
    lam = max(1.0, (2 + margin) / D.min(), (2 + margin) / np.linalg.norm(Y, axis=1).min())
    return Y * lam

def valid(Y):
    D = np.linalg.norm(Y[:, None] - Y[None], axis=2) + np.eye(len(Y)) * 9
    return D.min() >= 2 - 1e-9 and np.linalg.norm(Y, axis=1).min() >= 2 - 1e-9

def descend(Y):
    Y = repair(Y)
    v, G = cell(Y)
    if v is None:
        return None, Y
    eta = 0.02
    for _ in range(maxsteps):
        g = G / max(1e-12, np.abs(G).max())
        Yn = repair(Y - eta * g)
        vn, Gn = cell(Yn)
        if vn is not None and vn < v - 1e-12:
            Y, v, G = Yn, vn, Gn
            eta = min(eta * 1.3, 0.2)
        else:
            eta *= 0.5
            if eta < 1e-6:
                break
    return v, Y

def unit(v):
    return v / np.linalg.norm(v, axis=-1, keepdims=True)

def start(k):
    if k == "d4":          # perturbed D4 shell
        U = unit(ROOTS + rng.normal(scale=rng.uniform(0.03, 0.25), size=ROOTS.shape))
        return U * rng.uniform(2.0, 2.15, size=(24, 1))
    if k == "d4hole":      # D4 with roots removed and centres inserted near deep holes / gaps
        keep = rng.permutation(24)[: 24 - rng.integers(1, 5)]
        U = list(ROOTS[keep])
        holes = np.vstack([np.eye(4), -np.eye(4), np.array(list(itertools.product((.5, -.5), repeat=4)))])
        for h in holes[rng.permutation(24)[: rng.integers(1, 6)]]:
            U.append(unit(h + rng.normal(scale=0.1, size=4)))
        U = np.array(U)
        return U * rng.uniform(2.0, 2.4, size=(len(U), 1))
    if k == "c600":        # subsets of the 600-cell vertices (binary icosahedral group)
        phi = (1 + 5 ** .5) / 2
        V = [np.eye(4)[i] * s for i in range(4) for s in (1, -1)]
        V += [np.array(s) / 2 for s in itertools.product((1, -1), repeat=4)]
        base = np.array([phi, 1, 1 / phi, 0]) / 2
        for p in itertools.permutations(range(4)):
            sgn = sum(1 for a in range(4) for b in range(a + 1, 4) if p[a] > p[b]) % 2
            if sgn: continue
            for s in itertools.product((1, -1), repeat=3):
                w = np.zeros(4); w[list(p)] = base * np.array(list(s) + [1])
                V.append(w)
        V = np.unique(np.round(np.array(V), 12), axis=0)
        n = rng.integers(22, 31)
        U = V[rng.permutation(len(V))[:n]] @ np.linalg.qr(rng.normal(size=(4, 4)))[0]
        return U * rng.uniform(2.0, 2.3, size=(n, 1))
    n = rng.integers(22, 31)  # dense random shell
    return unit(rng.normal(size=(n, 4))) * rng.uniform(2.0, 2.25, size=(n, 1))

kinds = ["d4", "d4hole", "c600", "shell"] if kind == "all" else [kind]
t0 = time.time(); best = {}
for s in range(nstarts):
    k = kinds[s % len(kinds)]
    v, Y = descend(start(k))
    if v is None:
        continue
    r = np.linalg.norm(Y, axis=1)
    ncont = int((r < 2 + 1e-6).sum()); nnear = int((r < 2 * np.sqrt(2)).sum())
    if k not in best or v < best[k][0]:
        best[k] = (v, ncont, nnear, len(Y))
    if v < 8.0:
        np.save(f"cand_{seed}_{s}.npy", Y)
        print(f"!! start {s} kind {k}: volume {v:.9f} < 8, saved", flush=True)
print(f"seed {seed}: {nstarts} starts, {time.time()-t0:.0f}s; best per kind (vol, contacts, within 2sqrt2, N):")
for k, val in best.items():
    print(f"   {k:7s} {val[0]:.6f}  {val[1]:2d} {val[2]:2d} {val[3]:2d}")

