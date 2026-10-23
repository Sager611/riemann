#!/usr/bin/env python3
"""Numerical check of Banks 2610.31402, Prop. mu: for every Q-subspace U of K (dim r) and every set T
of t embeddings, rank(u^T) >= t r / n_K; hence mu(Gamma) = n_K/s.  Fields K = Q(sqrt d, xi_l), k = Q,
S = {1} u (G \ G_F), F = Q(sqrt d).  numpy only, < 100 MB, < 60 s."""
import itertools, random
import numpy as np

def setup(d, l):
    sq = np.sqrt(complex(d))
    xi = np.exp(2j * np.pi / l)
    emb = [(e, a) for e in (1, -1) for a in range(1, l)]           # G = {+-1} x (Z/l)^*
    basis = [(b, j) for b in (0, 1) for j in range(l - 1)]          # sqrt(d)^b xi^j
    M = np.array([[(e * sq) ** b * xi ** (a * j) for (b, j) in basis] for (e, a) in emb])
    return emb, basis, M

def rank(X):
    return int(np.linalg.matrix_rank(X, tol=1e-7)) if X.size else 0

random.seed(7)
total_viol = 0
for (d, l) in [(13, 5), (-7, 5), (5, 7), (-3, 7), (2, 7)]:
    emb, basis, M = setup(d, l)
    nK = len(emb)
    S = [0] + [i for i, (e, a) in enumerate(emb) if e == -1]        # identity + coset G \ G_F
    s = len(S)
    assert rank(M) == nK
    viol = 0; trials = 0; worst = 9.0
    # (a) coordinate subspaces U_J (all subsets for l = 5, sampled for l = 7)
    Js = (list(itertools.chain.from_iterable(itertools.combinations(range(nK), r) for r in range(1, nK + 1)))
          if nK <= 8 else [tuple(random.sample(range(nK), random.randint(1, nK))) for _ in range(4000)])
    for J in Js:
        r = len(J); dS = rank(M[np.ix_(S, list(J))]); trials += 1
        if r * s > nK * dS + 1e-9: viol += 1
        if dS < s: worst = min(worst, (nK - r) / (s - dS))
    # (b) random sparse rational subspaces, random T (any embeddings) and T = S
    for _ in range(6000):
        r = random.randint(1, nK)
        C = np.array([[random.choice([0, 0, 0, 1, -1, 2]) for _ in range(nK)] for _ in range(r)], dtype=float)
        if rank(C) < r: continue
        U = M @ C.T                                                 # columns u_i^E (all embeddings)
        for T in (S, random.sample(range(nK), random.randint(1, nK))):
            t = len(T); dT = rank(U[T, :]); trials += 1
            if r * t > nK * dT + 1e-9: viol += 1
            if T is S and dT < s: worst = min(worst, (nK - r) / (s - dT))
    # (c) subfields E = K^{H'} for all subgroups H' (G abelian): rank of S-rows on E vs [E:Q]
    G = emb
    def mul(g, h): return (g[0] * h[0], (g[1] * h[1]) % l)
    elems = list(G)
    subgroups = set()
    for gens in itertools.chain.from_iterable(itertools.combinations(elems, k) for k in range(0, 3)):
        Hs = {(1, 1)}
        changed = True
        while changed:
            changed = False
            for x in list(Hs):
                for g in gens:
                    y = mul(x, g)
                    if y not in Hs: Hs.add(y); changed = True
        subgroups.add(frozenset(Hs))
    for Hs in subgroups:
        # trace vectors: for basis beta, (sum_{h in H'} (iota o h)(beta))_iota
        idx = {g: i for i, g in enumerate(emb)}
        V = np.array([[sum(M[idx[mul(g, h)], c] for h in Hs) for c in range(nK)] for g in emb])
        r = rank(V); dS = rank(V[S, :]); trials += 1
        if r * s > nK * dS + 1e-9: viol += 1
        if dS < s: worst = min(worst, (nK - r) / (s - dS))
    total_viol += viol
    print(f"K=Q(sqrt{d}, xi_{l}): n_K={nK} s={s} n_K/s={nK/s:.4f}  trials={trials}  violations={viol}"
          f"  min tested (n_K-r)/(s-d)={worst:.4f}  #subgroups={len(subgroups)}")
print("total violations:", total_viol)

