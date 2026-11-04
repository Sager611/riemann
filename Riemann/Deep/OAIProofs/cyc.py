"""Isotypic support of a cyclic module C[S_N] w for w = x (x) y (x) ... (one vector per layer),
by exact class sums  A(mu) = sum_{g of type mu} prod_layers <x, g x>  and
||e_gamma w||^2 = (f^gamma / N!) sum_mu chi^gamma(mu) A(mu).
S^gamma occurs in C[S_N] w  iff  sum_mu chi^gamma(mu) A(mu) != 0 (it is then > 0).
Independent of every construction in the audited papers: plain enumeration of S_N."""
import itertools
import numpy as np
import mn

def cycle_type(p):
    n = len(p); seen = [False]*n; ct = []
    for i in range(n):
        if not seen[i]:
            l = 0; j = i
            while not seen[j]:
                seen[j] = True; j = p[j]; l += 1
            ct.append(l)
    return tuple(sorted(ct, reverse=True))

def layer_arrays(terms, N):
    """terms: list of (word tuple of ints >=0 of length N, int coeff)"""
    base = 1 + max(max(w) for w, c in terms)
    W = np.array([w for w, c in terms], dtype=np.int64)
    C = np.array([c for w, c in terms], dtype=np.int64)
    pw = base ** np.arange(N, dtype=np.int64)
    dense = np.zeros(base**N, dtype=np.int64)
    np.add.at(dense, W @ pw, C)
    return W, C, pw, dense

def class_sums(layers, N, batch=4000):
    """layers: list of term lists. returns dict mu -> exact int A(mu)"""
    arrs = [layer_arrays(t, N) for t in layers]
    tmax = max(len(t) for t in layers)
    batch = max(20, min(batch, 6_000_000 // (tmax * N)))   # keeps the (T,B,N) int64 array under ~50 MB
    P = list(mn.partitions(N))
    idx = {mu: k for k, mu in enumerate(P)}
    A = np.zeros(len(P), dtype=np.int64)
    it = itertools.permutations(range(N))
    while True:
        chunk = list(itertools.islice(it, batch))
        if not chunk: break
        G = np.array(chunk, dtype=np.int64)            # (B, N): (g.x)(word) = x(word o g)
        vals = np.ones(len(chunk), dtype=np.int64)
        for (W, C, pw, dense) in arrs:
            PW = W[:, G]                                # (T, B, N)  word o g
            codes = PW @ pw                             # (T, B)
            vals *= (C[:, None] * dense[codes]).sum(axis=0)
        cls = np.array([idx[cycle_type(p)] for p in chunk], dtype=np.int64)
        np.add.at(A, cls, vals)
    return {mu: int(A[idx[mu]]) for mu in P}

def support(layers, N):
    A = class_sums(layers, N)
    P = list(mn.partitions(N))
    norms = {}
    for gam in P:
        s = sum(mn.chi(gam, mu) * a for mu, a in A.items() if a)
        norms[gam] = s
    return norms
