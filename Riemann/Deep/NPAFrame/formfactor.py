#!/usr/bin/env python3
"""npaframe: form factor E|p_k|^2 (k = 1..2N) and type weights of saved laws; usage: formfactor.py file.npz ..."""
import sys
import numpy as np
for fn in sys.argv[1:]:
    d = np.load(fn, allow_pickle=True); w = d['sw']/d['sw'].sum(); M = d['sm']; TH = d['sth']; N = M.shape[1]
    k = np.arange(1, 2*N+1)
    P = np.einsum('ci,cik->ck', M, np.exp(1j*TH[:, :, None]*k[None, None, :]))
    S = (w[:, None]*np.abs(P)**2).sum(0)
    nd = (M == 2).sum(1); nt = (M == 3).sum(1)
    print(f"{fn}: N={N} support={len(w)} simple={float(w @ ((M == 1).sum(1)/N)):.6f} mean doubles={float(w @ nd):.3f} triples={float(w @ nt):.3f}")
    print("  E|p_k|^2/k: " + " ".join(f"{k_}:{s/k_:.3f}" for k_, s in zip(k, S)))
